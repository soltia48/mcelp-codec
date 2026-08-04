//! The encoder's noise-analysis front end.
//!
//! This is the glue between the spectral analysis and the voice-activity
//! detector: the magnitude spectrum is brought back to unnormalised units,
//! folded into critical bands, and compared against the two long-term noise
//! references to produce the handful of features the detector votes on.
//!
//! [`Suppressor`] then chains that decision into the noise floors and the
//! per-band weights, and [`Frontend`] is the whole of it behind one call.

use crate::analysis::BINS;
use crate::bands::{self, BANDS};
use crate::fixed::{acc, hi, sat, shift};
use crate::vad::{self, Features, Frame};

/// Frames of silence the analysis sits out at start-up before it trusts its
/// own long-term estimates.
pub const WARMUP: i16 = 11;

/// Renormalise the magnitude spectrum out of the analysis headroom.
pub fn unnormalise(magnitude: &[i16; BINS], headroom: i16) -> [i16; BINS] {
    let mut out = [0i16; BINS];
    for (o, &m) in out.iter_mut().zip(magnitude.iter()) {
        *o = hi(sat(shift(acc((m as i64) << 16), -(headroom as i32))));
    }
    out
}

/// Build the detector's features from one frame's spectrum.
///
/// `fast` and `slow` are the two long-term band-energy estimates; the frame is
/// described by how far the current band levels sit above each of them, both
/// band by band and as a whole.
pub fn features(
    spectrum: &[i16; BINS],
    fast: &[i64; BANDS],
    slow: &[i64; BANDS],
    voicing: i16,
) -> Frame {
    let levels = bands::levels(spectrum);
    let summary = bands::summary(&levels);

    let mut frame = Frame {
        summary,
        voicing,
        ..Frame::default()
    };
    for (i, energies) in [fast, slow].into_iter().enumerate() {
        let tracked = bands::tracked_levels(energies);
        let excess = bands::excess(&levels, &tracked);
        frame.reference[i] = Features {
            tracked: bands::summary(&tracked),
            excess_sum: vad::excess_sum(&excess),
            difference: bands::summary_difference(summary, bands::summary(&tracked)),
            variance: bands::snr_statistics(&excess).1,
        };
    }
    frame
}

use crate::noise::NoiseEstimate;
use crate::shaping::{self, HOP, Levels, Overlap, SPAN, Scale, Smoother};
use crate::vad::{Decision, Vad};
use crate::weights::{self, Depths, Mix, Shaping};

/// Everything the noise-suppression front end carries between frames.
pub struct Suppressor {
    /// Frames still to sit out before the analysis trusts itself.
    warmup: i16,
    detector: Vad,
    fast: NoiseEstimate,
    slow: NoiseEstimate,
    levels: Levels,
    depths: Depths,
    shaping: Shaping,
    scale: Scale,
    smoother: Smoother,
    overlap: Overlap,
}

impl Default for Suppressor {
    fn default() -> Self {
        Suppressor {
            warmup: WARMUP,
            detector: Vad::default(),
            fast: NoiseEstimate::fast(),
            slow: NoiseEstimate::slow(),
            levels: Levels::default(),
            depths: Depths::default(),
            shaping: Shaping::default(),
            scale: Scale::default(),
            smoother: Smoother::default(),
            overlap: Overlap::default(),
        }
    }
}

impl Suppressor {
    /// Run one half-frame through the whole path and return the shaping signal.
    ///
    /// `spectrum` is the complex transform of the analysis window and
    /// `magnitude` its magnitude; both come straight from the front end.
    pub fn run(
        &mut self,
        spectrum: &mut crate::analysis::Spectrum,
        magnitude: &[i16; BINS],
        headroom: i16,
        voicing: i16,
    ) -> [i16; HOP] {
        // Voice activity, unless the analysis is still warming up.
        let decision = if self.warmup >= 0 {
            self.warmup -= 1;
            Decision {
                active: false,
                fast: false,
                score: -4,
            }
        } else {
            let scaled = unnormalise(magnitude, headroom);
            let frame = features(&scaled, &self.fast.energy, &self.slow.energy, voicing);
            self.detector.run(&frame)
        };

        // Long-term noise, then everything that is derived from it.
        self.fast.update(magnitude, headroom, decision.fast);
        self.slow.update(magnitude, headroom, decision.active);

        let level = weights::overall_level(&self.fast.energy);
        let (band_level, sum) =
            self.levels
                .run(magnitude, &self.slow.energy, decision.score, headroom);
        let (mut depth, weight) = weights::depths(sum);
        let (mut mix, _) = Mix::run(&band_level, sum, level, decision.score, decision.active);
        self.depths.refine(
            &mut depth,
            &mut mix,
            &self.fast.energy,
            &band_level,
            &weight,
            decision.active,
        );
        let (band_weight, mut band_scale) =
            self.shaping
                .finish(&mut mix, &band_level, level, sum, decision.active);
        self.scale.smooth(&mut band_scale, decision.active);

        let floor = shaping::noise_floor(
            magnitude,
            &self.fast.energy,
            &band_scale,
            &band_weight,
            headroom,
            level,
            sum,
        );

        // Shape the spectrum and take it back to the time domain.
        shaping::suppress(spectrum, magnitude, &floor, headroom);
        let mut span = shaping::inverse_transform(&spectrum.re, &spectrum.im);
        self.smoother.run(&mut span);
        self.overlap.emit(&span, headroom)
    }
}

/// Compile-time check that the shaping span is the size this module expects.
const _: () = assert!(SPAN == 256);

use crate::analysis::{self, Analysis, Direction, WINDOW};
use crate::preprocess::{self, HALF, Highpass};

/// The encoder's analysis front end: mu-law samples in, shaping signal out.
#[derive(Default)]
pub struct Frontend {
    input: Highpass,
    analysis: Analysis,
    suppressor: Suppressor,
}

impl Frontend {
    /// Expand and condition one frame's worth of mu-law samples.
    ///
    /// The high-pass runs across the whole frame; everything after it works on
    /// half-frames.
    pub fn condition(&mut self, frame: &[u8; preprocess::FRAME]) -> [i16; preprocess::FRAME] {
        let mut linear = [0i16; preprocess::FRAME];
        for (l, &b) in linear.iter_mut().zip(frame.iter()) {
            *l = preprocess::ulaw_to_linear(b);
        }
        self.input.run(&mut linear);
        linear
    }

    /// Run one conditioned half-frame through the analysis and the suppressor.
    pub fn process(&mut self, block: &[i16; HALF]) -> [i16; HOP] {
        let (mut window, headroom) = self.analysis.window(block);
        let voicing = analysis::voicing(&window);
        self.analysis.preemphasise(&mut window, headroom);

        let mut spectrum = analysis::deinterleave(&window);
        analysis::fft(&mut spectrum);
        analysis::unpack_real(&mut spectrum, Direction::Forward);
        let magnitude = analysis::magnitude(&spectrum);

        self.suppressor
            .run(&mut spectrum, &magnitude, headroom, voicing)
    }
}

/// Compile-time check that the analysis window and the shaping span agree.
const _: () = assert!(WINDOW == SPAN);
