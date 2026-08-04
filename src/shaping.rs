//! The spectrally shaped signal: applying the suppression and getting back to
//! the time domain.
//!
//! Overlap handling is what the
//! module is named after.  The shaping path transforms a whole 256-sample
//! analysis span at a time but the encoder consumes 160 samples per half-frame,
//! so consecutive spans have to be stitched together.  Rather than a full
//! overlap-add the reference cross-fades a 40-sample seam and takes the rest of
//! each span straight, which is enough because the two spans differ only where
//! the analysis windows do.
//!
//! The rest of the module is the shaping itself: the inverse transform and its
//! leaky integrator, the per-band noise scale and its smoother
//!, and the gain the
//! whole span is finally taken through.

use crate::analysis::{BINS, Direction, FFT_POINTS, Spectrum, fft, negate, ratio, unpack_real};
use crate::bands::BANDS;
use crate::fixed::{acc, exp, hi, low, mulr, norm, round, sat, scale32, shift};
use crate::tables::{SCALE_SMOOTH_ACTIVE as SMOOTH_ACTIVE, SCALE_SMOOTH_IDLE as SMOOTH_IDLE};

/// Samples the transform produces per frame.
pub const SPAN: usize = 256;
/// Samples emitted per half-frame.
pub const HOP: usize = 160;
/// Length of the cross-faded seam.
pub const SEAM: usize = 40;
/// Where in the span the emitted block starts.
const START: usize = 96;
/// Step the cross-fade weights move by, and the total they always sum to.
const WEIGHT_STEP: i16 = 512;
const WEIGHT_TOTAL: i16 = 20480;
/// Gain that turns the weight sum back into unity.
const WEIGHT_GAIN: i16 = 26214;

/// Carry-over between consecutive spans.
#[derive(Clone, Copy)]
pub struct Overlap {
    /// The seam of the previous span, kept in that span's scale.
    tail: [i16; SEAM],
    /// The headroom the tail is scaled by; `None` before the first span.
    scale: Option<i16>,
}

impl Default for Overlap {
    fn default() -> Self {
        Overlap {
            tail: [0; SEAM],
            scale: None,
        }
    }
}

impl Overlap {
    /// Stitch one span onto the previous one and emit a half-frame.
    ///
    /// `headroom` is the shift the analysis normalised this frame by; the
    /// result is handed back in unnormalised units.
    pub fn emit(&mut self, span: &[i16; SPAN], headroom: i16) -> [i16; HOP] {
        // Bring the previous seam into this span's scale before mixing.
        if let Some(previous) = self.scale {
            let step = headroom - previous;
            for v in self.tail.iter_mut() {
                *v = hi(norm((*v as i64) << 16, step));
            }
        }

        let mut out = [0i16; HOP];
        let seam = &span[START - SEAM..START];
        for i in 0..SEAM {
            let rising = WEIGHT_STEP * i as i16;
            let falling = WEIGHT_TOTAL - rising;
            let mut a = sat((seam[i] as i64) * (rising as i64) * 2);
            let b = sat((self.tail[i] as i64) * (falling as i64) * 2);
            a = acc(((hi(a) as i64) << 16) + b);
            let mixed = sat((WEIGHT_GAIN as i64) * (hi(a) as i64) * 2);
            out[i] = hi(sat(shift((hi(mixed) as i64) << 16, 1)));
        }
        out[SEAM..].copy_from_slice(&span[START..START + HOP - SEAM]);
        self.tail.copy_from_slice(&span[SPAN - SEAM..]);
        self.scale = Some(headroom);

        for v in out.iter_mut() {
            *v = hi(norm((*v as i64) << 16, -headroom));
        }
        out
    }
}

/// Coefficient of the leaky integrator that follows the inverse transform.
const LEAK: i16 = 26214;

/// First-order integrator run along the span.
///
/// The inverse transform leaves a differentiated signal behind, so the span is
/// run through `y[n] = 0.8 y[n-1] + x[n]` before it is stitched together.  The
/// recursion carries across spans, which is why it needs state of its own.
#[derive(Clone, Copy, Default)]
pub struct Smoother {
    state: i16,
}

impl Smoother {
    /// Filter one span in place.
    pub fn run(&mut self, span: &mut [i16; SPAN]) {
        for v in span.iter_mut() {
            let leaked = round(sat((LEAK as i64) * (self.state as i64) * 2));
            self.state = hi(sat(acc(leaked + ((*v as i64) << 16))));
            *v = self.state;
        }
    }
}

/// Halved unit gain the two halves of the transform are scaled by, and the
/// left shift that undoes it.
const HALF: i16 = 16384;
const REGAIN: i32 = 2;

/// Turn a shaped 129-bin spectrum back into 256 time samples.
///
/// The same machinery that packs a 256-point real transform into 129 bins runs
/// in reverse, and the 128-point complex transform that follows it is the
/// ordinary forward one with the imaginary half conjugated on the way in.  What
/// comes out is the real sequence interleaved across the two halves.
pub fn inverse_transform(re: &[i16; BINS], im: &[i16; BINS]) -> [i16; SPAN] {
    let mut s = Spectrum { re: *re, im: *im };
    unpack_real(&mut s, Direction::Inverse);
    for v in s.im[..FFT_POINTS].iter_mut() {
        *v = negate(*v);
    }
    fft(&mut s);

    let mut out = [0i16; SPAN];
    for i in 0..FFT_POINTS {
        out[2 * i] = rescale(s.re[i], HALF);
        out[2 * i + 1] = rescale(s.im[i], -HALF);
    }
    out
}

/// Halve a sample, keep only the high word, then shift it back up.
fn rescale(v: i16, gain: i16) -> i16 {
    let halved = sat((v as i64) * (gain as i64) * 2);
    hi(sat(shift((hi(halved) as i64) << 16, REGAIN)))
}

/// Apply the noise-suppression gain to the complex spectrum.
///
/// Each bin is scaled by how much of it the noise floor accounts for: a bin
/// that sits entirely under the floor keeps its full value, one well above it
/// is cut down to the floor's share, and a bin with no energy at all is
/// zeroed.  What comes out is an estimate of the noise alone, which is what the
/// inverse transform turns back into a time-domain shaping signal.
///
/// The reference also mirrors the result into bins 129..255 to complete the
/// Hermitian spectrum, but nothing downstream reads those, so that is left out.
pub fn suppress(
    spectrum: &mut Spectrum,
    magnitude: &[i16; BINS],
    floor: &[i16; BINS],
    headroom: i16,
) {
    for i in 0..BINS {
        let gain = if magnitude[i] <= 0 {
            0
        } else if magnitude[i] <= floor[i] {
            32767
        } else {
            ratio(floor[i], magnitude[i])
        };
        spectrum.re[i] = scaled(spectrum.re[i], gain, headroom);
        spectrum.im[i] = scaled(spectrum.im[i], gain, headroom);
    }
}

/// One bin: scale, undo the analysis headroom, round.
fn scaled(v: i16, gain: i16, headroom: i16) -> i16 {
    let product = sat((v as i64) * (gain as i64) * 2);
    hi(sat(acc(norm(product, headroom) + 32768)))
}

/// Smallest share of a bin that survives, chosen from how loud the frame is.
const FLOOR_GAINS: [i16; 3] = [5194, 7337, 10361];
/// Level above which the mildest gain is used unconditionally.
const LOUD_LEVEL: i16 = 15360;
/// Thresholds on the secondary level measure.
const QUIET: i16 = 6144;
const MODERATE: i16 = 12288;

/// Build the per-bin noise floor the suppression gain works from.
///
/// Two estimates of the floor are computed for every bin and the larger wins:
/// spectral subtraction, which takes what is left after the band's tracked
/// noise level is removed and weights it per band, and a plain fraction of the
/// bin itself, which stops the floor from ever collapsing to nothing.
pub fn noise_floor(
    magnitude: &[i16; BINS],
    energy: &[i64; BANDS],
    scale: &[i16; BANDS],
    weight: &[i16; BANDS],
    headroom: i16,
    level: i16,
    secondary: i16,
) -> [i16; BINS] {
    let gain = if level >= LOUD_LEVEL {
        FLOOR_GAINS[2]
    } else if secondary < QUIET {
        FLOOR_GAINS[0]
    } else if secondary < MODERATE {
        FLOOR_GAINS[1]
    } else {
        FLOOR_GAINS[2]
    };

    let mut out = [0i16; BINS];
    let mut bin = 0usize;
    for band in 0..BANDS {
        let first = crate::tables::BAND_EDGES[2 * band] as usize;
        let last = crate::tables::BAND_EDGES[2 * band + 1] as usize;

        // The band's tracked noise, brought into the frame's scale.
        let raised = norm(energy[band], headroom);
        let noise = hi(shift(mulr(scale[band], hi(raised)), -4)) as i64;

        for _ in 0..=(last - first) {
            let m = magnitude[bin] as i64;
            let subtracted = norm(mulr(weight[band], hi(acc((m - noise) << 16))), -headroom);
            let scaled = norm(m << 16, -headroom);
            let fraction = sat((gain as i64) * (hi(scaled) as i64) * 2);
            out[bin] = hi(acc((hi(subtracted) as i64) << 16).max(acc((hi(fraction) as i64) << 16)));
            bin += 1;
        }
    }
    out
}

/// Per-band smoothing coefficients: during speech the scale is held much
/// closer to its previous value than during silence.
/// Total weight the two coefficients of a pair add up to.
const SCALE_UNITY: i16 = 16384;

/// One-pole smoother applied to the per-band noise scale.
///
/// The scale a frame computes is mixed with the one the last frame settled on,
/// so a single odd frame cannot swing the suppression.
#[derive(Clone, Copy, Default)]
pub struct Scale {
    previous: [i16; BANDS],
}

impl Scale {
    /// Smooth this frame's scale in place.
    pub fn smooth(&mut self, scale: &mut [i16; BANDS], active: bool) {
        let table = BANDS * !active as usize;
        for (band, (current, previous)) in
            scale.iter_mut().zip(self.previous.iter_mut()).enumerate()
        {
            let c = crate::tables::BAND_SMOOTH[table + band];
            let mut b = sat((c as i64) * (*current as i64) * 2);
            b = sat(acc(b + (*previous as i64) * ((SCALE_UNITY - c) as i64) * 2));
            let v = hi(sat(acc(sat(shift(b, 1)) + 32768)));
            *current = v;
            *previous = v;
        }
    }
}

/// How much of the frame is assumed to be noise, indexed by the voice-activity
/// score offset by four: a confident speech frame contributes nothing.
/// Per-band rate at which the running blend state follows that assumption.
/// Smoothing coefficients along the band axis; the idle set moves faster.
use crate::bands::DB_PER_OCTAVE;
/// Ceiling on a band level, and the level at which it is reached.
const LEVEL_CLIP: i16 = 7936;
const LEVEL_MAX: i16 = 31744;
/// Weight the band levels are summed with to give the frame's second measure.
const LEVEL_WEIGHT: i16 = 1638;

/// Per-band level of the mean magnitude against the tracked noise.
///
/// The core of it is a ratio: how far each band's mean magnitude sits above a
/// blend of the slow noise estimate and the magnitude itself, in decibels.  How
/// much of the blend comes from which side is a running state that follows the
/// voice-activity score, so during speech the comparison is made against the
/// noise estimate alone and during silence it slides towards the frame itself,
/// which drives the levels towards zero.
#[derive(Clone, Copy, Default)]
pub struct Levels {
    /// Blend state, carried across bands and across frames.
    state: i64,
    /// Last frame's smoothed level per band.
    previous: [i16; BANDS],
}

impl Levels {
    /// Returns the twenty band levels and their weighted sum.
    pub fn run(
        &mut self,
        magnitude: &[i16; BINS],
        slow: &[i64; BANDS],
        score: i16,
        headroom: i16,
    ) -> ([i16; BANDS], i16) {
        let confidence =
            crate::tables::CONFIDENCE[(crate::tables::CONFIDENCE_ZERO as i16 + score) as usize];
        let mut level = [0i16; BANDS];

        for band in 0..BANDS {
            if slow[band] == 0 {
                level[band] = 32767;
                continue;
            }
            let first = crate::tables::BAND_EDGES[2 * band] as usize;
            let last = crate::tables::BAND_EDGES[2 * band + 1] as usize;
            let mut sum = 0i64;
            for &m in magnitude[first..=last].iter() {
                sum = sat(acc(sum + ((m as i64) << 16)));
            }
            let mean = scale32(sum, crate::tables::BAND_INV_WIDTH_HALVED[band]);

            // Slide the blend state towards this frame's confidence.
            let rate = crate::tables::STATE_RATE[band];
            let mut s = sat(shift(scale32(self.state, rate), -5));
            s = sat(acc(s + (confidence as i64) * ((32767 - rate) as i64) * 2));
            self.state = sat(shift(s, -5));
            let mix = hi(s);

            let blend = sat(acc(
                scale32(slow[band], 32767 - mix) + norm(scale32(mean, mix), 5 - headroom)
            ));

            // Line both sides up, divide, and take the logarithm.
            let blend_exp = exp(blend);
            let mean_exp = exp(mean) - 1;
            let shift_out = headroom as i64 - blend_exp as i64 + 16 + mean_exp as i64;
            let ratio = divide32(norm32(mean, mean_exp), norm32(blend, blend_exp));

            let log = crate::bands::log2(ratio);
            let combined = acc(((hi(log) as i64 - shift_out) << 16) + (log & 0xffff));
            let db = scale32(acc(shift(combined, -1) - (5 << 16)), DB_PER_OCTAVE);
            // The clip test looks at the whole accumulator, not just its high
            // word, so a level a fraction above the ceiling still clips.
            level[band] = if acc(sat(shift(db, 12)) - ((LEVEL_CLIP as i64) << 16)) > 0 {
                LEVEL_MAX
            } else {
                hi(sat(shift(db, 14)).max(0))
            };
        }

        self.smooth(&mut level, score);
        let mut total = 0i64;
        for &l in level.iter() {
            total = sat(acc(total + (l as i64) * (LEVEL_WEIGHT as i64) * 2));
        }
        (level, hi(round(total)))
    }

    /// Smooth the levels forwards along the band axis and then backwards in
    /// time, band by band.
    fn smooth(&mut self, level: &mut [i16; BANDS], score: i16) {
        let table = if score >= 0 {
            SMOOTH_ACTIVE
        } else {
            SMOOTH_IDLE
        };
        for band in 0..BANDS - 1 {
            let c = crate::tables::SCALE_SMOOTH[table + band];
            let mut b = sat((level[band] as i64) * ((32767 - c) as i64) * 2);
            b = round(sat(acc(b + (level[band + 1] as i64) * (c as i64) * 2)));
            level[band + 1] = hi(b);
        }
        for band in (0..BANDS).rev() {
            let c = crate::tables::SCALE_SMOOTH[table + band - 1];
            let mut b = sat((self.previous[band] as i64) * ((32767 - c) as i64) * 2);
            b = round(sat(acc(b + (level[band] as i64) * (c as i64) * 2)));
            level[band] = hi(b);
            self.previous[band] = hi(b);
        }
    }
}

/// Normalise a 32-bit value by a shift.
fn norm32(v: i64, amount: i32) -> i64 {
    sat(shift(acc(v), amount))
}

/// `numerator / denominator` for two normalised 32-bit values.
///
/// A fifteen-step restoring division gives a first reciprocal, one Newton step
/// refines it against the full 32-bit denominator, and the numerator is then
/// multiplied by the refined value.
fn divide32(numerator: i64, denominator: i64) -> i64 {
    let seed = low(crate::fixed::restoring_divide(
        16383i64 << 16,
        hi(denominator),
        15,
    ));
    let product = scale32(denominator, seed);
    let error = sat(acc(sat(shift(16384i64 << 16, 8)) - product));
    let refined = sat(shift(scale32(error, seed), 0));
    mul32(numerator, refined)
}

/// Full 32x32 product, keeping the upper half and shifted back up by two.
fn mul32(x: i64, y: i64) -> i64 {
    let xl = ((x as u32 as u16) >> 1) as i64;
    let yl = ((y as u32 as u16) >> 1) as i64;
    let mut b = shift(sat((hi(x) as i64) * yl * 2), -16);
    b = sat(acc(b + shift(sat((hi(y) as i64) * xl * 2), -16)));
    let mut a = sat(shift(b, 1));
    a = sat(acc(a + (hi(x) as i64) * (hi(y) as i64) * 2));
    sat(shift(a, 2))
}
