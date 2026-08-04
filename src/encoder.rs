//! The encoder's frame loop.
//!
//! One frame of 320 mu-law samples goes in and one set of transport fields
//! comes out.  The frame is dealt with in two halves: each is expanded,
//! high-passed, put through the noise suppression and the shaping filters, and
//! what comes out of that is what the rest of the encoder analyses.
//!
//! The line spectrum is taken once per frame over a window that reaches back
//! into the previous one, quantised, and its two fields written; the remaining
//! twelve fields come from the two half-frames' subframe searches.

use crate::bitstream::{self, FIELDS, Params};
use crate::codebook;
use crate::convolve::filter;
use crate::excitation::{cap_gain, carry_gain, excite, interpolation_gains, update_memories};
use crate::fixed::acc;
use crate::frontend::Frontend;
use crate::gain::{GainState, align, search as gain_search};
use crate::lpc::{
    ANALYSIS, ORDER, autocorrelate, lag_window, levinson, line_spectrum, split_polynomials,
};
use crate::lsf_weight::{quantise, weights};
use crate::lsp::{self, LsfIndices, LsfState, Quantiser};
use crate::mode::Mode;
use crate::pitch;
use crate::pitch_search::{
    SEARCH_SPAN, adaptive_gain, best_closed_lag, cross, first as first_window, interpolate,
    lag_correlations, measures, open_loop_lag, periodic, prescale, search_range, track_peak,
};
use crate::postfilter::inverse_filter;
use crate::preprocess::{FRAME, HALF};
use crate::pulses;
use crate::shaping::HOP;
use crate::synth::synthesis;
use crate::weight_lpc::{Sharpness, expand, midpoint};
use crate::weighting::{Highpass, Smooth, Tilt};

/// Samples of the previous frame the analysis window reaches back for.
const REACH: usize = ANALYSIS - 2 * HALF;

/// Shaping samples the subframe filters see: eighty carried over from the
/// previous half-frame followed by this one's.
const HISTORY: usize = CARRY + HOP;
/// Samples of the previous half-frame kept in front of it.
const CARRY: usize = 80;
/// Where the first subframe's filtered window starts inside that.
const WINDOW: usize = 40;
/// Samples a subframe covers.
use crate::SUBFRAME as SUB;

/// The line spectrum the interpolation starts from: evenly spaced frequencies,
/// which is what a flat spectrum comes to.
const FLAT: [i16; ORDER] = [
    31441, 27568, 21460, 13612, 4663, -4663, -13612, -21460, -27568, -31441,
];

/// The three filters the shaping signal goes through.
#[derive(Default)]
struct Shaping {
    highpass: Highpass,
    smooth: Smooth,
    tilt: Tilt,
}

impl Shaping {
    fn run(&mut self, block: &mut [i16]) {
        self.highpass.run(block);
        self.smooth.run(block);
        self.tilt.run(block);
    }
}

/// Everything the encoder carries between frames.
pub struct Encoder {
    frontend: Frontend,
    shaping: Shaping,
    /// The tail of the previous frame's shaping signal, which the analysis
    /// window starts on.
    tail: [i16; REACH],
    /// The quantiser's predictor memory, kept in step with what the decoder
    /// will reconstruct.
    lsf: LsfState,
    /// Set while the line spectrum may be interpolated across the two
    /// half-frames.
    interpolating: bool,
    /// The previous frame's line spectrum, which the interpolation halves
    /// towards.
    previous: [i16; ORDER],
    /// What the first half-frame is weighted against: the two frames' midpoint
    /// while interpolation runs, or this frame's own line spectrum when it has
    /// just been allowed again.
    interpolated: [i16; ORDER],
    /// The line spectrum the last half-frame was coded with, which the next
    /// one's first subframe interpolates from.
    last_half: [i16; ORDER],
    /// The two subframe line spectra of each half-frame, in the cosine domain
    /// and as frequencies.
    subframe_lsp: [[[i16; ORDER]; 2]; 2],
    subframe_lsf: [[[i16; ORDER]; 2]; 2],
    /// The first two reflection coefficients, which the weighting filter's
    /// shape is decided from.
    reflection: [i16; 2],
    /// State of the weighting filter's shape decision.
    sharpness: Sharpness,
    /// Its numerator and denominator factors, two subframes per half-frame.
    factors: [([i16; 2], [i16; 2]); 2],
    /// The weighting filter itself: numerator and denominator for each of the
    /// two subframes of each half-frame.
    weighting: [[([i16; ORDER + 1], [i16; ORDER + 1]); 2]; 2],
    /// The shaping signal with the previous half-frame's tail in front of it.
    history: [i16; HISTORY],
    /// Memory of the weighting filter's denominator, carried between
    /// subframes.
    weighting_memory: [i16; ORDER],
    /// The weighted signal the open-loop search runs over: history followed by
    /// this half-frame's three passes.
    search: [i16; SEARCH_SPAN],
    /// Impulse response of the weighted synthesis filter, per half-frame.
    impulse: [[[i16; SUB]; 2]; 2],
    /// The prediction filter the decoder will have, which the encoder's own
    /// synthesis has to match.
    decoded_lpc: [i16; ORDER + 1],
    /// The same ladder as the weighting filter's, but climbing towards the
    /// line spectrum the decoder will rebuild rather than the analysis one.
    decoded_previous: [i16; ORDER],
    decoded_last_half: [i16; ORDER],
    synthesis_lpc: [[[i16; ORDER + 1]; 2]; 2],
    /// The open-loop lag and the closed-loop bracket it gives, per half-frame.
    open_loop: [i16; 2],
    bracket: [(i16, i16); 2],
    /// The prediction filter the first half-frame's second subframe is
    /// weighted with: the analysis one when this frame's line spectrum is
    /// taken whole, the interpolated one otherwise.
    second_lpc: [i16; ORDER + 1],

    /// Past excitation followed by room for this half-frame's two subframes.
    /// The closed loop reads its own history back out of it.
    excitation: [i16; EXCITATION],
    /// Memory of the filter that takes the residual back to speech: the last
    /// ten samples of the difference between the input and what the decoder
    /// will reconstruct.
    error_memory: [i16; ORDER],
    /// Memory of the weighting denominator inside the subframe loop.
    residual_memory: [i16; ORDER],
    /// Memory of the synthesis filter that reconstructs the speech.
    speech_memory: [i16; ORDER],
    /// Each subframe's search target, kept for comparison.
    subframe_target: [[[i16; SUB]; 2]; 2],
    /// The four most recent open-loop correlation peaks, which decide whether
    /// a subframe counts as periodic.
    peaks: [i64; PEAKS],
    /// The lag the second subframe of a half-frame codes against.
    first_lag: i16,
    /// The per-position correlations against the adaptive codebook, which the
    /// fixed codebook search refreshes a mode at a time.
    adaptive_correlations: [i16; SUB],
    /// The two gains a subframe's periodic extension may use, and the clamped
    /// adaptive gain the next subframe tests against.
    extension: (i16, i16),
    carried_gain: i16,
    /// The lag the last half-frame finished on, which decides whether the
    /// pitch has carried over, and the one this half-frame is building.
    previous_lag: i16,
    half_lag: i16,
    /// How each subframe is coded, and how its line spectrum is quantised.
    mode: Mode,
    quantiser: Quantiser,
    /// The gain history the code gain is predicted against, kept in step with
    /// what the decoder will have.
    gain_state: GainState,
    /// The prediction filter each subframe's weighting was built from, which
    /// the quantiser mode is also decided from.
    weighting_lpc: [[[i16; ORDER + 1]; 2]; 2],
    /// The fields the subframe searches fill in.
    field: [i16; FIELDS],
}

/// Open-loop peaks the periodicity test looks back over.
const PEAKS: usize = 4;
/// Lags the closed-loop search correlates over beyond the bracket itself: four
/// below and five above, for the fractional interpolation to lean on.
const LAGS_OVER: usize = 9;
/// Where the bracket starts inside that.
const LAG_MARGIN: usize = 4;
/// Impulse response taps the closed-loop convolution keeps.
const REACH_TAPS: usize = 40;
/// Lag above which the first subframe stops sending a fraction.
const COARSE_LAG: i16 = 84;

/// Excitation the closed loop keeps: enough history for the longest lag,
/// followed by the half-frame being coded.
const EXCITATION: usize = HISTORY_EXC + HOP;
/// Where this half-frame's excitation starts inside it.
const HISTORY_EXC: usize = 154;

impl Default for Encoder {
    fn default() -> Self {
        Encoder {
            frontend: Frontend::default(),
            shaping: Shaping::default(),
            tail: [0; REACH],
            lsf: LsfState::default(),
            interpolating: false,
            previous: FLAT,
            interpolated: [0; ORDER],
            last_half: FLAT,
            subframe_lsp: [[[0; ORDER]; 2]; 2],
            subframe_lsf: [[[0; ORDER]; 2]; 2],
            reflection: [0; 2],
            sharpness: Sharpness::default(),
            factors: [([0; 2], [0; 2]); 2],
            weighting: [[([0; ORDER + 1], [0; ORDER + 1]); 2]; 2],
            history: [0; HISTORY],
            weighting_memory: [0; ORDER],
            search: [0; SEARCH_SPAN],
            impulse: [[[0; SUB]; 2]; 2],
            decoded_lpc: [0; ORDER + 1],
            decoded_previous: FLAT,
            decoded_last_half: FLAT,
            synthesis_lpc: [[[0; ORDER + 1]; 2]; 2],
            open_loop: [0; 2],
            bracket: [(0, 0); 2],
            second_lpc: [0; ORDER + 1],
            excitation: [0; EXCITATION],
            error_memory: [0; ORDER],
            residual_memory: [0; ORDER],
            speech_memory: [0; ORDER],
            subframe_target: [[[0; SUB]; 2]; 2],
            peaks: [0; PEAKS],
            first_lag: 0,
            adaptive_correlations: [0; SUB],
            extension: (0, 0),
            carried_gain: 0,
            previous_lag: 0,
            half_lag: 0,
            mode: Mode::default(),
            quantiser: Quantiser::default(),
            gain_state: GainState::default(),
            weighting_lpc: [[[0; ORDER + 1]; 2]; 2],
            field: [0; FIELDS],
        }
    }
}

impl Encoder {
    /// An encoder in its reset state.
    pub fn new() -> Self {
        Self::default()
    }

    // The accessors below expose intermediate state so that each stage can be
    // replayed against the reference.  They are not part of the
    // library's surface and may change with the internals.

    /// Each subframe's synthesis filter.
    #[doc(hidden)]
    pub fn synthesis_lpc(&self) -> &[[[i16; ORDER + 1]; 2]; 2] {
        &self.synthesis_lpc
    }

    /// The impulse response the searches convolve with.
    #[doc(hidden)]
    pub fn impulse(&self) -> &[[[i16; SUB]; 2]; 2] {
        &self.impulse
    }

    /// The weighting filters themselves.
    #[doc(hidden)]
    pub fn weighting(&self) -> &[[([i16; ORDER + 1], [i16; ORDER + 1]); 2]; 2] {
        &self.weighting
    }

    /// Shape one frame and return the two half-frames' shaping signals.
    fn shape(&mut self, frame: &[u8; FRAME]) -> [[i16; HOP]; 2] {
        let linear = self.frontend.condition(frame);
        let mut out = [[0i16; HOP]; 2];
        for (half, signal) in out.iter_mut().enumerate() {
            let mut block = [0i16; HALF];
            block.copy_from_slice(&linear[half * HALF..(half + 1) * HALF]);
            *signal = self.frontend.process(&block);
            self.shaping.run(signal);
        }
        out
    }

    /// The line spectrum of one frame, taken over a window that reaches back
    /// into the previous one.
    fn line_spectrum(&mut self, shaped: &[[i16; HOP]; 2]) -> ([i16; ORDER], [i16; ORDER + 1]) {
        let mut window = [0i16; ANALYSIS];
        window[..REACH].copy_from_slice(&self.tail);
        for (half, signal) in shaped.iter().enumerate() {
            let at = REACH + half * HOP;
            window[at..at + HOP].copy_from_slice(signal);
        }
        self.tail.copy_from_slice(&shaped[1][HOP - REACH..]);

        let mut correlation = autocorrelate(&window);
        lag_window(&mut correlation);
        let (a, k) = levinson(&correlation);
        self.reflection = [k[0], k[1]];
        let (sum, difference) = split_polynomials(&a);
        (line_spectrum(&sum, &difference), a)
    }

    /// Each subframe's search target.
    #[doc(hidden)]
    pub fn subframe_target(&self) -> &[[[i16; SUB]; 2]; 2] {
        &self.subframe_target
    }

    /// One subframe of the closed loop.
    ///
    /// The target the two codebooks are searched against is the weighted
    /// signal with everything the previous subframes already accounted for
    /// taken out.  It is built by putting the shaping signal through the
    /// filter the decoder will have, running the residual back through that
    /// filter starting from the error the last subframe left, and weighting
    /// what comes out.
    fn subframe(&mut self, half: usize, sub: usize) {
        let offset = sub * SUB;
        let start = CARRY - WINDOW + offset;
        let synthesis_filter = self.synthesis_lpc[half][sub];

        let mut residual = [0i16; SUB];
        inverse_filter(
            &synthesis_filter,
            &self.history[start - ORDER..start + SUB],
            &mut residual,
        );

        // The residual is left in the excitation buffer, where the closed-loop
        // search is about to read it: a lag shorter than the subframe reaches
        // into the part of the buffer this subframe has not filled in yet.
        self.excitation[HISTORY_EXC + offset..HISTORY_EXC + offset + SUB]
            .copy_from_slice(&residual);

        let mut resynthesised = [0i16; SUB];
        let mut memory = self.error_memory;
        synthesis(
            &synthesis_filter,
            &residual,
            &mut resynthesised,
            &mut memory,
        );

        let (numerator, denominator) = self.weighting[half][sub];
        let mut window = [0i16; ORDER + SUB];
        window[..ORDER].copy_from_slice(&self.error_memory);
        window[ORDER..].copy_from_slice(&resynthesised);
        let mut weighted = [0i16; SUB];
        inverse_filter(&numerator, &window, &mut weighted);

        let mut target = [0i16; SUB];
        let mut memory = self.residual_memory;
        synthesis(&denominator, &weighted, &mut target, &mut memory);
        self.subframe_target[half][sub] = target;

        // The impulse response the searches convolve with is this subframe's
        // own: its weighting numerator through the synthesis filter and then
        // through the weighting denominator, both starting from rest.
        let mut fed = [0i16; SUB];
        fed[..ORDER + 1].copy_from_slice(&numerator);
        let mut once = [0i16; SUB];
        let mut rest = [0i16; ORDER];
        synthesis(&synthesis_filter, &fed, &mut once, &mut rest);
        let mut twice = [0i16; SUB];
        let mut rest = [0i16; ORDER];
        synthesis(&denominator, &once, &mut twice, &mut rest);
        self.impulse[half][sub] = twice;

        let lag = self.closed_loop_lag(half, sub, &target);
        let field = Params::subframe_base(half, sub);
        self.field[field + bitstream::LAG] = if sub == 0 {
            pitch::encode_absolute(&lag)
        } else {
            pitch::encode_relative(&lag, self.first_lag)
        };
        if sub == 0 {
            self.first_lag = lag.integer;
        }

        // The adaptive codebook's contribution, and how much of the target it
        // accounts for.
        let at = HISTORY_EXC + offset;
        pitch::predict(&mut self.excitation, at, &lag);
        let adaptive: [i16; SUB] = self.excitation[at..at + SUB].try_into().unwrap();
        let filtered = filter(&adaptive, &self.impulse[half][sub]);
        let (rounded, correlation, against) = cross(&target, &filtered);
        let (mut gain, energy) = adaptive_gain(&filtered, rounded, against);
        cap_gain(
            i16::from(periodic(lag.integer, lag.frac, &self.peaks)),
            &mut gain,
        );

        // The two gains the coming subframe's periodic extension may use.  The
        // first is always a half; the second only survives while the pitch is
        // carrying over.
        let reflection = lsp::reflection_coefficients(&synthesis_filter)[0];
        self.extension = interpolation_gains(
            reflection,
            self.carried_gain,
            lag.integer,
            self.previous_lag,
        );

        let coding = self.mode.choose(&target, &filtered, gain);
        let weights = self.quantiser.choose(
            &self.weighting_lpc[half][sub],
            &target,
            &self.subframe_lsf[half][sub],
        );

        let chosen = pulses::search(
            &pulses::Search {
                target: &target,
                contribution: &filtered,
                impulse: &self.impulse[half][sub],
                lag: lag.integer,
                fraction: lag.frac,
                extension: self.extension,
                coding,
                weights: weights as usize,
            },
            &mut self.adaptive_correlations,
        );
        self.field[field + bitstream::CODE] = chosen.index;

        // The innovation itself is built the way the decoder will build it,
        // sharpening included: which of the two gains that uses depends on the
        // part of the codebook the index landed in.
        let decoded = codebook::decode(chosen.index as u16);
        let mut innovation = decoded.code;
        if lag.integer < SUB as i16 {
            let sharpen = if decoded.class == 5 || (1..3).contains(&decoded.class) {
                self.extension.1
            } else {
                self.extension.0
            };
            pitch::sharpen(&mut innovation, &lag, sharpen);
        }
        let shaped = filter(&innovation, &self.impulse[half][sub]);
        let taken = measures(&shaped, &target, &filtered);

        let (scale, gain_exponent) = crate::gain::predict_code_gain(&self.gain_state, &innovation);
        let mut mantissa = [
            energy.mantissa,
            correlation.mantissa,
            taken[0].mantissa,
            taken[1].mantissa,
            taken[2].mantissa,
        ];
        let exponent = [
            energy.exponent,
            correlation.exponent,
            taken[0].exponent,
            taken[1].exponent,
            taken[2].exponent,
        ];
        let shifts = align(&mut mantissa, &exponent, gain_exponent);
        let (entry, pitch_gain, code_gain) = gain_search(scale, &mantissa, &shifts, gain_exponent);
        self.field[field + bitstream::GAIN] = entry as i16;

        // The gain history the next subframe predicts against is the decoder's,
        // so it is kept by running the decoder's own routine.
        crate::gain::decode(&mut self.gain_state, entry as i16, &innovation);
        self.carried_gain = carry_gain(pitch_gain);

        // What the decoder will excite the synthesis filter with.
        excite(
            &mut self.excitation[at..at + SUB],
            &innovation,
            pitch_gain,
            code_gain,
        );
        track_peak(&mut self.peaks, lag.integer, pitch_gain);

        let mut reconstructed = [0i16; SUB];
        synthesis(
            &synthesis_filter,
            &self.excitation[at..at + SUB],
            &mut reconstructed,
            &mut self.speech_memory,
        );

        let tail = SUB - ORDER;
        let (error, residue) = update_memories(
            &self.history[start + tail..start + SUB],
            &reconstructed[tail..],
            &target[tail..],
            &filtered[tail..],
            &shaped[tail..],
            (pitch_gain, code_gain),
        );
        self.error_memory = error;
        self.residual_memory = residue;
        self.half_lag = lag.integer;
    }

    /// Search the excitation history for the lag that best predicts the target.
    ///
    /// The correlations are taken over a short bracket around the open-loop
    /// estimate, with four lags either side so that the fractional
    /// interpolation has something to lean on.  Only the first subframe of a
    /// half-frame codes a fraction, and only while the lag is short enough for
    /// one to be worth sending.
    fn closed_loop_lag(&mut self, half: usize, sub: usize, target: &[i16; SUB]) -> pitch::Lag {
        // The first subframe of a half-frame searches around the open-loop
        // estimate; the second searches the ten lags its own field can code,
        // which is the window the first one's choice opens up.
        let (low, high) = if sub == 0 {
            self.bracket[half]
        } else {
            let window = first_window(self.first_lag, 0);
            (window.low, window.high)
        };
        let lags = (high - low) as usize + LAGS_OVER;
        let at = HISTORY_EXC + sub * SUB;
        let start = at + LAG_MARGIN - low as usize;

        // The convolution windows the impulse response to forty taps and
        // reads the excitation as far as it needs to, which is the other way
        // round from how the fixed codebook's shapes are filtered.
        let mut filtered = filter(
            &self.impulse[half][sub][..REACH_TAPS],
            &self.excitation[start..start + SUB],
        );
        let down = prescale(&mut filtered);

        let past: Vec<i16> = (0..lags).map(|k| self.excitation[start - 1 - k]).collect();
        let correlations = lag_correlations(
            target,
            &mut filtered,
            &past,
            &self.impulse[half][sub],
            down,
            lags,
        );

        let integer = best_closed_lag(&correlations[LAG_MARGIN..], low, high);
        if sub == 0 && integer > COARSE_LAG {
            return pitch::Lag { integer, frac: 0 };
        }

        // The five sub-sample positions, keeping the last of any equals.
        let centre = LAG_MARGIN + (integer - low) as usize;
        let mut best = interpolate(&correlations, centre, -2);
        let mut frac = -2;
        for f in -1..=2 {
            let here = interpolate(&correlations, centre, f);
            if acc(here - best) >= 0 {
                best = here;
                frac = f;
            }
        }

        // The outermost positions belong to the neighbouring integer lag.
        match frac {
            -2 => pitch::Lag {
                integer: integer - 1,
                frac: 1,
            },
            2 => pitch::Lag {
                integer: integer + 1,
                frac: -1,
            },
            _ => pitch::Lag { integer, frac },
        }
    }

    /// Encode one frame of 320 mu-law samples into an 18-byte transport frame.
    ///
    /// This is [`frame`](Self::frame) followed by the bit packing, and is what
    /// [`Decoder::decode`](crate::Decoder::decode) reads back.
    pub fn encode(&mut self, frame: &[u8; FRAME]) -> [u8; bitstream::FRAME_BYTES] {
        bitstream::to_bytes(&bitstream::pack(&self.frame(frame)))
    }

    /// Encode one frame into its fourteen parameter fields: the two line
    /// spectrum fields plus the twelve the two half-frames' subframe searches
    /// contribute.
    ///
    /// [`encode`](Self::encode) is the same thing packed into bytes; this is
    /// the form to use when the fields themselves are of interest.
    pub fn frame(&mut self, frame: &[u8; FRAME]) -> Params {
        let shaped = self.shape(frame);
        let (lsp, a) = self.line_spectrum(&shaped);
        let lsf = lsp::lsp_to_lsf(&lsp);

        let chosen = quantise(&lsf, &weights(&lsf), &self.lsf.history);
        let (first, second) = chosen.fields();
        // Track what the decoder will make of those indices, so that the
        // predictor memory the next frame quantises against is the same one.
        // The interpolation test is on the filter the decoder will actually
        // have, so the fields are taken back through the decoder's own path.
        let rebuilt = lsp::decode(&mut self.lsf, LsfIndices::unpack(first, second));
        let decoded = lsp::lsf_to_lsp(&rebuilt);
        self.decoded_lpc = lsp::lsp_to_lpc(&decoded);
        let reflection = lsp::reflection_coefficients(&self.decoded_lpc);
        let (interpolating, take_current) =
            lsp::interpolation_control(&reflection, self.interpolating);
        self.interpolated = if take_current {
            lsp
        } else {
            midpoint(&self.previous, &lsp)
        };
        self.second_lpc = if take_current {
            a
        } else {
            lsp::lsp_to_lpc(&self.interpolated)
        };
        self.previous = lsp;
        self.interpolating = interpolating;

        // The same ladder on the decoded spectrum, for the synthesis filters.
        let decoded_interpolated = if take_current {
            decoded
        } else {
            midpoint(&self.decoded_previous, &decoded)
        };
        self.decoded_previous = decoded;
        for (half, &current) in [decoded_interpolated, decoded].iter().enumerate() {
            let first = midpoint(&self.decoded_last_half, &current);
            self.synthesis_lpc[half] = [lsp::lsp_to_lpc(&first), lsp::lsp_to_lpc(&current)];
            self.decoded_last_half = current;
        }

        // Each half-frame codes its first subframe against the midpoint of what
        // the last one used and its own line spectrum, and its second against
        // its own.
        for (half, &current) in [self.interpolated, lsp].iter().enumerate() {
            let first = midpoint(&self.last_half, &current);
            self.subframe_lsp[half] = [first, current];
            self.subframe_lsf[half] = [lsp::lsp_to_lsf(&first), lsp::lsp_to_lsf(&current)];
            self.last_half = current;
            self.factors[half] = self.sharpness.choose(
                &self.reflection,
                &[&self.subframe_lsf[half][0], &self.subframe_lsf[half][1]],
            );

            // Each subframe's prediction filter, pulled in twice: once by the
            // numerator factor and once by the denominator's.
            let (numerator, denominator) = self.factors[half];
            for sub in 0..2 {
                // The second subframe is weighted with a filter that does not
                // come from its own line spectrum: the first half-frame uses
                // the interpolated one and the second the analysis one.
                let lpc = match (sub, half) {
                    (0, _) => lsp::lsp_to_lpc(&self.subframe_lsp[half][0]),
                    (_, 0) => self.second_lpc,
                    (_, _) => a,
                };
                self.weighting[half][sub] =
                    (expand(&lpc, numerator[sub]), expand(&lpc, denominator[sub]));
                self.weighting_lpc[half][sub] = lpc;
            }

            // The half-frame's shaping signal, behind the tail of the last one.
            self.history.copy_within(HISTORY - CARRY.., 0);
            self.history[CARRY..].copy_from_slice(&shaped[half]);

            // Three passes over the half-frame: the first subframe's filter
            // over a window straddling the boundary, then the second's over the
            // next eighty samples and the forty after that.
            let mut memory = self.weighting_memory;
            for (pass, &(start, len)) in [
                (CARRY - WINDOW, SUB),
                (CARRY + WINDOW, SUB),
                (CARRY + WINDOW + SUB, SUB / 2),
            ]
            .iter()
            .enumerate()
            {
                let sub = usize::from(pass > 0);
                let (numerator, denominator) = &self.weighting[half][sub];
                let mut through_storage = [0i16; SUB];
                let through = &mut through_storage[..len];
                inverse_filter(
                    numerator,
                    &self.history[start - ORDER..start + len],
                    through,
                );
                // The first two passes carry the denominator's memory forward;
                // the third runs the shorter routine, which leaves it alone.
                let mut out_storage = [0i16; SUB];
                let out = &mut out_storage[..len];
                let mut carried = memory;
                synthesis(denominator, through, out, &mut carried);
                if pass < 2 {
                    memory = carried;
                }
                // The searches run over the weighted signal, so each pass is
                // laid down behind the last.
                let at = SEARCH_SPAN - 2 * SUB - SUB / 2 + pass * SUB;
                let room = (SEARCH_SPAN - at).min(out.len());
                self.search[at..at + room].copy_from_slice(&out[..room]);
            }
            self.weighting_memory = memory;

            self.open_loop[half] = open_loop_lag(&self.search);
            self.bracket[half] = search_range(self.open_loop[half]);
            self.search.copy_within(HOP.., 0);

            for sub in 0..2 {
                self.subframe(half, sub);
            }
            // The lag the coming half-frame tests its own against is the one
            // this half-frame finished on, not the previous subframe's.
            self.previous_lag = self.half_lag;

            // What the next half-frame's closed loop looks back into.
            self.excitation.copy_within(HOP.., 0);
        }

        self.field[0] = first;
        self.field[1] = second;
        Params {
            field: self.field,
            suppress: false,
        }
    }
}
