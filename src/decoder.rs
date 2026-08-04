//! The decoder: everything it carries between frames, and the frame loop
//! that ties the stages together.
//!
//! A frame is decoded half at a time.  Each half re-derives its two subframe
//! LPC sets by interpolating LSPs, builds and synthesises the two subframes,
//! then postfilters them and runs the result through the output filter.  The
//! synthesised-speech and excitation buffers are only one half-frame long and
//! are slid along as each half completes.

use crate::bitstream::{self, Params};
use crate::codebook;
use crate::fixed::{acc, hi, sat, shift};
use crate::gain::{self, GainState};
use crate::lsp::{self, LsfState, ORDER};
use crate::pitch::{self, Lag, SUBFRAME};
use crate::postfilter;
use crate::synth;
use crate::ulaw;

/// Samples per frame.
pub const FRAME: usize = 4 * SUBFRAME;
/// Excitation history kept before the current half-frame.
pub const EXC_HISTORY: usize = 154;
/// Half a frame.
pub const HALF: usize = 2 * SUBFRAME;
/// Residual history the pitch postfilter searches over.
pub const RES_HISTORY: usize = 152;

/// Everything that has to survive from one frame to the next.
#[derive(Clone)]
pub struct Decoder {
    /// LSF quantiser memory.
    pub(crate) lsf: LsfState,
    /// Gain quantiser memory.
    pub(crate) gain: GainState,
    /// LSP set of the previous half-frame, one end of the interpolation.
    pub(crate) prev_lsp: [i16; ORDER],
    /// `1/A(z)` filter memory.
    pub(crate) synth_mem: [i16; ORDER],
    /// Excitation: `EXC_HISTORY` past samples followed by the current half.
    pub(crate) exc: [i16; EXC_HISTORY + HALF],
    /// Synthesised speech: ten past samples followed by the current *half*
    /// frame, which is all the buffer the reference keeps.
    pub(crate) speech: [i16; ORDER + HALF],
    /// `1/A(z/g1)` memory of the postfilter's re-synthesis.
    pub(crate) pf_synth_mem: [i16; ORDER],
    /// Sample immediately before the postfiltered block, for tilt compensation.
    pub(crate) pf_last: i16,
    /// Smoother state of the postfilter's automatic gain control.
    pub(crate) agc_gain: i16,
    /// Output reconstruction filter.
    pub(crate) output_filter: crate::postfilter::OutputFilter,
    /// Pitch lag of the half-frame's first subframe, which the long-term
    /// postfilter searches around.
    pub(crate) reference_lag: i16,
    /// Last subframe lag of the previous half-frame.
    pub(crate) last_half_lag: i16,
    /// LPC residual the postfilter works on: history plus the current subframe.
    pub(crate) residual: [i16; RES_HISTORY + SUBFRAME],
    /// Sharpening gain used by fixed-codebook classes 0, 3 and 4.
    pub(crate) sharpen_fixed: i16,
    /// Sharpening gain used by classes 1, 2 and 5; doubles as a voicing flag.
    pub(crate) sharpen_voiced: i16,
    /// Pitch lag carried across frames for erasure concealment.
    pub(crate) conceal_lag: i16,
    /// Pitch lag of the previous subframe, for the lag-continuity test.
    pub(crate) prev_lag: i16,
    /// Latest non-zero long-term postfilter decision, read by the next
    /// half-frame's excitation build.
    pub(crate) voiced: i16,
    /// Whether LSP interpolation is currently disabled.
    pub(crate) no_interp: i16,
    /// State of the 32-bit LCG used to fill erased frames.
    pub(crate) rng: u32,
}

impl Default for Decoder {
    /// The reset state.
    fn default() -> Self {
        let mut prev_lsp = [0i16; ORDER];
        prev_lsp.copy_from_slice(&crate::tables::LSP_MEAN[..ORDER]);
        Decoder {
            lsf: LsfState::default(),
            gain: GainState::default(),
            prev_lsp,
            synth_mem: [0; ORDER],
            exc: [0; EXC_HISTORY + HALF],
            speech: [0; ORDER + HALF],
            pf_synth_mem: [0; ORDER],
            pf_last: 0,
            agc_gain: 16384,
            output_filter: crate::postfilter::OutputFilter::default(),
            reference_lag: 60,
            last_half_lag: 60,
            residual: [0; RES_HISTORY + SUBFRAME],
            sharpen_fixed: 3277,
            sharpen_voiced: 3277,
            conceal_lag: 60,
            prev_lag: 60,
            voiced: 60,
            no_interp: 0,
            rng: 21845,
        }
    }
}

/// Lower clamp on the pitch gain when it is reused for pitch sharpening, Q14.
const SHARPEN_MIN: i16 = 3277;
/// Upper clamp on the same, Q14.
const SHARPEN_MAX: i16 = 13017;
/// The "definitely voiced" sharpening gain, Q15.
const SHARPEN_VOICED: i16 = 16384;
/// Multiplier of the pseudo-random generator used to fill erased frames.
const RNG_MULTIPLIER: u16 = 31821;
/// Increment of the same generator.
const RNG_INCREMENT: i16 = 13849;

impl Decoder {
    /// A decoder in its reset state.
    pub fn new() -> Self {
        Self::default()
    }

    /// Decode one transport frame into 320 mu-law samples.
    ///
    /// Returns `None` for an in-band reset frame, which resets the decoder
    /// instead of producing audio; every other frame yields a full frame of
    /// samples, including one flagged as erased, which is concealed.
    pub fn decode(&mut self, frame: &[u8; bitstream::FRAME_BYTES]) -> Option<[u8; FRAME]> {
        let linear = self.decode_linear(frame)?;
        let mut out = [0u8; FRAME];
        ulaw::frame_from_linear(&linear, &mut out);
        Some(out)
    }

    /// The same, as linear 16-bit samples rather than mu-law.
    ///
    /// The codec is mu-law end to end, so [`decode`](Self::decode) is what the
    /// reference produces; this returns the samples one step earlier, before
    /// the output companding, for callers that want linear PCM.
    pub fn decode_linear(&mut self, frame: &[u8; bitstream::FRAME_BYTES]) -> Option<[i16; FRAME]> {
        let words = bitstream::canonicalize(frame);
        if bitstream::is_reset(&words) {
            *self = Decoder::default();
            return None;
        }
        let params = bitstream::unpack(&words);
        Some(self.decode_params(&params))
    }

    /// Decode from already unpacked parameters into linear samples.
    pub fn decode_params(&mut self, params: &Params) -> [i16; FRAME] {
        let (first_half_lsp, second_half_lsp) = self.decode_spectrum(params);

        let mut linear = [0i16; FRAME];
        for half in 0..2 {
            let half_lsp = if half == 0 {
                first_half_lsp
            } else {
                second_half_lsp
            };
            let subframe_lsp = self.decode_half(params, half, &half_lsp);

            // Postfilter each subframe, then the half as a whole.
            self.voiced = 0;
            for (sub, lsp) in subframe_lsp.iter().enumerate() {
                let at = sub * SUBFRAME;
                let lag = self.postfilter_subframe(
                    lsp,
                    at,
                    &mut linear[half * HALF + at..half * HALF + at + SUBFRAME],
                );
                if lag != 0 {
                    self.voiced = lag;
                }
            }

            // The last ten samples of the half become the next half's history.
            let tail = ORDER + HALF - ORDER;
            self.speech.copy_within(tail..tail + ORDER, 0);

            self.output_filter
                .run(&mut linear[half * HALF..(half + 1) * HALF]);
        }

        linear
    }

    /// Test hook: run the spectral stage in isolation.
    #[doc(hidden)]
    pub fn dump_spectrum(&mut self, params: &Params) -> ([i16; ORDER], [i16; ORDER]) {
        self.decode_spectrum(params)
    }

    /// Test hook: run one half-frame in isolation.
    /// Test hook: the persistent state, formatted for diffing against the
    /// instrumented reference.
    #[doc(hidden)]
    pub fn dump_state(&self) -> String {
        format!(
            "voiced={} gains={:?} erasures={} lag={} rng={:08x}\nEXC {:?}\nSPEECH {:?}",
            self.voiced,
            self.gain.gains,
            self.gain.erasures,
            self.conceal_lag,
            self.rng,
            &self.exc[..],
            &self.speech[..],
        )
    }

    #[doc(hidden)]
    pub fn dump_half(&mut self, params: &Params, half: usize, lsp: &[i16; ORDER]) {
        self.decode_half(params, half, lsp);
    }

    /// Decode the frame's spectral parameters.
    ///
    /// Returns the LSP sets the two half-frames should aim at.
    fn decode_spectrum(&mut self, params: &Params) -> ([i16; ORDER], [i16; ORDER]) {
        let lsf = if params.suppress {
            lsp::decode_suppressed(&mut self.lsf)
        } else {
            let idx = lsp::LsfIndices::unpack(params.field[0], params.field[1]);
            lsp::decode(&mut self.lsf, idx)
        };
        self.lsf.prev_lsf = lsf;
        let current = lsp::lsf_to_lsp(&lsf);

        let mut mid = lsp::midpoint(&self.prev_lsp, &current);
        let reflection = lsp::reflection_coefficients(&lsp::lsp_to_lpc(&current));
        let (flag, force_current) = lsp::interpolation_control(&reflection, self.no_interp != 0);
        if force_current {
            mid = current;
        }
        self.no_interp = flag as i16;
        (mid, current)
    }

    /// Build and synthesise one half-frame.
    ///
    /// Returns the per-subframe LSP and LPC sets, which the postfilter needs.
    fn decode_half(
        &mut self,
        params: &Params,
        half: usize,
        half_lsp: &[i16; ORDER],
    ) -> [[i16; ORDER]; 2] {
        let (lsp_a, lsp_b, lpc_a, lpc_b) = lsp::interpolate_pair(&self.prev_lsp, half_lsp);
        self.prev_lsp = *half_lsp;
        let subframe_lsp = [lsp_a, lsp_b];
        let subframe_lpc = [lpc_a, lpc_b];

        for (sub, subframe_lpc) in subframe_lpc.iter().enumerate() {
            let at = EXC_HISTORY + sub * SUBFRAME;
            let field = params.subframe(half, sub);

            // ── pitch lag ────────────────────────────────────────────────
            let lag = if params.suppress {
                let l = Lag {
                    integer: self.conceal_lag,
                    frac: 0,
                };
                self.conceal_lag = (self.conceal_lag + 1).min(143);
                l
            } else if sub == 0 {
                let l = pitch::decode_absolute(field.lag);
                self.conceal_lag = l.integer;
                l
            } else {
                let l = pitch::decode_relative(field.lag, self.prev_lag);
                self.conceal_lag = l.integer;
                l
            };
            if sub == 0 {
                self.reference_lag = lag.integer;
            }
            self.prev_lag = lag.integer;

            // ── adaptive codebook ────────────────────────────────────────
            pitch::predict(&mut self.exc, at, &lag);

            // ── voicing test, which selects the sharpening gain ──────────
            self.sharpen_voiced = self.sharpen_fixed;
            let reflection = lsp::reflection_coefficients(subframe_lpc);
            self.sharpen_fixed = SHARPEN_VOICED;
            if acc(shift((reflection[0] as i64) << 16, 3) - (16384i64 << 16)) > 0 {
                self.sharpen_voiced = 0;
            } else if acc(shift((self.sharpen_voiced as i64) << 16, 1) - (13107i64 << 16)) > 0 {
                let lag_acc = (lag.integer as i64) << 16;
                let floor = acc(lag_acc - (self.last_half_lag as i64) * 24576 * 2);
                let ceiling = acc(lag_acc - shift(acc((self.last_half_lag as i64) * 24576 * 2), 1));
                if floor > 0 && ceiling < 0 {
                    self.sharpen_voiced = SHARPEN_VOICED;
                }
            }

            // ── fixed codebook ───────────────────────────────────────────
            let code_index = if params.suppress {
                self.next_random()
            } else {
                field.code as u16
            };
            let innovation = codebook::decode(code_index);
            let mut code = innovation.code;
            if lag.integer < SUBFRAME as i16 {
                let sharpen = if innovation.class == 5 || (1..3).contains(&innovation.class) {
                    self.sharpen_voiced
                } else {
                    self.sharpen_fixed
                };
                pitch::sharpen(&mut code, &lag, sharpen);
            }

            // ── gains ────────────────────────────────────────────────────
            let gains = if params.suppress {
                gain::decode_suppressed(&mut self.gain)
            } else {
                gain::decode(&mut self.gain, field.gain, &code)
            };
            self.sharpen_fixed = gains.pitch.clamp(SHARPEN_MIN, SHARPEN_MAX);

            // ── excitation ───────────────────────────────────────────────
            self.combine_excitation(params.suppress, at, &code, gains);

            // ── short-term synthesis ─────────────────────────────────────
            let mut out = [0i16; SUBFRAME];
            let Decoder { exc, synth_mem, .. } = self;
            synth::synthesis(subframe_lpc, &exc[at..at + SUBFRAME], &mut out, synth_mem);
            let dst = ORDER + sub * SUBFRAME;
            self.speech[dst..dst + SUBFRAME].copy_from_slice(&out);
        }

        // Slide the excitation history along by one half-frame.
        self.exc.copy_within(HALF.., 0);
        self.last_half_lag = self.prev_lag;
        subframe_lsp
    }

    /// Mix the adaptive and fixed contributions into the excitation.
    ///
    /// While frames are being received normally both contributions are summed;
    /// during an erasure only one of them survives, chosen by whether the last
    /// long-term postfilter found a pitch.
    fn combine_excitation(
        &mut self,
        suppressed: bool,
        at: usize,
        code: &[i16; SUBFRAME],
        gains: gain::Gains,
    ) {
        let voiced = self.voiced != 0;
        for (&innovation, excitation) in code.iter().zip(self.exc[at..at + SUBFRAME].iter_mut()) {
            let value = if !suppressed {
                let fixed = shift(acc((gains.code as i64) * (innovation as i64) * 2), 1);
                let adaptive = acc((gains.pitch as i64) * (*excitation as i64) * 2);
                sat(acc(shift(acc(fixed + adaptive), 1) + (1 << 15)))
            } else if voiced {
                sat(acc(shift(
                    acc((gains.pitch as i64) * (*excitation as i64) * 2),
                    1,
                ) + (1 << 15)))
            } else {
                sat(acc(shift(
                    acc((gains.code as i64) * (innovation as i64) * 2),
                    2,
                ) + (1 << 15)))
            };
            *excitation = hi(value);
        }
    }

    /// Postfilter one subframe and emit it.
    ///
    /// Returns the pitch lag the long-term postfilter settled on, or zero when
    /// it decided to stay out of the way.
    fn postfilter_subframe(
        &mut self,
        subframe_lsp: &[i16; ORDER],
        at: usize,
        out: &mut [i16],
    ) -> i16 {
        let denominator = postfilter::weighted_lpc(subframe_lsp, &crate::tables::PF_LSF_WEIGHT1);
        let numerator = postfilter::weighted_lpc(subframe_lsp, &crate::tables::PF_LSF_WEIGHT2);

        // Residual of the synthesised speech through the numerator.
        let speech = &self.speech[at..at + ORDER + SUBFRAME];
        let mut residual = [0i16; SUBFRAME];
        postfilter::inverse_filter(&numerator, speech, &mut residual);
        self.residual[RES_HISTORY..].copy_from_slice(&residual);

        // Long-term postfilter.
        let mut filtered = [0i16; SUBFRAME];
        let lag = crate::ltp::filter(&self.residual, self.reference_lag, &mut filtered);

        // Short-term postfilter: gain normalisation, then re-synthesis.
        let probe = postfilter::probe_response(&denominator, &numerator);
        let tilt = postfilter::spectral_tilt(&probe);
        let mut block = filtered;
        postfilter::normalise_gain(&probe, &mut block);

        let previous = self.pf_last;
        let mut resynth = [0i16; SUBFRAME];
        synth::synthesis(&denominator, &block, &mut resynth, &mut self.pf_synth_mem);
        self.pf_last = self.pf_synth_mem[ORDER - 1];

        let mut compensated = [0i16; SUBFRAME];
        postfilter::compensate_tilt(tilt, previous, &resynth, &mut compensated);

        let Decoder {
            speech, agc_gain, ..
        } = self;
        let speech_block = &speech[ORDER + at..ORDER + at + SUBFRAME];
        postfilter::agc(speech_block.try_into().unwrap(), &mut compensated, agc_gain);
        out.copy_from_slice(&compensated);

        // Slide the residual history along.
        self.residual.copy_within(SUBFRAME.., 0);
        lag
    }

    /// One step of the linear congruential generator that stands in for a lost
    /// codebook index.
    fn next_random(&mut self) -> u16 {
        // The two halves are multiplied separately; both products pick up the
        // fractional-mode doubling, which the right shift on the low half then
        // takes back out.
        let upper = ((self.rng >> 16) as u16 as i64) * (RNG_MULTIPLIER as i64) * 2;
        let lower = ((self.rng as u16) as i64) * (RNG_MULTIPLIER as i64) * 2;
        let mixed = acc(shift(acc(upper), 15) + shift(lower, -1) + (RNG_INCREMENT as i64));
        self.rng = mixed as u32;
        self.rng as u16
    }
}
