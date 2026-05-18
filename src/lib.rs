//! Mitsubishi CELP speech codec.

pub mod arith;
pub mod bitstream;
pub mod fcb;
pub mod gain;
pub mod lpc_analysis;
pub mod lsf;
pub mod pitch;
pub mod postfilter;
pub mod state;
pub mod synth;
pub mod tables;
pub mod ulaw;

pub use bitstream::{ControlFrame, FIELD_WIDTHS, NUM_FIELDS, canonicalize_payload, is_reset_frame};

pub const MCELP_FRAME_BYTES: usize = 18;
pub const PCM_FRAME_BYTES: usize = 320;
pub const HALF_FRAME_SAMPLES: usize = PCM_FRAME_BYTES / 2;
pub const SUBFRAME_SAMPLES: usize = HALF_FRAME_SAMPLES / 2;
pub const LPC_ORDER: usize = 10;

/// Self-healing watchdog: per-frame count of postfilter output samples
/// saturated at ±i16::MAX above which the decoder assumes it has been
/// driven into a divergent state by corrupt input (e.g. a frame that
/// slipped past upstream integrity checks). 50% of the frame — orders of
/// magnitude above what natural speech can produce.
const WATCHDOG_SATURATION_THRESHOLD: usize = 160;
/// Self-healing watchdog: number of subsequent frames to force the
/// `suppress` concealment path after a trigger, giving the postfilter
/// IIR delay line time to leak out its saturation limit cycle. The
/// hold extends naturally if the symptom keeps recurring.
const WATCHDOG_SUPPRESS_HOLD_FRAMES: u8 = 8;

/// Public entry point that decodes one frame: 18 bytes → 320 bytes of μ-law.
pub struct Codec {
    dec: state::DecoderState,
}

impl Default for Codec {
    fn default() -> Self {
        Self::new()
    }
}

impl Codec {
    pub fn new() -> Self {
        Self {
            dec: state::DecoderState::new(),
        }
    }

    /// Decode an 18-byte MCELP frame into 320 bytes of μ-law PCM.
    ///
    /// On detecting an in-band reset, reset internal state and return `None`.
    pub fn decode_frame(
        &mut self,
        frame: &[u8; MCELP_FRAME_BYTES],
    ) -> Option<[u8; PCM_FRAME_BYTES]> {
        if bitstream::is_reset_frame(frame) {
            self.dec.reset();
            return None;
        }
        let canon = bitstream::canonicalize_payload(frame);
        let ctrl = ControlFrame::unpack(&canon);

        // Self-healing watchdog: if a previous frame tripped the saturation
        // detector, force the suppress concealment path for this frame too.
        let watchdog_engaged = self.dec.suppress_hold_counter > 0;
        if watchdog_engaged {
            self.dec.suppress_hold_counter -= 1;
        }
        let effective_suppress = ctrl.suppress || watchdog_engaged;

        // 1. LSF → LSP
        let lsf_idx = lsf::LsfIndices::from_control(&ctrl);
        let lsf_q15 = lsf::decode_lsf_direct_mode(&lsf_idx, &mut self.dec.lsf_history);
        let lsp_q15 = lsf::lsf_to_lsp(&lsf_q15);

        // 2. Block-level LSP coeffs
        let prev_frame_lsp_at_entry = self.dec.prev_frame_lsp;
        let block_0_coeffs = lsf::block_average_lsp(&prev_frame_lsp_at_entry, &lsp_q15);
        self.dec.prev_frame_lsp = lsp_q15;

        // 3. LPC analysis pipeline:
        //    stability-check pass → block-0 autocorrelation → block-1 autocorrelation.
        //    Each pass carries the response_shaper dword scratch buffer forward.

        // Step 3a: stability-check pass — coeff_seed = lsp + 4 trailing zeros
        let mut stability_seed_window = [0i16; 14];
        stability_seed_window[..10].copy_from_slice(&lsp_q15);
        let (stability_autocorr, post_stability_buf) =
            lpc_analysis::build_autocorrelation_with_state(
                &prev_frame_lsp_at_entry,
                &stability_seed_window,
                &self.dec.response_shaper_buffer,
            );
        self.dec.response_shaper_buffer = post_stability_buf;

        // Step 3b: decide block-0 LPC override (stable && prev_flag == 0)
        let stability_lpc: [i16; 10] = std::array::from_fn(|i| stability_autocorr[12 + i]);
        let stability_reflection = lpc_analysis::levinson_recursion(&stability_lpc);
        let stable_now = lpc_analysis::is_lpc_stable(&stability_reflection);
        let override_block_0 = stable_now && self.dec.prev_lpc_stability_flag == 0;
        self.dec.prev_lpc_stability_flag = if stable_now { 1 } else { 0 };

        let block_0_coeffs_eff: [i16; 10] = if override_block_0 {
            lsp_q15
        } else {
            block_0_coeffs
        };

        // Step 3c: block-0 autocorrelation pass (produces sub 0/1 LPC)
        let mut block0_seed_window = [0i16; 14];
        block0_seed_window[..10].copy_from_slice(&block_0_coeffs_eff);
        block0_seed_window[10..14].copy_from_slice(&lsp_q15[0..4]);
        let (autocorr_block0, post_block0_buf) = lpc_analysis::build_autocorrelation_with_state(
            &prev_frame_lsp_at_entry,
            &block0_seed_window,
            &self.dec.response_shaper_buffer,
        );
        self.dec.response_shaper_buffer = post_block0_buf;
        let sub0_lpc: [i16; 10] = std::array::from_fn(|i| autocorr_block0[1 + i]);
        let sub1_lpc: [i16; 10] = std::array::from_fn(|i| autocorr_block0[12 + i]);
        let clamp_sub0 = lpc_analysis::levinson_recursion(&sub0_lpc)[0];
        let clamp_sub1 = lpc_analysis::levinson_recursion(&sub1_lpc)[0];

        // Step 3d: block-1 autocorrelation pass (produces sub 2/3 LPC)
        let block1_seed_window = stability_seed_window;
        let (autocorr_block1, post_block1_buf) = lpc_analysis::build_autocorrelation_with_state(
            &block_0_coeffs_eff,
            &block1_seed_window,
            &self.dec.response_shaper_buffer,
        );
        self.dec.response_shaper_buffer = post_block1_buf;
        let sub2_lpc: [i16; 10] = std::array::from_fn(|i| autocorr_block1[1 + i]);
        let sub3_lpc: [i16; 10] = std::array::from_fn(|i| autocorr_block1[12 + i]);
        let clamp_sub2 = lpc_analysis::levinson_recursion(&sub2_lpc)[0];
        let clamp_sub3 = lpc_analysis::levinson_recursion(&sub3_lpc)[0];

        let lpc_per_sub = [sub0_lpc, sub1_lpc, sub2_lpc, sub3_lpc];
        let clamp_inputs = [clamp_sub0, clamp_sub1, clamp_sub2, clamp_sub3];

        // 5. Decode pitch lags for 4 subframes (block 0 = sub 0/1, block 1 = sub 2/3)
        // Block 0: F[2] sub 0 abs, F[5] sub 1 diff
        // Block 1: F[8] sub 2 abs, F[11] sub 3 diff
        let pl0 =
            pitch::decode_subframe_lag(&mut self.dec.pitch_state_block0, 0, ctrl.fields[2] as i16);
        let pl1 =
            pitch::decode_subframe_lag(&mut self.dec.pitch_state_block0, 80, ctrl.fields[5] as i16);
        let pl2 =
            pitch::decode_subframe_lag(&mut self.dec.pitch_state_block1, 0, ctrl.fields[8] as i16);
        let pl3 = pitch::decode_subframe_lag(
            &mut self.dec.pitch_state_block1,
            80,
            ctrl.fields[11] as i16,
        );
        let pitch_lags = [pl0, pl1, pl2, pl3];

        // FCB index + gain phase fields
        let fcb_indices: [i16; 4] = [
            ctrl.fields[3] as i16,  // sub 0
            ctrl.fields[6] as i16,  // sub 1
            ctrl.fields[9] as i16,  // sub 2
            ctrl.fields[12] as i16, // sub 3
        ];
        let gain_phases: [i16; 4] = [
            ctrl.fields[4] as i16,  // sub 0
            ctrl.fields[7] as i16,  // sub 1
            ctrl.fields[10] as i16, // sub 2
            ctrl.fields[13] as i16, // sub 3
        ];

        // 6. Per-subframe synthesis
        let mut synth_pcm = [0i16; 320]; // Q15 PCM
        // gain pipeline persistent state (frame-local)
        let mut threshold = self.dec.gain_threshold;
        let mut counter = self.dec.gain_counter;
        let mut history = self.dec.gain_history;
        let mut lag_cursor = self.dec.lag_cursor;
        let mut fcb_gain_state = self.dec.prev_fcb_gain;

        for i in 0..4 {
            // 6a. Lag selection
            // - suppress=0: Use pitch_lags[i] already decoded by decode_subframe_lag,
            //   and update lag_cursor with selected_lag.
            // - suppress=1: Do not call ce3e; use lag_cursor as the selected lag.
            //   sub_lag is 0, cursor is +1 for the next sub (cap 143).
            let (lag_int, sub_lag) = if effective_suppress {
                let lag = lag_cursor;
                lag_cursor = ((lag_cursor as i32) + 1).min(143) as i16;
                (lag, 0i16)
            } else {
                let pl = pitch_lags[i];
                lag_cursor = pl.lag;
                (pl.lag, pl.sub_lag)
            };

            // 6b. Adaptive codebook (pitch_adaptive_codebook)
            let (fract, lag_adjust) = pitch::decode_lag_fract(sub_lag);
            let effective_lag = (lag_int as i32).wrapping_add(lag_adjust as i32);
            let write_offset = state::PAST_EXCITATION_OUTPUT_BASE; // = 160
            let sub_in_block = if i % 2 == 0 { 0 } else { 80 };
            let sub_write_offset = write_offset + sub_in_block;
            let base_offset_signed = (sub_write_offset as i32) - effective_lag;

            let mut v = [0i16; 80];
            // Same bounds check as rust_float: out-of-range returns zeros.
            if base_offset_signed >= 9
                && (base_offset_signed + 80 + 10) as usize <= state::PAST_EXCITATION_LEN
                && (sub_write_offset + 80) <= state::PAST_EXCITATION_LEN
            {
                pitch::pitch_adaptive_codebook(
                    &mut self.dec.past_excitation,
                    base_offset_signed as usize,
                    sub_write_offset,
                    fract,
                    &mut v,
                );
            }

            // 6c. Fixed codebook (dispatch + short or main path)
            let dispatch = fcb::fcb_dispatch_lag_class(fcb_indices[i]);
            // pitch_enhance gain — via update_synth_control. `clamp_input` is the
            // output of levinson_recursion (= already computed from each sub's autocorrelation window).
            let clamp_input = clamp_inputs[i];
            let pitch_enhance_gain = synth::compute_pitch_enhance_gain(
                fcb_indices[i],
                lag_int,
                clamp_input,
                self.dec.prev_lag,
                self.dec.prev_pitch_gain,
            );
            let c = if dispatch.lag_class == 5 {
                fcb::fcb_short_path(
                    dispatch.clamped_fcb_index,
                    lag_int,
                    sub_lag,
                    pitch_enhance_gain,
                )
            } else {
                fcb::fcb_main_path(
                    dispatch.lag_class,
                    dispatch.clamped_fcb_index,
                    lag_int,
                    sub_lag,
                    pitch_enhance_gain,
                )
            };

            // 6d. Gain — normal orchestrate or suppress decay
            // Branches based on suppress flag:
            // - Non-zero: gain_suppress_decay (threshold +1, multiply gains by decay factor)
            // - Zero:   Normal orchestrate (= gain_orchestrate_codec)
            let (pitch_gain, fcb_gain) = if effective_suppress {
                let decay =
                    gain::gain_suppress_decay(threshold, self.dec.prev_pitch_gain, fcb_gain_state);
                threshold = decay.threshold_out;
                // For the suppress=1 path, the update logic for counter/history
                // is separate. In the minimal implementation, the current state is maintained.
                (decay.pitch_gain_out, decay.fcb_gain_out)
            } else {
                let gain_out = gain::gain_orchestrate_codec(&gain::GainOrchestrateCodecInput {
                    phase_word: gain_phases[i],
                    threshold_in: threshold,
                    counter_in: counter,
                    candidate: c,
                    history,
                });
                history = gain_out.history_out;
                threshold = gain_out.threshold_out;
                counter = gain_out.counter_out;
                (gain_out.pitch_gain, gain_out.fcb_gain)
            };
            // commit_pitch_gain — used by the next sub's compute_pitch_enhance_gain.
            self.dec.prev_pitch_gain = pitch_gain;
            fcb_gain_state = fcb_gain;

            // 6e. Excitation mix
            // - watchdog engaged: hard-mute (e = 0). The bitstream suppress
            //   path's decayed fcb_gain is still close to i16::MAX, so under
            //   continuous corruption it can keep saturating mix_excitation
            //   and re-feeding the postfilter limit-cycle. Watchdog mode is
            //   panic mode: drain every IIR.
            // - bitstream suppress=1: pitch component zeroed, fcb scaled.
            // - normal: pitch_gain × v + fcb_gain × c.
            let e = if watchdog_engaged {
                [0i16; 80]
            } else if ctrl.suppress {
                synth::mix_excitation(&v, &c, 0, fcb_gain)
            } else {
                synth::mix_excitation(&v, &c, pitch_gain, fcb_gain)
            };
            // Write back to past_excitation
            for n in 0..80 {
                self.dec.past_excitation[sub_write_offset + n] = e[n];
            }

            // 6f. LPC synthesis filter (1/A(z))
            // Pass cells [1..11] / [12..22] of the autocorrelation window from
            // build_autocorrelation_with_state directly as the 10-tap LPC
            // (= Levinson is already done inside response_shaper).
            let lpc10 = &lpc_per_sub[i];
            let (s_hat, new_history) =
                synth::lpc_synthesis_filter(&e, lpc10, &self.dec.lpc_synth_history);
            self.dec.lpc_synth_history = new_history;
            // Copy synthesis output to PCM buffer
            let pcm_offset = i * 80;
            for n in 0..80 {
                synth_pcm[pcm_offset + n] = s_hat[n];
            }

            // 6g. Block-end past_excitation shift (after sub 1 and sub 3)
            if i % 2 == 1 {
                // Shift [160..320] → [0..160], clear [160..320]
                for k in 0..160 {
                    self.dec.past_excitation[k] = self.dec.past_excitation[160 + k];
                }
                for k in 160..320 {
                    self.dec.past_excitation[k] = 0;
                }
                // commit_half_frame_lag — prev_lag used by the next block's pitch_enhance_gain calculation.
                self.dec.prev_lag = lag_int;
            }
        }
        self.dec.gain_history = history;
        self.dec.gain_threshold = threshold;
        self.dec.gain_counter = counter;
        self.dec.lag_cursor = lag_cursor;
        self.dec.prev_fcb_gain = fcb_gain_state;

        // 7. Postfilter (per half-frame × 2)
        let mut postfiltered = [0i16; 320];
        for half in 0..2 {
            let pre: [i16; 160] = std::array::from_fn(|k| synth_pcm[half * 160 + k]);
            let pf = postfilter::postfilter_apply(
                &pre,
                &mut self.dec.postfilter_history,
                &mut self.dec.postfilter_delay,
            );
            for k in 0..160 {
                postfiltered[half * 160 + k] = pf[k];
            }
        }

        // Watchdog: if the postfilter output is heavily saturated, the
        // decoder is in a divergent state (corrupt frame not caught
        // upstream → postfilter IIR saturation limit-cycle). The IIR is
        // self-sustaining once its delay line holds ±i32::MAX dword pairs
        // — even zero input from the synthesis filter cannot break it via
        // the linear dynamics alone — so we must do two things on
        // trigger: (1) zero the postfilter delay/history to break the
        // cycle, and (2) hold the watchdog suppress for the next several
        // frames so the LPC synth filter, past-excitation, and any
        // residual postfilter state all drain cleanly.
        let saturated = postfiltered
            .iter()
            .filter(|&&v| v == i16::MAX || v == i16::MIN)
            .count();
        if saturated >= WATCHDOG_SATURATION_THRESHOLD {
            self.dec.suppress_hold_counter = WATCHDOG_SUPPRESS_HOLD_FRAMES;
            self.dec.postfilter_history = [0; 6];
            self.dec.postfilter_delay = [0; 12];
        }

        // 8. PCM → μ-law
        let mut out = [0u8; PCM_FRAME_BYTES];
        for i in 0..PCM_FRAME_BYTES {
            out[i] = ulaw::linear_i16_to_ulaw(postfiltered[i]);
        }
        Some(out)
    }
}
