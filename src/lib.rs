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
        let mut threshold = 0i16;
        let mut counter = 0i16;
        let mut history = self.dec.gain_history;

        for i in 0..4 {
            let block_phase = if i % 2 == 0 { 0i16 } else { 80 };
            let _ = block_phase;
            let lag_int = pitch_lags[i].lag;
            let sub_lag = pitch_lags[i].sub_lag;

            // 6a. Adaptive codebook (pitch_adaptive_codebook)
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

            // 6b. Fixed codebook (dispatch + short or main path)
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

            // 6c. Gain pipeline orchestrator
            let gain_out = gain::gain_orchestrate_codec(&gain::GainOrchestrateCodecInput {
                phase_word: gain_phases[i],
                threshold_in: threshold,
                counter_in: counter,
                candidate: c,
                history,
            });
            let pitch_gain = gain_out.pitch_gain;
            let fcb_gain = gain_out.fcb_gain;
            history = gain_out.history_out;
            threshold = gain_out.threshold_out;
            counter = gain_out.counter_out;
            // commit_pitch_gain — used by the next sub's compute_pitch_enhance_gain.
            self.dec.prev_pitch_gain = pitch_gain;

            // 6d. Excitation mix
            let e = synth::mix_excitation(&v, &c, pitch_gain, fcb_gain);
            // Write back to past_excitation
            for n in 0..80 {
                self.dec.past_excitation[sub_write_offset + n] = e[n];
            }

            // 6e. LPC synthesis filter (1/A(z))
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

            // 6f. Block-end past_excitation shift (after sub 1 and sub 3)
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

        // 8. PCM → μ-law
        let mut out = [0u8; PCM_FRAME_BYTES];
        for i in 0..PCM_FRAME_BYTES {
            out[i] = ulaw::linear_i16_to_ulaw(postfiltered[i]);
        }
        Some(out)
    }
}
