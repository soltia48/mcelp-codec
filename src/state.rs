//! Codec-wide persistent state.

use crate::lsf::LsfHistory;
use crate::pitch::PitchLagState;

/// Length of the past-excitation sliding buffer (= previous 160 + current 160 samples).
/// Layout:
///   `[0..160]`   = previous block's output
///   `[160..240]` = current block's sub a (phase=0) write region
///   `[240..320]` = current block's sub b (phase=80) write region
/// At the end of each block, shift `[160..320]` down to `[0..160]`.
pub const PAST_EXCITATION_LEN: usize = 320;
/// Offset on the past-excitation buffer where sub writes start (= start of the current block's region).
pub const PAST_EXCITATION_OUTPUT_BASE: usize = 160;

#[derive(Clone, Debug)]
pub struct DecoderState {
    /// 3-slot rolling history for LSF decoding. Shifted+inserted each frame.
    pub lsf_history: LsfHistory,
    /// Previous-frame LSP (post-`lsf_to_lsp` output). Used as the `prev` input
    /// of `block_average_lsp` for block 0. Initial value: `tables::lsf::INIT_LSP_TEMPLATE`.
    pub prev_frame_lsp: [i16; 10],
    /// block_coeffs of the previous half-frame (= `prev` argument of `subframe_lsp_interp`).
    pub prev_block_coeffs: [i16; 10],
    /// Sliding past-excitation buffer for the adaptive codebook (`pitch_adaptive_codebook`).
    pub past_excitation: [i16; PAST_EXCITATION_LEN],
    /// Inter-subframe history for `lpc_synthesis_filter` (last 10 samples).
    pub lpc_synth_history: [i16; 10],
    /// Pitch-lag state for block 0 (sub 0/1).
    pub pitch_state_block0: PitchLagState,
    /// Pitch-lag state for block 1 (sub 2/3).
    pub pitch_state_block1: PitchLagState,
    /// 5-cell past-gain history used by the gain orchestrator for prediction.
    /// codec init: `[-17254, -17254, -17254, -17254, 0]`.
    pub gain_history: [i16; 5],
    /// Suppress-path threshold of the gain orchestrator. codec init: 0.
    pub gain_threshold: i16,
    /// Tail decay counter of the gain orchestrator. codec init: 0.
    pub gain_counter: i16,
    /// Postfilter history shift-register state (6 cells).
    /// `[0..3]` = forward history, `[3..6]` = reverse history. codec init: all zero.
    pub postfilter_history: [i16; 6],
    /// Postfilter delay-line window (12 cells).
    /// `[0..6]` = forward delay, `[6..12]` = reverse delay. codec init: all zero.
    pub postfilter_delay: [i16; 12],
    /// Last selected_lag of the previous half-frame. Used by `compute_pitch_enhance_gain`.
    /// codec init: 60.
    pub prev_lag: i16,
    /// Post-orchestrate pitch_gain of the previous subframe. Used by
    /// `compute_pitch_enhance_gain`. codec init: 3277.
    pub prev_pitch_gain: i16,
    /// `response_shaper` dword scratch buffer. Called 3 times per frame
    /// (stability-check pass → block 0 → block 1); state is carried across passes.
    /// codec init: all zero.
    pub response_shaper_buffer: [i16; 35],
    /// LPC stability flag of the previous frame. Used as the condition
    /// `stable && prev_flag == 0` to overwrite block-0 coeffs with the LSP.
    pub prev_lpc_stability_flag: i16,
}

impl Default for DecoderState {
    fn default() -> Self {
        Self::new()
    }
}

impl DecoderState {
    pub fn new() -> Self {
        Self {
            lsf_history: LsfHistory::new(),
            prev_frame_lsp: crate::tables::lsf::INIT_LSP_TEMPLATE,
            prev_block_coeffs: [0; 10],
            past_excitation: [0; PAST_EXCITATION_LEN],
            lpc_synth_history: [0; 10],
            pitch_state_block0: PitchLagState::default(),
            pitch_state_block1: PitchLagState::default(),
            gain_history: [-17254, -17254, -17254, -17254, 0],
            gain_threshold: 0,
            gain_counter: 0,
            postfilter_history: [0; 6],
            postfilter_delay: [0; 12],
            prev_lag: 60,
            prev_pitch_gain: 3277,
            response_shaper_buffer: [0; 35],
            prev_lpc_stability_flag: 0,
        }
    }

    pub fn reset(&mut self) {
        *self = Self::new();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_state_has_template_lsp_and_zero_excitation() {
        let s = DecoderState::new();
        assert_eq!(s.prev_frame_lsp, crate::tables::lsf::INIT_LSP_TEMPLATE);
        assert!(s.past_excitation.iter().all(|&v| v == 0));
        assert!(s.lpc_synth_history.iter().all(|&v| v == 0));
    }

    #[test]
    fn reset_returns_to_initial_state() {
        let mut s = DecoderState::new();
        s.past_excitation[0] = 1234;
        s.lpc_synth_history[5] = -5678;
        s.reset();
        assert!(s.past_excitation.iter().all(|&v| v == 0));
        assert!(s.lpc_synth_history.iter().all(|&v| v == 0));
    }
}
