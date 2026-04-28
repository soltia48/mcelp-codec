# Decoder State

The decoder is a streaming, frame-stateful machine: each 18-byte frame
depends on roughly a frame's worth of state inherited from the previous
frame. This document enumerates that state and explains how each piece
flows between frames, blocks, and subframes. Source:
[src/state.rs](../src/state.rs), [src/lib.rs](../src/lib.rs).

## 1. The `DecoderState` struct

```rust
pub struct DecoderState {
    pub lsf_history:               LsfHistory,           // 3-slot LSF rolling history
    pub prev_frame_lsp:            [i16; 10],            // previous frame's LSP (Q15)
    pub prev_block_coeffs:         [i16; 10],            // previous block's LPC work coeffs
    pub past_excitation:           [i16; 320],           // sliding excitation buffer
    pub lpc_synth_history:         [i16; 10],            // last 10 samples of synth output
    pub pitch_state_block0:        PitchLagState,        // prev_lag for sub 0/1
    pub pitch_state_block1:        PitchLagState,        // prev_lag for sub 2/3
    pub gain_history:              [i16; 5],             // 5-cell past-gain ring
    pub gain_threshold:            i16,                  // tail-decay threshold
    pub gain_counter:              i16,                  // tail-decay countdown
    pub postfilter_history:        [i16; 6],             // forward+reverse 4-tap shift regs
    pub postfilter_delay:          [i16; 12],            // forward+reverse 3-section delays
    pub prev_lag:                  i16,                  // last-block selected lag (synth control)
    pub prev_pitch_gain:           i16,                  // last-subframe g_p (synth control)
    pub response_shaper_buffer:    [i16; 35],            // dword scratch carried across passes
    pub prev_lpc_stability_flag:   i16,                  // 0/1 — for block-0 override edge
}
```

The total state is around 800 bytes — small enough to copy cheaply, big
enough to make the codec strictly *stateful*.

## 2. Initial values

```rust
DecoderState::new() = {
    lsf_history          = LsfHistory::new()       // all 3 slots = INIT_LSF_TEMPLATE
    prev_frame_lsp       = INIT_LSP_TEMPLATE
    prev_block_coeffs    = [0; 10]
    past_excitation      = [0; 320]
    lpc_synth_history    = [0; 10]
    pitch_state_block0   = PitchLagState::default()  // prev_lag = 0
    pitch_state_block1   = PitchLagState::default()
    gain_history         = [-17254, -17254, -17254, -17254, 0]
    gain_threshold       = 0
    gain_counter         = 0
    postfilter_history   = [0; 6]
    postfilter_delay     = [0; 12]
    prev_lag             = 60
    prev_pitch_gain      = 3277
    response_shaper_buffer = [0; 35]
    prev_lpc_stability_flag = 0
}
```

A few non-zero seeds are worth calling out:

- `INIT_LSF_TEMPLATE` is a uniformly-spaced 10-LSF vector
  `[2340, 4679, ..., 23396]`, i.e. an even split of `[0, 0.7]` in Q15.
  When the decoder boots, the LSF history holds three copies of this
  template so that `predictive_combine` does not multiply garbage on
  the first frame.
- `INIT_LSP_TEMPLATE` is the cosine of that even split (computed
  offline and tabulated). The first frame's `prev_frame_lsp` is this
  template.
- `gain_history` is seeded with `-17254` repeated 4 times and a
  zero-padded fifth slot. This is a small negative bias that
  approximates the mean past-gain prediction state — without it, the
  first few frames mispredict drastically.
- `prev_lag = 60` and `prev_pitch_gain = 3277` are quiescent values
  that put the synth-control hysteresis in a stable mid-state on
  startup.

## 3. State flow per frame

Within a frame, state is updated in this order:

1. **LSF history** (`lsf_history.splice_new_vector`) — committed
   *during* `decode_lsf_direct_mode`, before stabilization, with the
   raw scratch vector (not the predicted output).
2. **Previous-frame LSP** (`prev_frame_lsp`) — committed at the end of
   stage 1 (after `lsf_to_lsp`), so that the *current* frame's LSP
   becomes the *next* frame's `prev_frame_lsp`.
3. **`response_shaper_buffer`** — written 3 times (once per pass) and
   each time the new value is committed back. The order of the three
   passes within a frame is significant.
4. **`prev_lpc_stability_flag`** — committed after the stability check.
5. **Per-subframe** (i = 0..3):
   - `past_excitation[write_offset + n]` ← `e[n]`
   - `lpc_synth_history` ← last 10 samples of `s_hat`
   - `prev_pitch_gain` ← `gain_out.pitch_gain`
6. **Per-block** (after sub 1 and sub 3):
   - shift `past_excitation[160..320]` → `[0..160]`, clear
     `[160..320]`
   - `prev_lag` ← this block's last subframe lag
7. **Per-frame** (after the last subframe):
   - `gain_history`, `gain_threshold`, `gain_counter` ← from the gain
     orchestrator's outputs
   - `postfilter_history`, `postfilter_delay` ← from
     `postfilter_apply` (in place)

The pitch states `pitch_state_block0` and `pitch_state_block1` are
mutated *during* `decode_subframe_lag` calls (they hold the per-block
`prev_lag` for the sub-1 / sub-3 differential decoding) and are not
explicitly committed.

## 4. Reset frame

The codec recognizes an in-band reset frame
(`(frame[17] & 0x0F) == 0x08` — see [architecture.md](architecture.md#33-reset-frame))
and reacts by calling `DecoderState::reset()`, which is
`*self = DecoderState::new()`. After reset the decoder behaves as if it
had just been instantiated. The reset frame itself produces no audio —
`Codec::decode_frame` returns `None`.

## 5. Buffer ownership

A few state fields *look* like scratch but are persisted across calls:

- **`response_shaper_buffer`**: this 35-cell `[i16]` is the dword
  scratch buffer of the LPC analysis. It is written 3 times per frame
  (stability pass, block-0 pass, block-1 pass) and the resulting
  buffer is *intentionally* preserved into the next frame, where it
  becomes the initial buffer of the next stability pass. Replacing
  this with a per-frame zero-init breaks bit exactness.
- **`prev_block_coeffs`**: declared in the struct, but the current
  decoder follows the simpler block-average rule and never reads this
  field. It is retained for compatibility with an alternative
  subframe-LSP-interpolation path (`subframe_lsp_interp`) that the
  reverse-engineering surface area exposed but the runtime decode
  path does not use.

## 6. Subframe / block / frame relationship summary

| Quantity                  | Cadence    | Updated when                                      |
| ------------------------- | ---------- | ------------------------------------------------- |
| `lsf_history`             | per frame  | during LSF decode (before stabilization)          |
| `prev_frame_lsp`          | per frame  | after `lsf_to_lsp`                                |
| `response_shaper_buffer`  | per pass   | after each `build_autocorrelation_with_state`     |
| `prev_lpc_stability_flag` | per frame  | after stability check                             |
| `past_excitation`         | per sub    | inside `pitch_adaptive_codebook` and after mix    |
| `past_excitation` (shift) | per block  | after subframe 1 and subframe 3                   |
| `lpc_synth_history`       | per sub    | after `lpc_synthesis_filter`                      |
| `prev_pitch_gain`         | per sub    | after `gain_orchestrate_codec`                    |
| `prev_lag`                | per block  | after subframe 1 and subframe 3                   |
| `gain_history`            | per sub    | inside `gain_history_update`; final commit per frame |
| `gain_threshold/counter`  | per sub    | inside `gain_tail_decay`; final commit per frame  |
| `postfilter_history/delay`| per half-frame | inside `postfilter_apply`                     |
| `pitch_state_block{0,1}`  | per sub    | inside `decode_subframe_lag`                      |
