# LPC: LSF → LSP → LPC

Mitsubishi CELP transmits the spectral envelope as a 10-th order **Line
Spectral Frequency (LSF)** vector that is split-VQ encoded with inter-frame
prediction. On the decoder side the frame's LSF is converted into LSPs and
then expanded into per-subframe LPC coefficients through a 3-pass autocorrelation
pipeline that the source calls `response_shaper`. This document describes
that chain end-to-end.

Source files:

- [src/lsf.rs](../src/lsf.rs) — LSF decoding and LSP conversion
- [src/lpc_analysis.rs](../src/lpc_analysis.rs) — autocorrelation + Levinson +
  stability check
- [src/tables/lsf.rs](../src/tables/lsf.rs) — VQ codebooks and predictor
- [src/tables/lsp.rs](../src/tables/lsp.rs) — cosine LUT

## 1. LSF decoding

The frame's first two bitstream fields (F[0], F[1]) carry the LSF index:

```rust
struct LsfIndices {
    mode:  bool,  // F[0] bit 7        — selects predictor / smoothing row
    seed:  u8,    // F[0] bits 0..6    — stage-1 codebook index (7 bits)
    upper: u8,    // F[1] bits 6..11   — stage-2 codebook index for k=0..4
    lower: u8,    // F[1] bits 0..5    — stage-2 codebook index for k=5..9
}
```

That is **8 + 12 = 20 bits** for the 10-D LSF, of which 1 bit selects the
predictor mode and 19 bits index the split-VQ. The names `upper` and
`lower` refer to the **bit position within F[1]**, not the LSP index
range.

### 1.1 Split-VQ codebook combination

`combine_codebook` (see [`combine_codebook`][lsf]) sums a 7-bit stage-1
vector with two 6-bit stage-2 vectors:

```
scratch[0..5]  = STAGE1_CODEBOOK[seed][0..5]  + STAGE2_CODEBOOK[upper][0..5]
scratch[5..10] = STAGE1_CODEBOOK[seed][5..10] + STAGE2_CODEBOOK[lower][5..10]
```

Tables are stored in Q15 in [src/tables/lsf.rs](../src/tables/lsf.rs)
(`STAGE1_CODEBOOK` is 128×10, `STAGE2_CODEBOOK` is 64×10).

### 1.2 Minimum-gap enforcement

After combining, `enforce_min_gap` is called twice — once with `min_gap = 10`
and once with `min_gap = 5` — to nudge each pair `(scratch[i], scratch[i+1])`
apart so that `scratch[i+1] − scratch[i] ≥ −min_gap`. The split between the
two passes (10 then 5) reproduces a two-stage smoothing schedule used by the
reference DSP.

### 1.3 Predictive combination

Inter-frame prediction is computed by `predictive_combine` against a 3-frame
LSF history (`LsfHistory`) that holds the **scratch** vectors of the past
3 frames (oldest at slot 2). For each `i ∈ 0..10`:

```
out[i] = ( scratch[i]    · SMOOTH[mode][i]
         + history[0][i] · PREDICTOR[mode][i]
         + history[1][i] · PREDICTOR[mode][i + 10]
         + history[2][i] · PREDICTOR[mode][i + 20]
         ) · 2 / 65536
```

(`× 2 / 65536` reproduces the DSP's `frct_mul=2` followed by a `>>16`.)

The two predictor modes (`mode ∈ {0, 1}`) provide alternative
predictor / smoother rows; the selection bit is the MSB of F[0]. Mode 0
typically corresponds to a less aggressively predictive setting and mode 1
to a more aggressive one.

After predictive combination the *raw scratch* (not the predicted output)
is spliced into the history with `LsfHistory::splice_new_vector`.

### 1.4 Stabilization

The predicted output is finalized by `stabilize_and_finalize`:

1. **Bubble-1 pass** swaps adjacent out-of-order entries.
2. **Lower clamp** `lsf[0] ≥ 41` (≈ 0.00125 in Q15).
3. **Forward minimum-gap = 321** to guarantee LSP separability.
4. **Upper clamp** `lsf[9] ≤ 25682` (≈ 0.7838 in Q15).

After this pass the LSF is monotonically increasing and bounded. The LSP
conversion that follows depends on these invariants.

## 2. LSF → LSP

`lsf_to_lsp` (see [src/lsf.rs](../src/lsf.rs)) implements the discrete
cosine transform that maps a normalized LSF (Q15, on the range `[0, 0.5)`)
to the corresponding LSP (= `cos(π · lsf)`):

```
scaled       = lsf · 20861 · 2     // i64
index        = clamp(scaled >> 24, 0, 63)
fraction_bits= (scaled >> 16) & 0xFF
lsp = COS_LUT_VALUE[index] + (COS_LUT_SLOPE[index] · fraction_bits) >> 12
```

The constant 20861 is `64 · 2^15 / π ≈ 20861` so that
`scaled >> 24` is the integer index into a 64-entry cosine LUT, while the
remaining 8 fraction bits (Q8) interpolate linearly between table cells
using a precomputed slope. Tables `COS_LUT_VALUE` and `COS_LUT_SLOPE` are
in [src/tables/lsp.rs](../src/tables/lsp.rs).

The decoder also retains a `prev_frame_lsp` so that `block_average_lsp`
(component-wise `(prev + curr) >> 1`) can produce the **block-0 work
coefficients**, which are the average of the previous frame's LSP and the
current frame's LSP. The block-1 work coefficients are simply the current
frame's LSP itself.

## 3. The 3-pass `response_shaper` pipeline

The interesting part of LPC analysis on this codec is that the synthesis
filter coefficients are **not** computed directly from the current frame's
LSP. Instead, `lib.rs` calls `build_autocorrelation_with_state` three times
per frame, with three different "windows", and the resulting coefficients
are split into per-subframe LPCs. Each pass mutates a shared 35-cell
scratch buffer (`response_shaper_buffer` in `DecoderState`), so the
sequence is order-dependent.

### 3.1 Pass A — stability check

Input window: `lsp_q15` followed by 4 zeros (`stability_seed_window =
[lsp_q15..., 0, 0, 0, 0]`).

```
(stability_autocorr, post_buf) = build_autocorrelation_with_state(
    prev_frame_lsp,            // "work" coefficients
    stability_seed_window,     // "seed" window
    response_shaper_buffer,
);
response_shaper_buffer ← post_buf
```

The decoder takes the **second half** of `stability_autocorr` (cells
`12..22`) as a putative LPC vector and runs 10 Levinson-Durbin steps on it
(`levinson_recursion`) to get the 10 reflection coefficients
$k_1, \dots, k_{10}$. Those are passed to `is_lpc_stable`, which returns
`true` iff $|k_1| < 1$ and a tail metric — formed by accumulating
`d14b_metric_residual(k_i)` across `i ∈ 2..10` as a 32-bit dword pair —
also clears the `STABILITY_REF` threshold.

The decision rule for the **block-0 LPC override** is:

```
override_block_0 = stable_now && prev_lpc_stability_flag == 0
```

i.e. an override fires only on the *rising edge* of stability. When
override fires, block-0 work coefficients are replaced by the raw
`lsp_q15` instead of the per-block average. `prev_lpc_stability_flag` is
then updated. This rising-edge gate is a hysteresis that prevents jitter
near the stability boundary.

### 3.2 Pass B — block-0 autocorrelation

Input window: `[block_0_coeffs_eff..., lsp_q15[0..4]]`.

```
(autocorr_block0, post_buf) = build_autocorrelation_with_state(
    prev_frame_lsp,
    block0_seed_window,
    response_shaper_buffer,
);
response_shaper_buffer ← post_buf
sub0_lpc = autocorr_block0[1..11]
sub1_lpc = autocorr_block0[12..22]
```

`build_autocorrelation_with_state` internally:

1. Forms a 14-cell **first_autocorr** = `[(work + seed[0..10])>>1, seed[0..4]]`.
2. Calls `response_shaper(first_autocorr, …)` → `mirror_output[0..11]`
   that is the 11-tap LPC for the first subframe of this pass.
3. Calls `response_shaper(seed_window, …)` → `mirror_output[0..11]` for
   the second subframe.
4. Concatenates them into a 22-cell result.

A single `response_shaper` is itself an analysis stage: two
`levinson_step` accumulations that fill the front and back halves of the
35-cell dword buffer, a **merge phase** that adds the symmetric pairs and
subtracts the antisymmetric pairs (sums on `[10, 8, 6, 4, 2]` and
differences on `[22, 20, 18, 16, 14]`, repeated 5 times), and a **mirror
phase** that emits the 11-cell LPC by reading sums and differences of
those halves (with `mirror[0] = 4096 = Q12 1.0` always). The mathematical
intuition is that the front half encodes a *symmetric* (P-polynomial-like)
half and the back half an *antisymmetric* (Q-polynomial-like) half, and
the mirror phase is a fast LSP→LPC reconstruction. The 11-cell output is
in Q12.

In addition, `levinson_recursion` is called on each subframe's 10-cell
LPC tail (= cells `1..11` and `12..22` of `autocorr_block0`) to extract
$k_1$ as `clamp_input` for that subframe — the gain pipeline uses it as
the `clamp_input` for the synth-control decision.

### 3.3 Pass C — block-1 autocorrelation

Input window: `stability_seed_window` (the same `[lsp_q15..., 0, 0, 0, 0]`
used in pass A, now with `block_0_coeffs_eff` as the work-coefficient
input).

```
(autocorr_block1, post_buf) = build_autocorrelation_with_state(
    block_0_coeffs_eff,
    stability_seed_window,
    response_shaper_buffer,
);
response_shaper_buffer ← post_buf
sub2_lpc = autocorr_block1[1..11]
sub3_lpc = autocorr_block1[12..22]
```

Conceptually this is the symmetric construction for block 1, with the
roles of "work" (= the block-0 coefficients) and "seed" (= the current
LSP) chosen so that block 1's first subframe is the average of block-0's
output and the current frame's LSP, and block 1's second subframe is the
current frame's LSP itself.

### 3.4 Why three passes?

The role of each pass:

| Pass | Window combination                      | Role                                            |
| ---- | --------------------------------------- | ----------------------------------------------- |
| A    | `prev_lsp ↔ [lsp                     0]`| Decide block-0 override (stability check)       |
| B    | `prev_lsp ↔ [block_0_coeffs   lsp[0..4]]`| LPC for sub 0 and sub 1                         |
| C    | `block_0_coeffs ↔ [lsp                0]`| LPC for sub 2 and sub 3                         |

The three passes share the 35-cell scratch buffer in a deterministic
forward order. Replacing any with a stand-alone analysis breaks bit
exactness because the dword scratch carries information between passes —
this is a property of the original DSP firmware that the Rust
implementation reproduces faithfully.

### 3.5 The reflection coefficients $k_i$ as `clamp_input`

For each subframe `i ∈ 0..4`, the value
`clamp_inputs[i] = levinson_recursion(sub_i_lpc)[0]` (i.e. the first
reflection coefficient $k_1$ obtained from a 10-step Levinson recursion on
the 10-cell autocorrelation tail) is passed into the synth-control
decision (see [synthesis.md](synthesis.md#synth-control)). $k_1$ is a
classical voicing indicator (`k_1 → -1` for highly voiced frames), so
the decoder uses it to decide whether to enable the pitch-periodicity
enhancement IIR for that subframe.

A separate utility `lpc_to_first_reflection` is also provided in
[src/lsf.rs](../src/lsf.rs); it implements the inverse Levinson (Schur
recursion) on a Q12 LPC vector. It is currently unused in the runtime
decode path but is kept for symmetry with encoder-side analysis and for
unit testing.

## 4. The Levinson-Durbin recursion

`levinson_recursion` implements the classical Levinson-Durbin recursion,
specialised for 10 steps and using the i40-accumulator semantics shared
with the rest of the codec. The recursion produces, in place, the
reflection coefficients $k_1, \dots, k_{10}$ given a 10-cell
autocorrelation tail.

The inner step computes

```
E_i      = E_{i-1} · (1 - k_i^2)
k_i      = -⟨a_{i-1}, R_i⟩ / E_{i-1}      (forward step)
a_i[m]   = a_{i-1}[m] - k_i · a_{i-1}[i-m]  for m ∈ 1..i  (reverse step)
```

with two specifics:

- **`q15_reciprocal`** computes `1 / E_{i-1}` in Q15 by an iterative
  `subc`-based division (15 steps of the DSP `subc` instruction
  implemented as `subc_step`, followed by a 32-bit-mask scrub and 40-bit
  sign extension). This reproduces the DSP firmware's reciprocal estimate
  bit-for-bit.
- The shift amount `asm` is derived from the upper 5 bits of the
  reciprocal's normalisation exponent, then sign-extended from 5 bits.

The Q15 reflection coefficients returned by `levinson_recursion` are
consumed both by the stability check (`is_lpc_stable`) and by the per-
subframe `clamp_input` lookup in the gain pipeline.

## 5. Per-subframe LPC summary

After the three passes complete, the decoder holds **four 10-cell Q12 LPC
vectors** — `sub0_lpc`, `sub1_lpc`, `sub2_lpc`, `sub3_lpc` — together with
**four Q15 first-reflection coefficients** — `clamp_inputs[0..4]`. From
this point the LPC analysis is finished; the synthesis pipeline does not
further transform these coefficients.

[lsf]: ../src/lsf.rs
