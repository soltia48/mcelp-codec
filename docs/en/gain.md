# Gain Quantization

The pitch gain $g_p$ and the FCB gain $g_c$ for each subframe are jointly
encoded into a **7-bit phase word** carried by F[4], F[7], F[10], F[13]
of the bitstream. Recovering them on the decode side is non-trivial: it
involves a phase-table lookup, a predictive refresh that walks several
fixed tables, an exponential history update of the prediction state, and a
threshold-driven tail-decay limiter.

This document walks through that pipeline. Source:
[src/gain.rs](../../src/gain.rs), [src/tables/gain.rs](../../src/tables/gain.rs).
The orchestrator entry point is `gain_orchestrate_codec`.

## 1. Inputs and outputs

### Inputs

```rust
struct GainOrchestrateCodecInput {
    phase_word:    i16,    // 7-bit field from the bitstream
    threshold_in:  i16,    // tail-decay threshold (frame-persistent)
    counter_in:    i16,    // tail-decay counter   (frame-persistent)
    candidate:     [i16; 80], // the FCB output c[n] for this subframe
    history:       [i16; 5],  // past-gain history (frame-persistent)
}
```

### Outputs

```rust
struct GainOrchestrateOutput {
    pitch_gain:     i16,   // g_p, fed into mix_excitation as pitch_gain
    fcb_gain:       i16,   // g_c, fed into mix_excitation as fcb_gain
    history_out:    [i16; 5],
    threshold_out:  i16,
    counter_out:    i16,
}
```

The `phase_word` is shared by both gains — that is, the encoder transmits
a single 7-bit index per subframe that selects one of 128 jointly
quantized `(g_p, g_c)` pairs. The other inputs come from per-frame state
that is updated subframe-by-subframe (see [state.md](state.md)).

## 2. Phase setup

`gain_phase_setup(phase_word, …)` decomposes the phase word into two
4-bit halves and uses each half to look up two cells (an "upper" and
"lower" pair) of `PHASE_TABLE`:

```
phase_first_idx  = (phase_word >> 4) << 1     // 0..14
phase_second_idx = ((phase_word & 0xF) << 1) + 16   // 16..46

phase_first_upper = PHASE_TABLE[phase_first_idx]
phase_first_lower = PHASE_TABLE[phase_first_idx + 1]
phase_second_upper= PHASE_TABLE[phase_second_idx]
phase_second_lower= PHASE_TABLE[phase_second_idx + 1]
```

`PHASE_TABLE` is a 48-cell quantization table from
[src/tables/gain.rs](../../src/tables/gain.rs). The **upper 3 bits** of
`phase_word` (i.e. `phase_word >> 4`, range 0..7) select a pair of cells
in `PHASE_TABLE[0..16]`, and the **lower 4 bits** (range 0..15) select a
pair of cells in `PHASE_TABLE[16..48]`. Conceptually the two halves carry
quantized values for $g_c$ and $g_p$ respectively, but they are not used
as-is — they go through the predictive refresh below.

`gain_phase_setup` then computes:

- `update_base = ((phase_first_upper + phase_first_lower) << 16) >> 16`
  — initial $g_p$ candidate.
- `saved_ah    = ((phase_second_upper + phase_second_lower) << 16) >> 16`
  — bookkeeping value passed into the history setup.

It also derives two i40 accumulators (`acc_main`, `acc_scratch`) seeded
with the phase-word bits to feed the predictive refresh.

## 3. Predictive refresh

`gain_predictive_refresh` is the dominant computation. It produces
`history_shift_init` (an integer shift count) and `history_sample_init`
(a Q15 sample) which together define the **excitation reference level**
the gain quantizer will normalise against. The 10 phases of the function
mirror the original DSP routine; the high-level structure is:

```
1. Energy of the FCB candidate:
       E = 10·2^17 + 2·Σ c[n]^2
   (the 10·2^17 is a small-energy floor)
2. Primary normalisation against PREDICT_NORMALIZE_TABLE
       — pick a lookup index from the normalised energy,
         interpolate between two table cells (Q15 frac).
3. Linear correction (subtract 19488<<19, add 18432<<20).
4. 4-tap inner product with the past-gain history:
       acc_main = Σ_{k=0..3} COEFF_TABLE[k] · history[k]
   where COEFF_TABLE = [5571, 4096, 2703, 1311] is a 4-pole AR predictor.
5. Combine + 5439-multiplier mix → history_shift_init.
6. Secondary normalisation against PREDICT_SECONDARY_TABLE
       (the same shape as primary, with different coefficients).
7. history_sample_init = hi16 of the post-secondary accumulator.
```

The two normalisation steps are implemented by `gain_normalize_primary`
and `gain_normalize_secondary`. Both are **piecewise-linear** approximations
of (different) nonlinear functions, with a 33-cell lookup table indexed
by the normalised mantissa of the i40 accumulator and a Q15 fractional
weight to interpolate between cells. The structure is:

```
acc          ← sft acc by 5 (Q?? → Q??)
t_exp        ← exp_acc(acc)              (Q40 normalize)
acc          ← norm_acc_with_t(acc, t_exp)
lookup_index ← upper bits of acc - 16384·1.0
interp_frac  ← lower bits of acc
result       ← table[lookup_index]
             - interp_frac · (table[lookup_index] - table[lookup_index+1])
```

`PREDICT_NORMALIZE_TABLE` and `PREDICT_SECONDARY_TABLE` (33 cells each)
are tabulated logarithm-like curves: the primary table is approximately
$\log_2$ of a 16-bit number scaled to Q15, and the secondary table is a
companion that reshapes the result for the second normalisation. They
live in [src/tables/gain.rs](../../src/tables/gain.rs).

`gain_orchestrate_codec` performs one *subtle* additional optimisation
over the symmetric `gain_orchestrate`: rather than calling the
normalisation functions and then propagating their results, it
**predicts** the lookup index from the input and pre-fetches
`PREDICT_NORMALIZE_TABLE[idx]` and `PREDICT_NORMALIZE_TABLE[idx+1]`
*outside* the normalisation. The function is then called with these two
table cells precomputed, and the same is done for the secondary
normalisation by replaying phases 2–9 of the predictive refresh once to
peek at the secondary input. This pattern is preserved in the Rust code
to remain bit-exact with the reverse-engineered behaviour.

## 4. History setup and update

`gain_history_setup` produces the **initial gain** that will become the
candidate for $g_c$ before tail decay:

```
restored_acc = saved_ah ∥ predict_acc_main_low_bits  (40-bit reconstruct)
shift_count  = history_shift_init - 7                (no-CMPT backoff)
acc          = sample · acc32_16(restored_acc) · 2
acc          = shift_acc40(acc, shift_count)
acc          = shift_acc40(acc, -1)
acc          = sat32(acc)
initial_gain = hi16(acc)
```

`gain_history_update` then refreshes the 5-cell past-gain history. It
runs another `gain_normalize_primary` (this time using a third lookup
prediction over `PREDICT_NORMALIZE_TABLE`), shifts the result into the
post-update Q15 accumulator, and computes a new history entry as a
constant-multiplier-mix of the high and low halves of the saturated
accumulator with `0x6054 = 24660` as the multiplier:

```
new_gain = ((scratch_lo·24660·2) >> 16 + scratch_hi·24660·2) >> 16
history_out = [new_gain, history[0], history[1], history[2]]
```

After the refresh, `history_out[0..4]` is `[new_gain, prev[0], prev[1],
prev[2]]` and `history_out[4] = history[4]` is preserved as a constant
slot — it carries a frame-startup constant (initialised to 0 in
`DecoderState::new`).

## 5. Tail decay

`gain_tail_decay` is a soft limiter that cuts the gains once the
suppress-state threshold has fired. Every subframe it does:

```
threshold_in -= 8                  (TAIL_THRESHOLD_STEP)
if threshold_in >= 0:
    counter ← TAIL_COUNTER_RELOAD = 4
threshold_out = 0

if counter ≤ 0:
    return (no decay)              ← unscaled (g_p, g_c) pass through

counter -= 1
t_value = -counter + 3             (TAIL_EXP_BIAS)
scale   = MIN_SYNTH_STATE · 2^t_value     (= 3277 << t_value)
update_base ·= scale          ← scaled g_p
initial_gain ·= scale          ← scaled g_c
return (decay applied, counter_out)
```

The effect is that once the threshold trips, the gains are
multiplicatively scaled down by a power of 2 that increases with
`t_value` (= more decay later in the tail), counted down for at most
`TAIL_COUNTER_RELOAD` = 4 subframes. The constant `MIN_SYNTH_STATE =
3277` is the minimum-state floor used in several places in the codec
(see [synthesis.md](synthesis.md)).

In the steady state (no suppress), `threshold_in` is 0 across all
subframes and `counter_out` stays 0, so the tail-decay block is a pass-
through and `(g_p, g_c) = (update_base, initial_gain)`.

## 6. Wiring

The four gain results emitted per subframe are wired in
[src/lib.rs](../../src/lib.rs) as:

```
pitch_gain   = tail.update_base_out   → mix_excitation
fcb_gain     = tail.initial_gain_out  → mix_excitation
threshold    = tail.threshold_out     → next subframe
counter      = tail.counter_out       → next subframe
history_out[5] → next subframe (frame-persistent at frame end)
prev_pitch_gain = pitch_gain          → next subframe's compute_pitch_enhance_gain
```

Notice the **commit ordering**: `prev_pitch_gain` is updated *after*
gain orchestration but *before* the next subframe's
`compute_pitch_enhance_gain` is called, which is essential because the
pitch-enhance gain depends on the previous pitch gain via the synth-
control hysteresis (see [synthesis.md](synthesis.md#synth-control)).

## 7. Initial state

At codec init `gain_history` is set to `[-17254, -17254, -17254, -17254,
0]` and `gain_threshold = gain_counter = 0`. The non-zero history
seeds the AR predictor with a small negative bias that approximates the
quiescent pitch-gain envelope — without it, the first few decoded frames
would mispredict severely until the history fills with real values.
