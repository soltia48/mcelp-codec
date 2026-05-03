# Excitation: Pitch and Fixed Codebooks

The excitation signal `e[n]` for each 80-sample subframe is the linear
combination of two contributions:

```
e[n] = pitch_gain · v[n] + fcb_gain · c[n]      (sum scaled by mix_excitation)
```

- `v[n]` is the **adaptive codebook** output — a fractional-pitch copy
  drawn from the past excitation buffer (long-term predictor).
- `c[n]` is the **fixed codebook** output — a sparse algebraic codeword
  optionally enhanced with a pitch-periodicity IIR.

This document covers the generation of `v` and `c`. Gain quantization (the
computation of `pitch_gain` and `fcb_gain`) is the subject of
[gain.md](gain.md). Source: [src/pitch.rs](../../src/pitch.rs),
[src/fcb.rs](../../src/fcb.rs), [src/tables/pitch.rs](../../src/tables/pitch.rs),
[src/tables/fcb_short.rs](../../src/tables/fcb_short.rs),
[src/tables/fcb_main.rs](../../src/tables/fcb_main.rs).

## 1. Pitch lag decoding

Each block carries an absolute pitch lag for its first subframe and a
differential lag for its second subframe. The bitstream layout is:

| Subframe | Field | Width | Coding              |
| -------- | ----- | ----: | ------------------- |
| sub 0    | F[2]  | 8     | absolute, block 0   |
| sub 1    | F[5]  | 5     | differential, block 0 |
| sub 2    | F[8]  | 8     | absolute, block 1   |
| sub 3    | F[11] | 5     | differential, block 1 |

The `PitchLagState` struct in `src/pitch.rs` carries `prev_lag` within a
block (sub 0 → sub 1 in block 0, sub 2 → sub 3 in block 1). Block 0's
state is independent of block 1's state.

### 1.1 Absolute decoding

`decode_lag_absolute` (8-bit input):

```
if phase < 197:
    lag     = (phase + 2) / 3 + 19            // integer lag, range 20..85
    sub_lag = phase - lag · 3 + 58            // fractional offset, signed
else:
    lag     = phase - 112                     // 85..143, integer-only
    sub_lag = 0
```

Two regimes are present: small lags (`< 85`) are coded with 1/3-sample
fractional resolution (so that `(lag, sub_lag)` together represents a
fractional pitch in steps of `1/3`); large lags (`≥ 85`) are integer-only.
The integer range covered is 20..143 samples (≈ 56 Hz to 400 Hz at 8 kHz).

### 1.2 Differential decoding

`decode_lag_differential` (5-bit input + previous-subframe lag):

```
anchor = max(prev_lag - 5, 20)
if anchor + 9 > 143: anchor = 134            // upper clamp
q  = (phase + 2) / 3
t1 = q - 1
lag     = anchor + t1                        // anchor-5..anchor+9
sub_lag = phase - 2 - t1 · 3
```

The differential lag spans `[anchor-1, anchor+10]` integer values, again
with 1/3-sample fractional resolution. The anchor is the previous lag
shifted down by 5 (or pinned to 134 if that would put us beyond the high
end) — the same shape that a typical CELP delta-pitch coder uses.

### 1.3 Fractional offset to interpolation tap

`decode_lag_fract(sub_lag)` converts the signed `sub_lag` into:

- a **fract index** `∈ {0, 1, 2}` selecting an interpolation filter
- a **lag adjustment** `lag_adjust ∈ {0, 1}` to be added to the integer lag

```
a = -sub_lag
lag_adjust = (a < 0) ? 1 : 0
a += lag_adjust · 3
return (a, lag_adjust)
```

The end result is that `(integer_lag, sub_lag)` maps to a half-open
fractional pitch and a phase index `0..2` that selects one of three
10-tap interpolation filters.

## 2. Adaptive codebook (fractional-pitch interpolation)

`pitch_adaptive_codebook` synthesises the 80 samples of `v[n]` by
interpolating between past excitation samples with a 10-tap polyphase
filter. The interpolation filters live in
[src/tables/pitch.rs](../../src/tables/pitch.rs):

```
INTERP_FILTER[3][10] in Q15:
    fract=0:  long, slowly-decaying lobes (full-resolution centre filter)
    fract=1:  +1/3-sample shift filter
    fract=2:  +2/3-sample shift filter
```

The forward filter `h` is `INTERP_FILTER[fract]` directly. The reverse
filter `h_rev` is `[h[1..10], 0]` if `fract == 0`, otherwise the *other*
fractional filter (`fract=1 → INTERP_FILTER[2]`,
`fract=2 → INTERP_FILTER[1]`). This produces a symmetric two-sided
convolution centred at the integer-lag tap.

For each output sample `n ∈ 0..80`:

```
acc = Σ_{k=0..9} past[base + n - k] · h[k]
    + Σ_{k=0..9} past[base + n + 1 + k] · h_rev[k]
v[n] = clamp_i16( (acc + (1 << 14)) >> 15 )

past[write_offset + n] = v[n]   ← self-feedback within the same call
```

`base` is `write_offset − effective_lag`, where
`effective_lag = lag + lag_adjust`. The crucial detail is the in-loop
**self-feedback**: each output sample is written back into the past-
excitation buffer at `write_offset + n`, so that subsequent samples of the
same subframe can refer to it. This makes very short pitch lags
(`lag < 80`) behave as a recursive pitch predictor — it generates the
periodic structure even when the lag is shorter than the subframe length.

`write_offset` is `state::PAST_EXCITATION_OUTPUT_BASE + sub_in_block`
(i.e. 160 for the first subframe in a block, 240 for the second). The
past-excitation buffer layout is described in [state.md](state.md).

A bounds-check around the call ensures `base ≥ 9` and `base + 80 + 10`
stays in range; if not, `v` is left at zeros (graceful degradation when
the pitch state is not yet fully populated).

## 3. Fixed codebook

The FCB is **algebraic**: each 16-bit FCB index is decoded into a small
set of unit-magnitude pulse positions, scaled to a fixed pulse magnitude,
and optionally periodicised by an IIR. Two paths are taken depending on
how large the index is:

```
fcb_index in {short path | main path 0 | 1 | 2 | 3 | 4}
            ↑
            └ DISPATCH_THRESHOLDS table selects the path (= "lag class")
```

### 3.1 Dispatch and lag class

`fcb_dispatch_lag_class(fcb_index)` (in [src/fcb.rs](../../src/fcb.rs))
classifies the 16-bit index against the 6-entry `DISPATCH_THRESHOLDS`:

```
DISPATCH_THRESHOLDS = [65432, 49048, 45568, 39168, 32768, 16384]
```

The mapping is:

| Range (unsigned)        | lag_class | Path        |
| ----------------------- | :-------: | ----------- |
| `[0, 16384)`            | 0         | main, k=0   |
| `[16384, 32768)`        | 1         | main, k=1   |
| `[32768, 39168)`        | 2         | main, k=2   |
| `[39168, 45568)`        | 3         | main, k=3   |
| `[45568, 49048)`        | 4         | main, k=4   |
| `[49048, 65432)`        | 5         | **short**   |
| `[65432, ∞)` (clamped)  | 5         | short       |

So the format reserves *6 lag classes* — five for the main algebraic
codebook with progressively narrower bias offsets and one for a tighter
"short" codebook. `clamped_fcb_index` is the index after clamping to the
`THRESHOLDS[0] - 1` upper bound.

### 3.2 Short path (κ = 5)

For lag class 5, `fcb_short_path` decodes the codeword into a pair of
14-bit subindices for two 5-pulse tracks (Track A and Track B):

```
delta = (clamped_index - CODEWORD_BIAS) mod 2^16
track_a = delta >> 7     // 7 bits
track_b = delta & 0x7F   // 7 bits
```

Each track row is a 5-element `i16` array of *signed* pulse positions
encoded as `(sign, |pos| - 1)`:

```
TRACK_A[track_a] = [b0, b1, b2, b3, b4]    // each |b| - 1 is a position 0..79
TRACK_B[track_b] = [b0, b1, b2, b3, b4]
```

`fcb_short_pulse_synthesis` deposits unit pulses of magnitude **8192**
(Q15 ≈ 0.25) at each track's positions. The pulses superpose where the
two tracks land on the same position. The full short codebook therefore
spans `128 · 128 = 16,384` codewords.

If the integer lag is shorter than a full subframe (`lag < 80`),
`pitch_enhance` is then applied (see §3.4).

### 3.3 Main path (κ ∈ 0..5)

For main lag classes 0..4, the codeword decodes into a sequence of 2 or
3 pulses placed at positions chosen from a *position codebook* indexed by
a sequence of small "recursion outputs". The decoding pipeline is:

```
fcb_main_pulse_decode → (bit_count, bit_decomposition[3], recursion_outputs[3])
fcb_main_pulse_template → 80-sample-per-block pulse template (amplitudes)
fcb_main_codebook_synth → final 80-sample c[n]
```

Per-class parameters (`BIT_COUNT`, `BIAS`, `T_MPYA_TABLE`, `T_MPY_TABLE`,
`PULSE_SEEDS`, `PULSE_DATA`, `PULSE_POSITION_CODEBOOK`) live in
[src/tables/fcb_main.rs](../../src/tables/fcb_main.rs). The decoding works
as follows:

1. **Effective index** — subtract a per-class `BIAS` (so that index 0 in
   each lag-class window starts the per-class sub-codebook at 0).
2. **Bit decomposition** — split the low `bit_count` bits of the
   effective index into individual sign bits `bit_decomposition[k] ∈
   {0, 1}`.
3. **Recursion** — the upper bits of the index drive a small recursion
   that, on each step `k`, computes:

   ```
   m1                = (T_MPYA[k] · a_hi) >> 17
   recursion_outputs[k] = a_hi - T_MPY[k] · m1
   a_hi              = m1
   ```

   This expresses the upper bits as a base-`T_MPY[k]` representation —
   it is a fixed-point quotient/remainder cascade selecting one of
   `T_MPY[k]` positions per pulse.
4. **Pulse template** — a 5-row template of 20 samples (each followed by
   60 zeros to fill an 80-sample block) is read from `PULSE_DATA`,
   starting at offset `PULSE_SEEDS[lag_class] · 20`. The number of
   blocks is `bit_count + 1`. The template gives the **amplitudes** of
   the pulses for this lag class (these are not unit pulses but a
   weighted, decaying envelope; see [src/tables/fcb_main.rs](../../src/tables/fcb_main.rs)).
5. **Codebook synthesis** — for each of the `bit_count` pulses,
   `PULSE_POSITION_CODEBOOK[lag_class · 120 + k · 40 + position]` gives
   the *starting position* in the 80-sample subframe, where `position`
   is `recursion_outputs[bit_count - 1 - k]`. The pulse template at
   `pulse_offset = k · 80` is then either added (if
   `bit_decomposition[k] == 1`) or subtracted (if `== 0`) into `c[]`,
   starting at that position. Saturation is applied per sample.

The codebook size for each main lag class is the gap between adjacent
`DISPATCH_THRESHOLDS` (i.e. the number of distinct `clamped_fcb_index`
values that route to that class):

| `lag_class` | `BIT_COUNT` | Pulse clusters | Codebook size |
| :--------: | :--------: | :----: | -----: |
| 0          | 3          | 3      | 16384  |
| 1          | 3          | 3      | 16384  |
| 2          | 2          | 2      |  6400  |
| 3          | 2          | 2      |  6400  |
| 4          | 2          | 2      |  3480  |

After codebook synthesis, `pitch_enhance` is again applied if `lag < 80`.

### 3.4 Pitch periodicity enhancement (`pitch_enhance`)

For both the short and main paths, when the integer pitch lag is smaller
than a full subframe (`lag < 80`), the FCB output `c[n]` is post-processed
by `pitch_enhance`. The role is to **inject the period of the pitch into
the FCB excitation**, so that voiced frames retain a periodic FCB
component even though the FCB codewords are themselves aperiodic.

The implementation is a single-tap recursive filter that re-uses the same
10-tap pitch-interpolation filter as the adaptive codebook:

```
lag_internal = lag - 10 + lag_adjust
copy buf[0..lag_internal] into a scratch line (offset 21)
for n in lag_internal..80:
    acc = Σ_k scratch[centre - k] · fwd[k]
        + Σ_k scratch[centre + 1 + k] · rev[k]
    new = sat16( (buf[n] << 16 + (acc · 2 · gain) >> 14 + 0x8000) >> 16 )
    buf[n]    = new
    scratch[…] = new       ← the new value is fed back into the scratch line
```

`gain` is the **pitch enhance gain** computed by
`compute_pitch_enhance_gain`, which dispatches on the lag class:

- Lag classes `{0, 3, 4}`: gain is always **enabled** (`= 16384 = Q15
  0.5`).
- Lag classes `{1, 2, 5}`: gain is the output of `update_synth_control`
  (see [synthesis.md](synthesis.md#synth-control)) which inspects the
  current pitch lag, the previous pitch lag, the first reflection
  coefficient $k_1$ from LPC analysis, and the previous pitch gain to
  decide whether to enable it.

When `gain = 0` the IIR becomes a pass-through (the recursion contributes
nothing) and the FCB output is left as the unmodified algebraic codeword.

## 4. Excitation mixing

After both `v` and `c` are formed, [`mix_excitation`][synth] in
`src/synth.rs` produces:

```
e[n] = sat16(((fcb_gain · c[n]) << 3 + (pitch_gain · v[n]) << 2 + 32768) >> 16)
```

— the 32-bit accumulator holds `8 · fcb_gain · c[n] + 4 · pitch_gain · v[n]`
plus a Q15 round bias, then is right-shifted to Q15 with i32 saturation.
The factors of 8 and 4 are fixed scaling factors that compensate the
Q-format choices of the gain pipeline.

The mixed excitation `e[n]` is written back to the past-excitation buffer
at `write_offset + n` (overwriting the temporary `v[n]` self-feedback
copy). It is also passed forward to the LPC synthesis filter (see
[synthesis.md](synthesis.md)).

[synth]: ../../src/synth.rs
