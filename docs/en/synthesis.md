# Synthesis, Postfilter, and μ-law Output

The remaining steps after excitation generation and gain quantization
are: (1) the synth-control decision used by both pitch enhancement and
the synthesis filter, (2) the LPC synthesis filter itself, (3) the
half-frame postfilter, and (4) the linear-to-μ-law conversion. This
document covers all four.

Source: [src/synth.rs](../../src/synth.rs),
[src/postfilter.rs](../../src/postfilter.rs),
[src/ulaw.rs](../../src/ulaw.rs),
[src/tables/postfilter.rs](../../src/tables/postfilter.rs),
[src/tables/ulaw.rs](../../src/tables/ulaw.rs).

## 1. <a name="synth-control"></a> Synth control (pitch-enhance gating)

`compute_pitch_enhance_gain` decides whether the pitch-periodicity
enhancement IIR (see [excitation.md](excitation.md#34-pitch-periodicity-enhancement-pitch_enhance))
should be enabled for the current subframe. Its inputs are:

- `fcb_index` — the FCB codeword for this subframe
- `current_lag` — this subframe's pitch lag
- `clamp_input` — the first reflection coefficient $k_1$ for this
  subframe (output of `levinson_recursion`; see
  [lpc.md](lpc.md#35-the-reflection-coefficients-k_i-as-clamp_input))
- `prev_lag` — previous half-frame's selected lag
- `prev_pitch_gain` — previous subframe's $g_p$

The dispatch on lag class is hard-wired:

| Lag class | Behavior                                          |
| --------- | ------------------------------------------------- |
| 0, 3, 4   | always **enable** (gain = 16384 = Q15 0.5)        |
| 1, 2, 5   | follow the `update_synth_control` decision        |

`update_synth_control` itself implements a 4-branch hysteresis. It first
clamps `prev_pitch_gain` to `[MIN_SYNTH_STATE, MAX_SYNTH_STATE] =
[3277, 13017]` to use as the **previous control state**, then decides:

```
if clamp_input is too large (≈ |k_1| > 0.5):
    return DISABLE                    (gain = 0)
elif prev_control_state ≤ SYNTH_CONTROL_GATE = 13107:
    return prev_control_state         (carry over)
elif current_lag ≤ prev_lag · 24576/32768:
    return prev_control_state         (slow lag drift → keep)
elif current_lag <  prev_lag · 24576/32768 · 2:
    return ENABLE                     (gain = 16384)
else:
    return prev_control_state         (carry over)
```

In words: the IIR is **disabled** when the LPC indicates a non-voiced
frame ($k_1$ too large), and otherwise it is enabled when the pitch lag
has *changed quickly* (jumped by between 1× and 2× the
24576/32768-scaled previous lag) — that is, the format expects the
encoder to drive a step change in pitch when the speaker switches from
unvoiced to voiced. For lag classes 0, 3, 4 — which correspond to FCB
ranges that the format reserves for strongly-voiced signals — the IIR is
unconditionally enabled.

The constants in this routine are:

```
SYNTH_CONTROL_ENABLE  = 16384   (Q15 0.5 — the active gain)
SYNTH_CONTROL_DISABLE = 0
SYNTH_CONTROL_GATE    = 13107   (≈ 0.4)
PREV_LAG_BLEND_GAIN   = 24576   (= 0.75)
MIN_SYNTH_STATE       = 3277
MAX_SYNTH_STATE       = 13017
```

## 2. Excitation mixing

`mix_excitation(v, c, g_p, g_c)` produces a Q15 excitation
`e[n] ∈ i16` for `n ∈ 0..80`:

```
acc = g_c · c[n] · 8 + g_p · v[n] · 4 + 32768       (i64)
acc = sat32(acc)
e[n] = sat16(acc >> 16)
```

The factors `· 8` and `· 4` are fixed scale corrections that align the
internal Q-formats of $g_c$ and $g_p$ to the synthesis filter's
expected Q15 input. The `+ 32768` is the round-to-nearest bias for the
final `>> 16`.

## 3. LPC synthesis filter ($1/A(z)$)

`lpc_synthesis_filter(excitation, lpc_coeffs, history)` (in
[src/synth.rs](../../src/synth.rs)) is a straight all-pole IIR:

```
buf[0..10]  = history (last 10 samples of the previous subframe's output)
for n in 0..80:
    mac = Σ_{k=0..9} lpc_coeffs[k] · buf[9 + n - k]   (Q12 × Q15 = Q27)
    acc = (excitation[n] << 16) - mac · 16 + 32768   (i64)
    acc = sat32(acc)
    buf[10 + n] = sat16(acc >> 16)
output[0..80]    = buf[10..90]
new_history[0..10] = buf[80..90]
```

The Q-format choices:

- `excitation[n]` is Q15. `excitation[n] << 16` lifts it into Q31.
- `lpc_coeffs[k]` is Q12. `mac` is Q27. `mac * 16` lifts it into Q31.
- The accumulator is i32-saturated, then right-shifted to Q15 with
  rounding.

Because the LPC coefficients come directly from the
`build_autocorrelation_with_state` mirror output (see [lpc.md](lpc.md)),
**no separate LSP→LPC conversion** is run on the synthesis path —
`response_shaper.mirror_phase` produces the Q12 LPC directly.

`lpc_coeffs[0]` (the implicit `a_0 = 1.0`) is **not** included in the
10-element vector passed to the synthesis filter; the filter loops over
$a_1 \dots a_{10}$ and adds the $a_0 \cdot e[n]$ term separately as
`excitation[n] << 16`.

The `history` carried out of one subframe becomes the `history` for the
next subframe. The decoder stores it as `lpc_synth_history` in
`DecoderState`.

## 4. Past-excitation buffer management

The 320-cell `past_excitation` buffer holds the last two blocks' worth of
excitation samples:

```
[0..160]   ← previous block's e[n]
[160..240] ← current block's sub a (phase 0)  e[n]
[240..320] ← current block's sub b (phase 80) e[n]
```

Each subframe writes its 80-sample excitation into the current block's
slot, then `pitch_adaptive_codebook` of the *next* subframe (or the next
block) reads from this buffer using a back-offset of `effective_lag`.

After every block — that is, after the second subframe (`i % 2 == 1`)
— [src/lib.rs](../../src/lib.rs) shifts `[160..320]` down to `[0..160]`,
clears `[160..320]`, and commits `prev_lag = lag_int` so the next block's
synth-control decision has the correct reference.

## 5. Postfilter

After all four subframes are synthesised the decoder holds a 320-sample
Q15 buffer `synth_pcm`. The postfilter is then applied per **half-frame**
(160 samples = block) by `postfilter_apply` in
[src/postfilter.rs](../../src/postfilter.rs).

The postfilter is a **forward + reverse IIR** that combines formant
emphasis with a tilt correction and a final output-gain trim. Per
sample:

```
input_acc         = sample << 16
b1 = input_step(input_acc, fwd_hist, INPUT_COEFF)  (4-tap MAC + 3-cell shift register)
b2 = feedback_correction(b1, delay_fwd, FEEDBACK_COEFF)  (3 IIR feedback sections)
shift_delay_line(delay_fwd, b2)
rounded = round_acc(b2)

b3 = input_step(rounded, rev_hist, INPUT_COEFF)
b4 = feedback_correction(b3, delay_rev, FEEDBACK_COEFF)
shift_delay_line(delay_rev, b4)

pf = output_sample(b4, OUTPUT_GAIN)
sample = (pf >> 16) & 0xffff
```

The fixed coefficients are:

```
INPUT_COEFF    = [7472, -22381,  22381, -7472]   (4-tap input MAC)
FEEDBACK_COEFF = [-23033, 21666,  -6815]         (3-tap feedback IIR)
OUTPUT_GAIN    = 19661                           (output trim ≈ 0.6 in Q15)
```

The postfilter holds two pieces of state:

- `postfilter_history[6]` — 3-cell forward + 3-cell reverse history
  shift registers
- `postfilter_delay[12]` — 6-cell forward + 6-cell reverse dword pair
  delay line (each `feedback_correction` reads three dword pairs)

The forward stage uses the first half of each, the reverse stage the
second half. Both halves are persisted in `DecoderState` across half-
frames (and, by extension, across frames).

`round_acc` and `output_sample` are micro-helpers that reproduce the DSP
firmware's saturation order (`<< 1`, `+ 0x8000`, `sat32`); they round
each intermediate step to the same 16-bit grid the original DSP would
have used.

## 6. μ-law conversion

The final 320 Q15 samples are converted to μ-law bytes by
`linear_i16_to_ulaw` (in [src/ulaw.rs](../../src/ulaw.rs)), implementing
the standard G.711 μ-law encoder:

```
m0 = sample if sample >= 0 else !sample        (one's-complement absolute value)
m  = min((m0 >> 2) + 33, 8191)                 (bias + 13-bit saturation)

seg_bits  = m >> 6                              (segment in 0..127)
seg_shift = 0 if seg_bits == 0 else floor(log2(seg_bits)) + 1
mantissa  = (m >> (seg_shift + 1)) & 0x0F
segment   = 7 - seg_shift

ulaw = (segment << 4) | (15 - mantissa)
if sample >= 0: ulaw |= 0x80
```

The inverse direction is a 256-entry LUT (`ULAW_TO_LIN`). It is
provided for symmetry but is **not used in the decode path** — the codec
only emits μ-law, never decodes it back to linear PCM as part of its own
operation.

A subtlety: the encoder collapses negative zero to positive zero
(`linear_i16_to_ulaw(0) == 0xFF`), which means the table-roundtrip is
*not* a strict inverse; bytes 127 and 255 both map to linear 0 and then
back to 0xFF. This is the standard G.711 negative-zero folding and
matches how the rest of the telephony stack interprets the format.

## 7. Why μ-law output?

The codec is intended to integrate with G.711 / μ-law-based telephony
endpoints, where the input/output of an audio path is already μ-law. By
emitting μ-law directly, the decoder avoids one extra linearization step
when feeding into a μ-law sink. The `Codec::decode_frame` API mirrors
this by returning `[u8; 320]` of μ-law rather than `[i16; 320]` of
linear PCM.
