# Mitsubishi CELP — Technical Documentation

This directory contains a technical description of the Mitsubishi CELP speech
codec (三菱CELP方式音声コーデック) as implemented by this crate. The
implementation was produced by reverse-engineering an existing decoder; this
documentation is therefore a description of the format and algorithm as
inferred from that reverse-engineered behavior, not a publication of an
official specification.

## At a glance

| Property              | Value                                                   |
| --------------------- | ------------------------------------------------------- |
| Sample rate           | 8 kHz, mono                                             |
| Output companding     | μ-law (G.711)                                           |
| Frame length          | 40 ms = 320 samples                                     |
| Compressed frame size | 18 bytes (144 bits, 138 bits payload + 1 control flag)  |
| Bit rate              | 18 B / 40 ms ≈ **3.6 kbit/s** (3.45 kbit/s of payload)  |
| LPC order             | 10                                                      |
| Subframes             | 4 × 80 samples (= 10 ms each), grouped into 2 blocks    |
| LSF quantization      | 2-stage split-VQ with inter-frame prediction (3 frames) |
| Adaptive codebook     | Fractional pitch (1/3-sample resolution), 10-tap interp |
| Fixed codebook        | Algebraic, 6 lag classes — 5 main paths + 1 short path  |
| Gain quantization     | 7-bit phase index + 4-tap predictive history            |
| Postfilter            | Forward + reverse formant/tilt IIR per half-frame       |

## Reading order

1. [architecture.md](architecture.md) — top-level decode pipeline, frame
   layout, the 14 control fields, reset frames.
2. [lpc.md](lpc.md) — the LSF → LSP → LPC chain, including the per-frame
   `response_shaper` autocorrelation pipeline and the stability check that
   gates the block‑0 LPC override.
3. [excitation.md](excitation.md) — pitch lag decoding, the fractional
   adaptive codebook, the algebraic fixed codebook (short / main paths),
   and the pitch‑periodicity enhancement IIR.
4. [gain.md](gain.md) — the gain‑quantization pipeline: phase setup,
   predictive refresh, history update, tail decay.
5. [synthesis.md](synthesis.md) — excitation mixing, the 1/A(z) synthesis
   filter, the half‑frame postfilter, and the final μ-law conversion.
6. [state.md](state.md) — the persistent decoder state and how it is
   carried across frames, blocks, and subframes.

## Source map

| Topic           | Module                                       |
| --------------- | -------------------------------------------- |
| Bitstream       | [src/bitstream.rs](../src/bitstream.rs)      |
| LSF / LSP       | [src/lsf.rs](../src/lsf.rs)                  |
| LPC analysis    | [src/lpc_analysis.rs](../src/lpc_analysis.rs)|
| Pitch           | [src/pitch.rs](../src/pitch.rs)              |
| Fixed codebook  | [src/fcb.rs](../src/fcb.rs)                  |
| Gain            | [src/gain.rs](../src/gain.rs)                |
| Synthesis       | [src/synth.rs](../src/synth.rs)              |
| Postfilter     | [src/postfilter.rs](../src/postfilter.rs)    |
| μ-law           | [src/ulaw.rs](../src/ulaw.rs)                |
| Fixed‑point ops | [src/arith.rs](../src/arith.rs)              |
| Decoder state   | [src/state.rs](../src/state.rs)              |
| Top level       | [src/lib.rs](../src/lib.rs)                  |

## Conventions

Throughout these documents:

- **Q-format.** `Qn` denotes a fixed-point format with `n` fractional bits.
  Most spectral state is kept in **Q15** (signed, range ≈ [-1, +1)), LPC
  coefficients are emitted in **Q12**, and the inner accumulator paths use
  **Q31 / Q40** (i.e. 32- and 40-bit signed accumulators implemented over
  `i64`). The fixed-point primitives live in [src/arith.rs](../src/arith.rs).
- **40-bit accumulator (i40).** The original DSP is a 40-bit-MAC architecture.
  Several arithmetic helpers (`shift_acc40`, `to_i40`, `exp_acc`,
  `norm_acc_with_t`, `acc32_16`) reproduce the i40 wrap/sign-extend
  semantics. The codec is bit-exact with the reverse-engineered reference
  precisely because of these.
- **Endianness.** The 144-bit frame is read MSB-first
  (see [bitstream](architecture.md#bitstream)).
- **Subframe / block.** A *subframe* is 80 samples. Two consecutive
  subframes form a *block* (= 160-sample half-frame). Each frame contains
  two blocks.
