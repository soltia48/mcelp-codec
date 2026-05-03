# Decoder Architecture and Bitstream

## 1. Frame timing

The codec runs at **8 kHz, mono**. Each compressed frame is **18 bytes** and
expands to **320 μ-law samples = 40 ms** of audio. Internally the 320-sample
frame is partitioned hierarchically:

```
frame   : 320 samples (40 ms)            ← unit of compressed I/O
 └─ block: 160 samples (20 ms)           ← half-frame, postfilter unit, LPC override unit
     └─ subframe: 80 samples (10 ms)     ← excitation / synthesis unit
```

A frame contains **2 blocks**, each containing **2 subframes** — 4 subframes
total. Several pipeline stages are organized "per-block" (LPC analysis,
pitch state) while excitation and synthesis are "per-subframe".

## 2. Top-level decode flow

The entry point is [`Codec::decode_frame`](../../src/lib.rs) in
[src/lib.rs](../../src/lib.rs). The decoder runs the following stages in order
on each 18-byte frame:

```
┌─────────────────────────────────────────────────────────────────────┐
│  18-byte input                                                      │
│     │                                                               │
│     ▼                                                               │
│  is_reset_frame? ──yes──► reset state, emit nothing                 │
│     │ no                                                            │
│     ▼                                                               │
│  canonicalize_payload          (move bit 143 → 138, clear 139..143) │
│     │                                                               │
│     ▼                                                               │
│  ControlFrame::unpack          (14 fields + suppress flag)          │
│     │                                                               │
│     ▼                                                               │
│  ┌── Stage 1: LSF decoding ───────────────────────────────┐         │
│  │  LsfIndices.from_control                               │         │
│  │  combine_codebook (split-VQ stage1+stage2)             │         │
│  │  enforce_min_gap × 2                                   │         │
│  │  predictive_combine (3-frame history)                  │         │
│  │  history.splice_new_vector                             │         │
│  │  stabilize_and_finalize (sort/clamp/min-gap)           │         │
│  │  lsf_to_lsp (cosine LUT)                               │         │
│  └────────────────────────────────────────────────────────┘         │
│     │                                                               │
│     ▼                                                               │
│  ┌── Stage 2: LPC analysis (3 response_shaper passes) ───┐          │
│  │  pass A: stability-check on lsp[0..10] | 0000         │          │
│  │           → reflection coeff k1..10                   │          │
│  │           → is_lpc_stable (decides block-0 override)  │          │
│  │  pass B: block-0 autocorrelation                      │          │
│  │           → sub0_lpc, sub1_lpc                        │          │
│  │  pass C: block-1 autocorrelation                      │          │
│  │           → sub2_lpc, sub3_lpc                        │          │
│  │  Each pass carries the 35-cell scratch buffer forward │          │
│  └───────────────────────────────────────────────────────┘          │
│     │                                                               │
│     ▼                                                               │
│  ┌── Stage 3: per-subframe loop (i = 0..4) ──────────────┐          │
│  │  decode_subframe_lag    (block 0/1 absolute or diff)  │          │
│  │  pitch_adaptive_codebook → v[80]                      │          │
│  │  fcb_dispatch_lag_class                               │          │
│  │  compute_pitch_enhance_gain                           │          │
│  │  fcb_short_path or fcb_main_path → c[80]              │          │
│  │  gain_orchestrate_codec → (g_p, g_c, history')        │          │
│  │  mix_excitation: e = 8·g_c·c + 4·g_p·v                │          │
│  │  lpc_synthesis_filter (1/A(z)) → s_hat[80]            │          │
│  │  past_excitation ← e                                  │          │
│  │  if i odd: shift past_excitation; commit prev_lag     │          │
│  └───────────────────────────────────────────────────────┘          │
│     │                                                               │
│     ▼                                                               │
│  postfilter_apply  × 2 half-frames                                  │
│     │                                                               │
│     ▼                                                               │
│  linear_i16_to_ulaw × 320                                           │
│     │                                                               │
│     ▼                                                               │
│  320-byte μ-law output                                              │
└─────────────────────────────────────────────────────────────────────┘
```

The remainder of this document covers the bitstream (Stage 0). The four
stages that follow are the subjects of [lpc.md](lpc.md),
[excitation.md](excitation.md), [gain.md](gain.md), and
[synthesis.md](synthesis.md).

## 3. Bitstream

### 3.1 Frame layout

Each frame is 18 bytes = 144 bits, read MSB-first. The payload is laid out
as 14 contiguous fields followed by a small control region:

```
bit  0                        137  138    139            143
     │  14 fields, MSB-first    │  S │ unused (4 bits) │  M │
     └──────────────────────────┴────┴─────────────────┴────┘
       138 bits of payload      sup.  (cleared on RX)  flag-src

S = bit 138 = "suppress flag" (canonicalized destination)
M = bit 143 = "frame flag source" (transport-level)
```

The 14 payload fields and their bit widths (see [`FIELD_WIDTHS`][1]):

| Idx | Width | Field                                  | Notes                  |
| --: | ----: | -------------------------------------- | ---------------------- |
| 0   | 8     | LSF mode (1 bit) + LSF stage-1 seed (7) | `mode = bit 7`        |
| 1   | 12    | LSF stage-2 upper (6) + lower (6)      | split-VQ index pair    |
| 2   | 8     | Pitch lag absolute, sub 0              |                        |
| 3   | 16    | FCB index, sub 0                       | algebraic codeword     |
| 4   | 7     | Gain phase word, sub 0                 | gain quantizer index   |
| 5   | 5     | Pitch lag differential, sub 1          | relative to sub 0      |
| 6   | 16    | FCB index, sub 1                       |                        |
| 7   | 7     | Gain phase word, sub 1                 |                        |
| 8   | 8     | Pitch lag absolute, sub 2              |                        |
| 9   | 16    | FCB index, sub 2                       |                        |
| 10  | 7     | Gain phase word, sub 2                 |                        |
| 11  | 5     | Pitch lag differential, sub 3          | relative to sub 2      |
| 12  | 16    | FCB index, sub 3                       |                        |
| 13  | 7     | Gain phase word, sub 3                 |                        |

These widths sum to 138. Symmetry across the two blocks is exact: F[2..7]
parameterises sub 0/1 (block 0); F[8..13] parameterises sub 2/3 (block 1).
LSF (F[0..1]) is sent once per frame; both blocks interpolate it locally.

### 3.2 Suppress flag and "flag-source" bit

Two bits in the control region carry transport-level state:

- **bit 143** ("flag-source", `M`) is the bit produced by the encoder's
  framer.
- **bit 138** ("suppress", `S`) is the canonicalized destination of the
  flag inside the decoder's payload model.

The receiver canonicalises the frame (see
[`canonicalize_payload`][1]) by:

1. copying bits `[0..139)` unchanged,
2. copying bit 143 → bit 138,
3. clearing bits `[139..143)`.

After canonicalisation the decoder reads `suppress` from bit 138, which is
exposed as `ControlFrame.suppress`. The current crate decodes the
non-suppress (`= 0`) path — the field is parsed and kept available to the
caller, but the suppress=1 alternative-gain path is not exercised by the
decoder yet.

### 3.3 Reset frame

A frame whose final byte ends in the low nibble `0x8` — i.e.
`(frame[17] & 0x0F) == 0x08` — is treated as an in-band **reset frame**.
This is detected by [`is_reset_frame`][1]. On reset the decoder reinitializes
all persistent state (LSF history, past-excitation buffer, LPC synthesis
history, pitch state, gain history, postfilter state, etc. — see
[state.md](state.md)) and emits no audio for that frame.

Because a reset is detected before canonicalisation, the reset code coexists
with the otherwise legal payload range — it is the *frame-as-marker* convention
of the format, not a payload value.

### 3.4 The `.mcelp` text container

The bundled `.mcelp` files (and the `decode` binary) use a simple, portable
text container: one frame per line, each line a 36-character hexadecimal
encoding of the 18-byte frame. The CLI uses [`parse_frame_hex`][1] to convert.
This is purely a convenience for command-line experiments; the codec library
itself operates on raw 18-byte slices.

[1]: ../../src/bitstream.rs
