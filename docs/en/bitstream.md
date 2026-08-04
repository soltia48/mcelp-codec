# Bit stream

One frame is 18 bytes — 144 bits — of which 139 are payload and the remaining
5 are framing. The payload is fourteen parameter fields totalling 138 bits,
followed by a single suppression flag.

## Field layout

Fields are laid out most significant bit first, with these widths:

```
[8, 12, 8, 16, 7, 5, 16, 7, 8, 16, 7, 5, 16, 7]
```

| field | bits | meaning |
|---|---|---|
| 0 | 8 | line spectrum: predictor mode in bit 7, first-stage index in bits 6–0 |
| 1 | 12 | line spectrum: second stage, two 6-bit halves |
| 2 | 8 | subframe 0 — pitch lag, absolute |
| 3 | 16 | subframe 0 — fixed-codebook index |
| 4 | 7 | subframe 0 — gain index |
| 5 | 5 | subframe 1 — pitch lag, relative to subframe 0 |
| 6, 7 | 16, 7 | subframe 1 — code and gain |
| 8–10 | 8, 16, 7 | subframe 2 — a second absolute lag, code, gain |
| 11–13 | 5, 16, 7 | subframe 3 — relative lag, code, gain |

The two halves of the frame are symmetric: each opens with an absolute pitch
lag and follows it with a delta. `Params::subframe(half, sub)` returns the
three fields belonging to one subframe by name, so no caller outside
`bitstream` needs to know the ordering.

Only the line spectrum is sent once per frame. Everything else is per
subframe, and the fixed-codebook index alone accounts for 64 of the 138
parameter bits.

## Suppression and reset

Bit 138 is the frame-suppression flag: it marks a frame the decoder should
conceal rather than decode. Bit 143 — the last bit of the last byte — carries
a transported copy of the same flag, and canonicalisation folds it down onto
bit 138 before anything else looks at the frame.

The top nibble of the last byte doubles as an in-band reset marker. It is
tested on the *canonicalised* payload, where those bits have already been
cleared, so in this transport the marker can never fire. The check is kept
because the two decoders have to agree bit for bit, including on paths that
are never taken.

## Container

The bundled examples use a plain text container: one frame per line, each line
36 hexadecimal characters encoding the 18 bytes.

```
27ef9c1000701f00070182000e03e000e000
...
```

`bitstream::parse_hex_line` reads one, `to_hex_line` writes one. Neither is
part of the codec proper — the format exists so that a bit stream can be kept
in a text file and diffed.

## Working below the frame

The library exposes the parameter level as well as the byte level:

```rust
use mcelp::bitstream;

let words = bitstream::canonicalize(&frame);   // 18 bytes → 9 words
let params = bitstream::unpack(&words);        // → 14 fields + suppress flag
let again = bitstream::to_bytes(&bitstream::pack(&params));
assert_eq!(again, frame);
```

`Encoder::frame` returns `Params` rather than bytes for the same reason: a
caller studying the coder usually wants the fields, and packing them is a
separate concern.
