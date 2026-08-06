# One frame, in detail

Everything else in these documents describes the coder in general. This one
takes a single real frame out of the bundled example and shows what is actually
in it.

The frame is **number 11 of `examples/male.mcelp`** — the loudest of the 75,
which makes it a clearly voiced one and therefore a good specimen. All values
below are what the decoder actually computes; you can reproduce them from the
public API.

## The 18 bytes

```
8bc0e96cbdf270a68e6723e8b575c7a6bb80
```

Unpacked into the fourteen fields:

| field | value | hex | meaning |
|---|---:|---|---|
| 0 | 139 | `008b` | LSF: predictor mode 1, first-stage index 11 |
| 1 | 3086 | `0c0e` | LSF: second stage, low half 14, high half 48 |
| 2 | 150 | `0096` | subframe 0 — pitch lag |
| 3 | −13345 | `cbdf` | subframe 0 — codebook index |
| 4 | 19 | `0013` | subframe 0 — gain index |
| 5 | 16 | `0010` | subframe 1 — lag delta |
| 6 | −22898 | `a68e` | subframe 1 — codebook index |
| 7 | 51 | `0033` | subframe 1 — gain index |
| 8 | 145 | `0091` | subframe 2 — pitch lag |
| 9 | −2982 | `f45a` | subframe 2 — codebook index |
| 10 | 93 | `005d` | subframe 2 — gain index |
| 11 | 14 | `000e` | subframe 3 — lag delta |
| 12 | 15669 | `3d35` | subframe 3 — codebook index |
| 13 | 110 | `006e` | subframe 3 — gain index |

The suppression flag is clear, so this is an ordinary frame.

Note that the codebook indices print as negative numbers. They are 16-bit
fields held in an `i16`; read as unsigned they are `0xcbdf`, `0xa68e`, `0xf45a`
and `0x3d35`, and it is the unsigned value that selects a class.

## The spectrum

Fields 0 and 1 give the quantiser indices:

```
mode = 1,  stage1 = 11,  stage2_low = 14,  stage2_high = 48
```

Dequantising and converting to line spectral frequencies gives, in Q13 radians
and in hertz:

| | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Q13 | 2813 | 3549 | 7166 | 9497 | 10438 | 12433 | 15565 | 16477 | 19464 | 20408 |
| Hz | 437 | 552 | 1114 | 1476 | 1622 | 1932 | 2419 | 2561 | 3025 | 3172 |

They are strictly ascending, as the representation requires. What makes them
readable is that **a closely spaced pair marks a resonance**:

| pair | spacing | implies a formant near |
|---|---:|---|
| 437 / 552 | 114 Hz | 494 Hz |
| 2419 / 2561 | 142 Hz | 2490 Hz |
| 3025 / 3172 | 147 Hz | 3099 Hz |

while the widely spaced pairs (362 Hz, 310 Hz) describe the flatter region
between them. A first formant near 500 Hz with a second in the 1300–1800 Hz
region is a typical male back vowel.

This is the property the quantiser's weighting exploits: the tightly spaced
lines are the ones that must be coded accurately, because moving them moves a
formant.

## The pitch

| subframe | field | decoded lag | as frequency |
|---|---:|---|---|
| 0 | 150, absolute | 69 + ⅓ | 115 Hz |
| 1 | 16, relative | 69 − ⅓ | 116 Hz |
| 2 | 145, absolute | 68 − ⅓ | 118 Hz |
| 3 | 14, relative | 67 + 0 | 119 Hz |

A lag of about 69 samples at 8 kHz is a fundamental near 116 Hz — a male
voice — and it drifts upward by two samples across the 40 ms of the frame,
which is an ordinary rate of pitch change in running speech.

Notice how cheap the second and fourth lags are. Each is 5 bits rather than 8,
and each lands within a sample or two of the lag before it, which is exactly
the assumption the differential coding is built on.

## The innovations

| subframe | index | class | structure |
|---|---|---:|---|
| 0 | `0xcbdf` | 5 | twin five-pulse |
| 1 | `0xa68e` | 3 | two pulses, 40 × 40 positions |
| 2 | `0xf45a` | 5 | twin five-pulse |
| 3 | `0x3d35` | 0 | three pulses, 8 × 16 × 16 |

Three different structures inside one frame, which is the point of having six
of them: the search picks per subframe, and a voiced frame is not uniformly
pulse-like from start to finish.

Across both bundled examples — 600 subframes — the classes are used like this:

| class | share |
|---:|---:|
| 0 | 33.0 % |
| 1 | 15.3 % |
| 2 | 8.8 % |
| 3 | 7.7 % |
| 4 | 4.3 % |
| 5 | 30.8 % |

The two extremes dominate: class 0 (three pulses on the finest track layout)
and class 5 (the dense twin codebook) together take almost two thirds. The
middle classes exist for what falls between.

## The gains

| subframe | index | first stage | second stage |
|---|---:|---:|---:|
| 0 | 19 | 1 | 3 |
| 1 | 51 | 3 | 3 |
| 2 | 93 | 5 | 13 |
| 3 | 110 | 6 | 14 |

Each index splits into a 3-bit row of the first stage and a 4-bit row of the
second; the two rows are summed component-wise to give the pitch gain and the
code-gain correction factor. The rising indices across the frame track the
signal getting louder — this is the loudest frame of the file, and it is still
rising at its end.

## What the decoder does with all this

For each of the four subframes, in order:

1. interpolate the line spectrum to get this subframe's `1/A(z)`;
2. read the past excitation at the lag to get the adaptive contribution;
3. rebuild the innovation from the 16-bit index;
4. dequantise the pitch gain, predict the code gain and correct it;
5. sum the two scaled contributions into the excitation;
6. run `1/A(z)` to get 80 samples of speech;
7. postfilter and emit.

Eighteen bytes in; 320 μ-law bytes out.

## Reproducing this

```rust
use mcelp::{bitstream, lsp};

let text = std::fs::read_to_string("examples/male.mcelp")?;
let frame = text.lines().nth(11).and_then(bitstream::parse_hex_line).unwrap();
let params = bitstream::unpack(&bitstream::canonicalize(&frame));

println!("{:?}", params.field);
println!("{:?}", lsp::LsfIndices::unpack(params.field[0], params.field[1]));
for half in 0..2 {
    for sub in 0..2 {
        println!("{:?}", params.subframe(half, sub));
    }
}
```

The lag fields need `pitch::decode_absolute` for the first subframe of each
half and `pitch::decode_relative` for the second; the class of a codebook
index is the largest entry of `tables::FCB_CLASS_BASE` it clears.
