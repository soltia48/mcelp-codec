# Architecture

## The three time units

Everything in the codec is organised around a nesting of three lengths:

```
frame        320 samples, 40 ms          one transport packet
  half-frame 160 samples                 one line spectrum, one analysis window
    subframe  80 samples                 one excitation search
```

Why three rather than two: the line spectrum is expensive to send, so it is
sent once per frame and *interpolated* to give each half-frame and each
subframe its own filter. The excitation is cheap to make but has to track the
signal closely, so it is searched four times per frame. The half-frame is the
unit in between — the encoder's spectral analysis and its noise suppression run
on 160 new samples at a time.

## Encoder

```
μ-law in (320)
  │
  ├─ preprocess ──────── expand to linear, high-pass
  │
  ├─ analysis ────────── 256-sample window, 128-point FFT, magnitude spectrum
  │    │
  │    ├─ bands ──────── fold 129 bins into 20 critical bands, convert to dB
  │    ├─ noise ──────── two long-term noise floors, updated on silent frames
  │    ├─ vad ────────── voice-activity score from band statistics
  │    └─ weights ────── per-band suppression depth
  │
  ├─ shaping ─────────── apply suppression, inverse transform, stitch spans
  │                      ↓ this cleaned signal is what everything below sees
  │
  ├─ lpc ─────────────── autocorrelation over 400 samples, Levinson-Durbin,
  │                      convert to line spectral frequencies
  ├─ lsf_weight ──────── quantise them  ─────────────────► fields 0, 1
  │
  └─ per subframe (×4):
       weighting ─────── build the perceptual weighting filter
       pitch_search ──── open-loop lag, then closed-loop refinement ► lag field
       mode ──────────── decide how this subframe is to be coded
       pulses ────────── search the fixed codebook ──────► code field
       gain ──────────── search the gain codebook ───────► gain field
       excitation ────── build the excitation, update the filter memories
```

The subframe loop is described in detail in
[encoder-loop.md](encoder-loop.md), and the weighting filter it is built
around in [weighting.md](weighting.md).

The encoder synthesises as it goes. After choosing an excitation it runs it
through the same synthesis filter the decoder will use, so that the next
subframe's search starts from exactly the state the decoder will be in. This
is what "analysis by synthesis" means in practice, and it is why an encoder
bug shows up as drift rather than as a single bad frame.

## Decoder

```
18 bytes in
  │
  ├─ bitstream ───────── unpack 14 fields
  ├─ lsp ─────────────── dequantise the line spectrum, interpolate per subframe
  │
  ├─ per subframe (×4):
  │    pitch ─────────── adaptive codebook: resample past excitation at the lag
  │    codebook ──────── fixed codebook: rebuild the innovation from the index
  │    gain ─────────── dequantise both gains, predict the code gain
  │    decoder ───────── sum the two contributions into the excitation
  │    synth ─────────── run 1/A(z) to get speech
  │
  ├─ postfilter ──────── per subframe: short-term + long-term postfilter,
  │                      tilt compensation, gain matching
  └─ output filter ───── final reconstruction filter, then compand to μ-law
```

## Where the time goes

Everything expensive in CELP is in the search, and the measurements bear that
out. Encoding and decoding the 3 seconds of `examples/male.orig.ulaw`, on one
core of a desktop CPU:

| | per 3 s of audio | real time |
|---|---:|---:|
| Encode | 37 ms | 81× |
| Decode | 2.5 ms | 1184× |

The encoder is about **15 times** the decoder's cost. That asymmetry is
inherent to the design rather than an artefact of this implementation: the
decoder evaluates one excitation, while the encoder evaluates several hundred
candidates per subframe and synthesises the winner as well.

Within the encoder, the fixed-codebook search dominates — it is the stage with
the largest candidate set and the one whose index costs the most bits. The
spectral front end is the next largest, since it runs an FFT and twenty bands
of statistics twice per frame.

Both figures are far above real time, so neither side is a practical
constraint on a modern machine; the ratio is what matters when reasoning about
where an optimisation would pay.

## What survives a frame

Both halves of the codec are stateful, and the state is what makes a stream
undecodable from the middle. The decoder carries:

| state | why |
|---|---|
| Line spectrum predictor memory | the quantiser is predictive: three past frames feed the next |
| Previous half-frame's line spectrum | one end of the interpolation |
| Synthesis filter memory | ten samples of `1/A(z)` |
| Excitation history | 154 samples, so a lag of up to 143 can reach back |
| Synthesised speech, residual history | the postfilter's inputs |
| Gain predictor memory | four past frames of code-gain energy |
| Postfilter memories, AGC state | continuity across subframes |
| Concealment lag, RNG state | for erased frames |

The encoder carries the decoder's state (it must, to stay in step) plus its
own: the noise floors, the voice-activity counters, the suppression scale, the
analysis window's history, and the open-loop peak track.

## Module map

| module | role |
|---|---|
| `bitstream` | frame packing and unpacking |
| `preprocess`, `ulaw` | companding and input conditioning |
| `analysis`, `bands`, `noise`, `vad`, `frontend`, `weights`, `shaping` | encoder front end |
| `lpc`, `lsf_weight` | prediction analysis and line spectrum quantisation |
| `lsp` | line spectrum representations, dequantisation, interpolation |
| `weighting`, `weight_lpc` | perceptual weighting |
| `pitch_search`, `convolve` | pitch search |
| `pitch` | adaptive codebook |
| `pulses`, `shape_pair`, `mode` | fixed-codebook search |
| `codebook` | fixed-codebook decoding |
| `gain` | gain codebook, both directions |
| `excitation` | assembling a subframe's excitation |
| `synth` | the synthesis filter |
| `postfilter`, `ltp` | decoder postfiltering |
| `encoder`, `decoder` | the frame loops that tie the stages together |
| `fixed` | the arithmetic primitives |
| `tables` | the constant tables |

Only `Encoder` and `Decoder` are meant to be called from outside. The stage
modules are public so that each can be tested on its own against reference
data, which is what the test suite does.
