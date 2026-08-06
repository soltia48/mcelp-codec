# Glossary

Speech coding has a dense vocabulary, and several terms mean something narrower
here than in general use. This lists them as they are used in this codebase.

### Adaptive codebook

The recent excitation, read at a delay of one pitch period. Not a fixed table —
it re-forms every subframe from the coder's own output, which is why it is
called adaptive. Also called the *long-term predictor*. → [pitch.md](pitch.md)

### Analysis by synthesis

Choosing coder parameters by synthesising the result of each candidate and
keeping the one closest to the input, rather than by measuring the input
directly. → [encoder-loop.md](encoder-loop.md)

### Bandwidth expansion

Multiplying predictor coefficient `aₖ` by `γᵏ` for some `γ < 1`. Pulls the
filter's poles towards the origin, broadening its resonances without moving
them. Used to build the weighting filter and the postfilter.

### Class

One of the six structures the 16-bit fixed-codebook index can select. Which
class an index belongs to is decided by which range it falls in.
→ [fixed-codebook.md](fixed-codebook.md)

### Closed loop

A search that synthesises candidates and measures the result. Expensive, and
therefore run only over a short range that an open-loop search has bracketed.

### Excitation

The signal that drives the synthesis filter — what the LPC model cannot
predict. In this codec it is the sum of a scaled adaptive contribution and a
scaled fixed-codebook contribution.

### Fixed codebook

The structured set the *innovation* is chosen from: sparse pulse patterns, one
per subframe, selected by a 16-bit index.
→ [fixed-codebook.md](fixed-codebook.md)

### Formant

A resonance of the vocal tract, appearing as a peak in the spectral envelope.
Speech intelligibility depends mostly on the first two or three.

### Frame / half-frame / subframe

320 / 160 / 80 samples. The frame is the transport unit, the half-frame the
analysis unit, the subframe the excitation unit.
→ [architecture.md](architecture.md)

### Innovation

The fixed codebook's contribution to the excitation — the "new" part, as
opposed to the periodic part the adaptive codebook supplies.

### Lag

The delay, in samples, at which the adaptive codebook reads the past
excitation. Equal to one pitch period for voiced speech. Range 20–143 here.

### LPC — linear predictive coding

Modelling a sample as a weighted sum of the preceding ones. The weights are the
*predictor coefficients*; the filter `1/A(z)` built from them is the *synthesis
filter*. → [lpc.md](lpc.md)

### LSF / LSP

Two equivalent alternative representations of the predictor. LSF are *line
spectral frequencies*, angles in radians; LSP are *line spectral pairs*, their
cosines. Both quantise and interpolate far better than the raw coefficients
because stability reduces to keeping them in ascending order.

### Open loop

A cheap search run directly on the signal, without synthesis. Used to bracket
the closed-loop pitch search.

### Perceptual weighting

Measuring coding error through a filter that de-emphasises the formant regions,
so the coder leaves its error where the ear will not hear it.
→ [weighting.md](weighting.md)

### Pitch sharpening

Feeding the innovation through the fractional-delay pitch predictor before use,
reinforcing periodicity the pulses alone would leave ragged. A decoder-side
operation.

### Postfilter

A decoder-side filter that deepens the spectral valleys, hiding quantisation
noise where the signal is too weak to mask it. → [synthesis.md](synthesis.md)

### Q notation

Fixed-point scaling. A Q*n* value is an integer to be read as a fraction scaled
by 2ⁿ. → [fixed-point.md](fixed-point.md)

### Reflection coefficient

A by-product of the Levinson-Durbin recursion, one per prediction order. The
first two summarise the gross spectral tilt and are what the weighting filter's
shape decision is made from.

### Residual

What is left after the inverse filter `A(z)` has removed everything the
short-term predictor can explain. The excitation is the coder's approximation
of it.

### Suppression

Two unrelated meanings, unfortunately, and both appear in the code:

- **Noise suppression** — the encoder's spectral cleaning front end.
  → [noise-suppression.md](noise-suppression.md)
- **Frame suppression** — the transport flag marking a frame the decoder should
  conceal rather than decode. → [bitstream.md](bitstream.md)

### Target

The signal a search tries to match: the perceptually weighted input with
everything the already-coded subframes account for subtracted out.
→ [encoder-loop.md](encoder-loop.md)

### Track

A subset of the positions within a subframe that one pulse of a multi-pulse
class may occupy. Interleaving tracks is what lets a few pulses cover the
subframe with a small index.

### VAD — voice activity detection

Deciding whether a frame is speech or background. Here it gates the noise-floor
updates and selects the suppression depth.
→ [noise-suppression.md](noise-suppression.md)
