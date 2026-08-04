# Short-term prediction

The short-term predictor models the spectral envelope — roughly, the shape the
vocal tract imposes. It is a 10th-order all-pole filter

```
          1                       1
H(z) = ------ = ---------------------------------
        A(z)     1 + a₁z⁻¹ + a₂z⁻² + … + a₁₀z⁻¹⁰
```

and one set of coefficients is sent per frame.

## Analysis

The encoder autocorrelates a **400-sample window** — longer than the 320-sample
frame, so it reaches back into the previous one. The window is asymmetric and
tapered; the overlap is what keeps the filter from jumping between frames.

The autocorrelations are then *lag-windowed*: each lag is multiplied by a
factor slightly below one, decreasing with lag. This widens the spectral peaks
a little and is what stops the Levinson-Durbin recursion from producing a
filter so sharp that it rings.

Levinson-Durbin turns the eleven autocorrelations into the eleven predictor
coefficients, and produces the ten reflection coefficients on the way. The
first two reflection coefficients are kept: the perceptual weighting filter's
shape is chosen from them.

Note that the analysis runs on the **noise-suppressed** signal, not on the
input. See [noise-suppression.md](noise-suppression.md).

## Line spectral frequencies

Predictor coefficients quantise badly — small errors can make the filter
unstable — so they are converted to *line spectral frequencies* first. `A(z)`
is split into a symmetric and an antisymmetric polynomial,

```
F₁(z) = A(z) + z⁻¹¹A(z⁻¹)      F₂(z) = A(z) − z⁻¹¹A(z⁻¹)
```

whose roots all lie on the unit circle and interleave. Those root angles are
the LSFs. Interleaving is the useful property: as long as the quantised values
stay in ascending order with a minimum spacing, the filter that comes back is
stable, whatever the quantisation error.

The roots are found by evaluating each polynomial on a grid of points around
the unit circle, watching for sign changes, and bisecting four times inside
each bracket.

Three representations appear in the code, and mixing them up is the most
common way to misread it:

| name | meaning | format |
|---|---|---|
| LSF | line spectral *frequencies*, in radians | Q13, π = 25736 |
| LSP | line spectral *pairs*, i.e. cos(LSF) | Q15 |
| LPC | direct-form coefficients | Q12, `a[0] = 4096` |

## The quantiser

Two stages, predictive:

```
LSF ──(subtract prediction)──► residual ──► stage 1: 128 entries × 10  ─► 7 bits
                                            stage 2: 64 entries, applied
                                                     to each half of the
                                                     vector separately  ─► 6+6 bits
```

The prediction is a third-order moving average over the three previous
frames' quantised residuals. Two predictor sets exist and the encoder picks
whichever gives less distortion, sending the choice as one bit — the mode bit
in field 0. Splitting the second stage into a low half and a high half means
64 entries buy what 4096 would otherwise cost.

Distortion is measured with per-coefficient weights: an LSF sitting close to
its neighbours describes a sharp formant and matters more, so it is weighted
up.

After dequantisation the decoder enforces the ordering the representation
depends on — a minimum gap between adjacent LSFs, plus clamps at both ends of
the range — and only then converts back to LPC.

## Interpolation

One line spectrum per frame is not enough for four subframes, so the rest are
interpolated. The interpolation is done in the LSP domain, where a linear
blend of two stable vectors is itself stable:

- the first half-frame is weighted against the **midpoint** of the previous
  frame's spectrum and this one's;
- the second half-frame uses this frame's spectrum whole;
- within a half-frame, the first subframe again uses a midpoint and the second
  the endpoint.

There is an escape hatch. If the two spectra differ sharply — a transient —
interpolating between them would describe a filter that matches neither. The
encoder tests the reflection coefficients of the decoded filter and can
disable interpolation for the frame, in which case both half-frames use the
current spectrum. The decoder makes the same test on the same values and so
reaches the same decision without being told.
