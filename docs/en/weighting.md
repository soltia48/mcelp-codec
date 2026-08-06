# Perceptual weighting

Every search in the encoder — pitch, fixed codebook, gain — asks the same
question: which candidate sounds closest to the input? The weighting filter is
what turns "sounds closest" into a number, and it is the single most important
reason a CELP coder at 3.6 kbit/s is listenable at all.

## Why weight at all

A coder that minimised plain squared error would spread its quantisation noise
evenly across the spectrum. That is the wrong thing to do, because the ear does
not hear noise evenly. Noise sitting under a formant is masked by the formant;
the same amount of noise in a spectral valley is plainly audible.

So the error is measured through a filter that *de-emphasises* the formant
regions before measuring. The coder is then free to leave more error where the
signal is loud, which is exactly where it will not be heard.

## The filter

```
        A(z/γ₁)
W(z) = ---------           γ₁ > γ₂
        A(z/γ₂)
```

Both halves come from the same predictor `A(z)` by **bandwidth expansion** —
multiplying coefficient `aₖ` by `γᵏ`:

```rust
out[i] = a[i] · γⁱ            // for i = 1 … 10
```

Multiplying by a factor below one pulls the filter's poles in towards the
origin, which broadens its resonances without moving them. The numerator,
with γ₁ close to one, keeps the formant structure nearly intact; the
denominator, with a much smaller γ₂, is a smeared version of the same
envelope. Their ratio is a filter whose response follows the formants but
with the peaks flattened — precisely the shape the masking argument asks for.

## Choosing the two factors

The factors are not constant. They are chosen **once per half-frame**, giving
one pair for each of its two subframes.

| state | γ₁ (numerator) | γ₂ (denominator) |
|---|---|---|
| sharp | 32113 ≈ 0.980 | from the line spectrum, 0.400 … 0.700 |
| flat | 30802 ≈ 0.940 | 19661 = 0.600, fixed |

In the **sharp** state the denominator is derived from how tightly the line
spectrum is bunched. The smallest gap between adjacent LSFs is taken, put
through a straight line, and clamped to the range above. Two lines close
together mean a sharp resonance, and a sharp resonance is one the weighting
filter has to follow more closely — so a small gap gives a larger γ₂ and a
narrower weighting.

In the **flat** state both factors are fixed. A spectrum without strong
resonances has nothing to follow, and adapting to noise would only make the
filter jitter between subframes.

## The sharp/flat decision

The state is decided from the **first two reflection coefficients** of the
half-frame's predictor. These two numbers summarise the gross tilt of the
spectrum: the first is essentially the normalised one-sample autocorrelation,
the second the two-sample one.

They are compared not directly but as **log area ratios**, an approximately
logarithmic transform that keeps the comparison meaningful as a coefficient
approaches unity. The implementation approximates the curve with three
straight segments above a knee, below which the coefficient is already close
enough to its own logarithm to be used as it stands.

The decision has **hysteresis** — two different threshold pairs, one to enter
the flat state and one to leave it:

```
entering flat:  area₀ ≤ −3113  and  area₁ ≥  881
leaving  flat:  area₀ <  −3564  and  area₁ >  1331
```

Hysteresis matters here for the same reason it matters in the voice-activity
detector. A spectrum sitting right at a single threshold would flip the
weighting filter's shape every half-frame, and a filter that changes shape
faster than the signal does introduces its own modulation.

One more refinement: the two subframes of a half-frame are judged slightly
differently. The first is judged on log areas **smoothed** against the
previous half-frame's, the second on this half-frame's alone. The filter
therefore lags slightly behind a changing spectrum at the start of a
half-frame and catches up by its end, rather than stepping.

## Where it is used

The weighting filter appears in three places, and they must agree:

1. **The target.** The signal the searches try to match is the weighted
   version of what remains to be coded — see [encoder-loop.md](encoder-loop.md).
2. **The impulse response.** Every candidate is compared after passing through
   the weighted synthesis filter, so the search precomputes that filter's
   impulse response once per subframe.
3. **The quantiser weights.** The line spectrum quantiser measures its own
   distortion under weights derived from the same filter.

The decoder never sees any of this. Weighting is entirely an encoder-side
measurement device; it shapes *which* parameters get sent, not what the
decoder does with them.
