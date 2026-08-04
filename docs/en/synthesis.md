# The decoder

## Synthesis

Once a subframe's excitation is assembled, reconstruction is one filter:

```
speech[n] = excitation[n] − Σ aₖ · speech[n−k]      k = 1 … 10
```

`1/A(z)`, run 80 samples at a time, carrying ten samples of memory between
subframes. The predictor is Q12 and the signal Q15, so the accumulated
prediction is shifted to line the two formats up before it is added.

That is the whole of the reconstruction. Everything else the decoder does is
either assembling the excitation — described in [pitch.md](pitch.md),
[fixed-codebook.md](fixed-codebook.md) and [gain.md](gain.md) — or cleaning up
after the fact.

## Postfiltering

At 3.6 kbit/s the quantisation noise in the reconstruction is audible,
particularly in the spectral valleys between formants where the signal itself
is weak. The postfilter does not remove that noise; it *hides* it, by
deepening the valleys where the ear will not miss the signal and where the
noise therefore stops being audible.

The chain, per subframe:

1. **Two weighted filters.** The subframe's LSFs are pulled towards a fixed
   bias vector by two different amounts, giving a numerator `A(z/g₂)` and a
   denominator `A(z/g₁)`. Deriving them in the LSF domain rather than by
   scaling the LPC coefficients keeps both stable.
2. **Residual.** The synthesised speech goes through `A(z/g₂)`.
3. **Long-term postfilter.** Applied to that residual — see below.
4. **Re-synthesis.** Back through `1/A(z/g₁)`.
5. **Gain normalisation.** So the short-term postfilter cannot change the
   level.
6. **Tilt compensation.** A first-order filter undoing the low-pass tilt the
   formant emphasis introduces; its gain is chosen so it cannot amplify.
7. **Energy matching.** The block is matched in energy to the unfiltered
   speech, through a smoothed automatic gain control, so the postfilter is
   inaudible as a level change.

### The long-term postfilter

This is the part that reinforces periodicity. On the residual:

1. the residual history is rescaled so the search has full precision;
2. the transmitted pitch lag is refined by ±1 sample, then over seven
   fractional phases in both directions, always maximising
   `correlation² / energy`;
3. the delayed signal is built a second time with a longer 16-tap fractional
   filter and whichever version correlates better is kept;
4. the delayed signal is mixed into the residual with a gain derived from the
   winning correlation — two thirds weight to the original unless the delayed
   signal's energy dominates its correlation with the residual.

The lag it settles on is remembered: it is the decoder's voicing flag, read by
the next half-frame's excitation build and by the concealment path.

## Erasure concealment

A frame flagged as suppressed carries no usable parameters. The decoder then:

- **holds the spectrum**, rebuilding the LSF predictor memory from the last
  good frame rather than from the received residual;
- **holds the lag**, using a concealment lag carried from the last good frame;
- **decays the gains** towards zero, deepening the attenuation the longer the
  erasure runs;
- **chooses one contribution, not both** — if the long-term postfilter last
  found a pitch, the adaptive contribution survives and the fixed one is
  dropped; if it did not, the reverse. Summing both would reinforce whichever
  was wrong.
- fills what remains from a pseudo-random generator, so a long erasure decays
  into shaped noise rather than into silence or a frozen buzz.

The decision to keep one contribution rather than both is the important one.
During voiced speech the periodic structure is what the ear tracks, and
extending it is convincing; during unvoiced speech there is no structure to
extend and repeating the adaptive contribution would produce a tone that was
never there.

## Output

The postfiltered signal goes through a final reconstruction filter — two
cascaded second-order sections — and is companded back to μ-law. Each 18-byte
frame becomes 320 μ-law bytes.

An in-band reset frame produces no audio at all: `Decoder::decode` returns
`None` and the decoder returns to its initial state.
