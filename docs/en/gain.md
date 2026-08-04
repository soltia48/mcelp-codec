# Gains

Each subframe's excitation is the sum of two scaled contributions:

```
excitation = gₚ · (adaptive codebook) + g_c · (fixed codebook)
```

Both gains come from a **single 7-bit index**. Sending them jointly rather
than separately exploits the fact that they are correlated — a strongly voiced
subframe wants a large pitch gain and a small code gain, and vice versa — so a
joint codebook covers the useful region far more efficiently than two
independent ones.

## The conjugate codebook

The 7 bits split into two stages:

```
index = [ 3 bits ][ 4 bits ]
           │         └── second stage: 16 rows × (pitch term, code term)
           └──────────── first stage:   8 rows × (pitch term, code term)
```

The two rows are summed component-wise. 8 × 16 = 128 combinations from a table
of only 24 rows — that is the "conjugate" structure, and it makes the search
cheap enough to run exhaustively.

The first component of the sum is the pitch gain, used directly. The second is
**not** the code gain. It is a correction factor applied to a *predicted* code
gain.

## Why the code gain is predicted

The fixed-codebook gain has to span the whole dynamic range of speech —
silence to a shout — which would cost many bits to code directly. But it is
strongly predictable: the energy the innovation needs this subframe is close
to what it needed last subframe, once the innovation's own energy is
accounted for.

So the coder predicts it:

1. take the innovation's energy for this subframe;
2. take a 4-tap moving average of the *prediction errors* of the last four
   subframes, kept in dB;
3. combine them to predict the gain this subframe should need;
4. send only the correction factor — how far the real answer sits from the
   prediction.

The correction spans a much smaller range than the gain itself, so 4 bits of
codebook carry it. The predictor memory is one of the pieces of decoder state
that makes a stream undecodable from the middle: get the four past errors
wrong and every subsequent code gain is wrong too.

Everything here is done in the logarithmic domain, which is why the module
carries `log2` and `2^x` tables. A gain is a multiplicative quantity; averaging
gains in dB is averaging the thing that is actually stationary.

## The encoder's search

The closed loop leaves behind five measures describing the subframe: two
energies and three correlations, relating the target, the adaptive
contribution and the fixed one. Every candidate's coding error can be written
as a weighted sum of five terms in the two gains —

```
p², p, v², v, p·v
```

— where `p` is the code gain and `v` the pitch gain. So the search does not
have to synthesise anything: it evaluates that quadratic for all 128 entries
and keeps the smallest.

The five measures arrive at different scales, since some are energies and
some correlations. They are first brought onto a common exponent: a measure
more than 16 bits below the largest has its mantissa shifted down and the rest
of the gap carried in a shift, and one more than 32 bits below is dropped
entirely.

## Erased frames

When a frame is lost, the decoder cannot read a gain index. Instead both gains
decay towards zero over successive lost frames, and the energy predictor's
memory is pulled down towards its reset floor. The attenuation deepens the
longer the erasure runs, so a burst fades out rather than either freezing on
the last good gain or cutting off abruptly.
