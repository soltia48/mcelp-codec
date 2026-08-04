# Long-term prediction — the adaptive codebook

Voiced speech is quasi-periodic: one pitch period looks much like the last.
The long-term predictor exploits that by building part of each subframe's
excitation from the excitation already produced one pitch period ago.

That past excitation, read at a delay, *is* the adaptive codebook. It is
"adaptive" because it is not a fixed table — it is the coder's own recent
output, so it re-forms itself every subframe.

## Lags

| | |
|---|---|
| Range | 20 to 143 samples (56 Hz to 400 Hz at 8 kHz) |
| Resolution | ⅓ sample up to about lag 85, integer above |
| Coding | 8 bits absolute, then 5 bits relative |

Fractional resolution matters more at short lags, where one sample is a large
fraction of a period; above lag 85 the integer grid is already fine enough,
so those codes are spent on reaching further back instead.

The first subframe of each half-frame sends an **absolute** lag in 8 bits. The
second sends a **delta** in 5 bits, covering a window of ten integer lags at ⅓
resolution centred near the first subframe's choice. The window is clamped so
it never runs past either end of the range.

A fractional lag means the past excitation has to be resampled. That is done
with a 31-tap interpolation filter read at one of three phases — the same
filter, indexed forwards for the taps ahead of the sample and backwards for
those behind.

## The two searches

Searching all 124 integer lags with full analysis-by-synthesis for every
subframe would dominate the encoder's cost, so the search is split.

### Open loop, once per half-frame

The open-loop search runs on the weighted signal directly — no synthesis. It
normalises the signal, then makes three passes over decreasing lag ranges:

```
(143 … 63)   (79 … 39)   (39 … 19)
```

Each pass keeps the lag whose correlation is largest once normalised by the
energy of the segment it was taken over. A bias favours the longest candidate,
so a shorter lag has to be clearly better to displace it. This is the standard
guard against **pitch halving** — twice the true period correlates almost as
well as the period itself, and picking the halved lag makes the excitation
buzz.

### Closed loop, once per subframe

The open-loop estimate only brackets the search. The closed loop then tries a
short run of lags around it — seven, clamped to the valid range — and for each
one builds the actual filtered contribution and measures how much of the
target it explains. Once the best integer lag is found, the fractional phases
around it are tried the same way.

The criterion throughout is `correlation² / energy`, compared between
candidates by cross-multiplication so no division is needed.

## Gain, and why it is capped

The adaptive contribution is scaled by a pitch gain, quantised jointly with the
fixed-codebook gain (see [gain.md](gain.md)). Before quantisation the encoder
clamps it: a periodic subframe is not allowed a gain above roughly 0.95.

The reason is stability of a different kind. The long-term predictor is a
feedback loop through the coder's own output. A gain at or above one makes the
loop ring, and because the decoder runs the same loop, the ringing survives
into the reconstruction and grows. Capping the gain costs a little prediction
accuracy on strongly voiced frames and buys unconditional stability.

## Periodic extension

When the lag is shorter than a subframe, the codebook vector runs off the end
of the history — there is no excitation one period ahead of where we are
writing. The missing tail is filled by repeating what has just been produced,
resampled at the fractional part of the lag by an eight-tap interpolator. The
extension feeds on itself: once more than one period has been produced, the
interpolator's window reaches into samples that same loop wrote.
