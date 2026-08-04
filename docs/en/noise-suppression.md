# The encoder front end

Before any prediction analysis happens, the encoder runs a full spectral noise
suppressor. This is unusual — most CELP coders of this vintage take their
input as given — and it matters: **everything downstream analyses the
suppressed signal, not the input.** The LPC analysis, the pitch search and the
codebook searches all see a cleaned version of the speech.

The front end works on **half-frames**: 160 new samples at a time, two passes
per frame.

## Chain

```
linear input
   │
   ├─ window ────── 256 samples: 96 of history + 160 new, tapered, pre-emphasised
   ├─ FFT ───────── 128-point, giving 129 magnitude bins
   ├─ bands ─────── fold into 20 critical bands, average, convert to dB
   │
   ├─ noise floors ─ two long-term estimates, one adapting about twice as fast
   ├─ features ───── how far this frame sits above each floor, per band and overall
   ├─ vad ────────── an integer speech/silence score from those features
   ├─ weights ────── a per-band suppression depth in dB
   │
   ├─ shaping ────── apply the per-band gains to the spectrum
   ├─ inverse FFT ── back to 256 time samples, then a leaky integrator
   └─ stitch ─────── cross-fade consecutive spans into a continuous signal
```

## Critical bands

The 129 bins are folded into 20 bands of increasing width, roughly following
the ear's frequency resolution: narrow at the bottom where pitch and the first
formant live, wide at the top where the ear integrates anyway. Each band's
bins are summed, divided by the band width, and converted to decibels.

Twenty numbers per half-frame is a small enough description to track over time
cheaply, and a detailed enough one to suppress selectively.

## Two noise floors

The suppressor keeps **two** long-term band-energy estimates rather than one,
differing in adaptation rate — one roughly twice as fast as the other. Each is
updated only on frames its own voice-activity decision called silent, so over
time each converges on the background rather than on the speech.

Two floors, because a single one cannot serve both purposes. A fast estimate
tracks a changing noise environment but is dragged upward by speech that slips
past the detector. A slow one is robust but lags a genuine change in the
noise. Running both and scoring the frame against each gives the detector two
opinions to reconcile.

## Voice activity

The detector turns each floor's description into a small integer score by
testing three features against thresholds:

- how far the current spectrum sits above the floor overall;
- the same per band, summarised;
- the variance of that per-band excess — noise is spectrally flat relative to
  its own average, speech is not.

The score then goes through a counter and a hangover. The hangover is the
important part: a single quiet frame in the middle of an utterance must not
flip the coder into silence mode, because that would let the noise floors
adapt to speech and would suppress the next frame's onset. The thresholds and
hangover lengths are themselves chosen per frame from how loud the frame is.

The second score additionally takes the open-loop voicing measure into
account, and it is the second one the rest of the encoder consumes.

## Suppression depth

How much of each band to take away comes from how speech-like the frame
looked. Three depth profiles exist — deepest for confident silence, shallowest
for confident speech — laid out across the bands with a fixed tilt, then
smoothed against the previous frame's depths so the gain does not jump between
half-frames.

The per-band gains are then applied to the spectrum and the result inverse
transformed. Because a 256-sample span is transformed but only 160 samples are
consumed per half-frame, consecutive spans have to be joined: rather than a
full overlap-add, a 40-sample seam is cross-faded and the rest of each span
taken straight. That is enough because neighbouring spans differ only where
their windows do.

## The trade

Suppression before analysis is a deliberate choice with a cost. Removing noise
from the signal the LPC analysis sees gives a cleaner spectral envelope and
stops the codebook searches from spending bits coding hiss. But any
suppression artefact is baked in — the decoder has no way to undo it — and
over-suppression in a band the detector misjudged shows up as a hole in the
reconstruction. The conservative depth profiles and the heavy smoothing are
what keep that in check.
