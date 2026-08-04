# The fixed codebook

What the adaptive codebook cannot explain is coded by the fixed codebook: a
sparse, mostly-pulse waveform chosen from a structured set. It is the largest
item in the frame — 16 bits per subframe, 64 of the 138 parameter bits.

## Six structures under one index

The 16-bit index is not a flat table lookup. Its range is divided into six
classes, and which class an index falls into decides how the rest of it is
read:

| class | index range | pulses | structure |
|---|---|---|---|
| 0 | 0x0000 – 0x3FFF | 3 | positions on 3 tracks, radices 8 × 16 × 16 |
| 1 | 0x4000 – 0x7FFF | 3 | the same, different tracks and shapes |
| 2 | 0x8000 – 0x98FF | 2 | 40 × 40 positions |
| 3 | 0x9900 – 0xB1FF | 2 | 40 × 40 positions |
| 4 | 0xB200 – 0xBF97 | 2 | 30 × 29 positions |
| 5 | 0xBF98 – 0xFF97 | 10 | twin five-pulse codebook |

A class is selected by finding the largest lower bound the index clears.
Indices at or above 0xFF98 are invalid and clamped.

Arithmetic check, class 0: 8 × 16 × 16 = 2048 position combinations × 2³ sign
combinations = 16384 = the size of its range. Every class is exactly full.

### Classes 0–4: multi-pulse

Subtracting the class base leaves a number whose low bits are one sign per
pulse and whose remainder is a **mixed-radix** number — the pulse position
indices, most significant digit first. Each position index selects an entry
from that track's position table (40 candidates per pulse slot), and each
pulse contributes not a single sample but a **20-sample shape** taken from a
table of twelve waveforms. Three shapes serve the three-pulse classes and two
each the rest.

Using shapes rather than unit impulses is what lets 2–3 pulses cover an
80-sample subframe convincingly: each pulse paints a short waveform rather
than a click.

### Class 5: the shape pair

The top quarter of the index space addresses something different — two
128-entry codebooks, each entry five signed unit pulse positions, with the
innovation being one entry drawn from each. 128 × 128 = 16384, matching the
range. This structure suits subframes that need dense, noise-like excitation
rather than a few strong pulses, and it is searched on its own terms and
offered to the subframe as one more candidate.

## The search

Trying all 65432 indices per subframe is out of the question. The search is
staged.

**1 — Mode.** Before searching, the encoder decides how the subframe should be
coded. If the adaptive codebook alone already covers about two thirds of the
target energy, the subframe is periodic and gets mode 0. Otherwise, if the
target has grown sharply against the previous subframe — an onset — it gets
mode 1, and so do the three subframes after it, so an attack is not cut short
halfway. Everything else is mode 2. The mode chooses which weights the
candidates are compared under.

**2 — Shortlist.** For each track, the target correlation is computed at every
candidate position and only the best few are kept: 3 or 6 per track for the
three-pulse classes, 10 for the two-pulse ones. The shortlist is built by
insertion, so it stays ordered.

**3 — Precompute.** Each pulse response is filtered through the weighted
synthesis filter once, and the cross-energies between every pair of
shortlisted positions are tabulated. After that, scoring a combination is
arithmetic on precomputed numbers rather than filtering.

**4 — Combine.** Every surviving combination of shortlisted positions is
scored by the same `correlation² / energy` ratio the pitch search uses,
compared by cross-multiplication.

**5 — Weigh.** Each class produces a winner and a score; the scores are
brought onto a common exponent and weighted by mode, and the best across all
classes wins. A penalty applies to a class whose pulses land far from where
the target's energy actually sits — measured as the point that splits the
subframe into two halves of equal rectified area.

**6 — Pack.** The winning positions and signs are folded back into a 16-bit
index by the inverse of the decomposition above.

## Decoding

Decoding is the cheap direction and needs no search: classify the index, split
out positions and signs, and lay the shapes into an 80-sample vector. A zero
sign bit subtracts the shape, a set one adds it.

The decoder then applies **pitch sharpening** — the innovation is fed through
the same fractional-delay predictor the adaptive codebook uses, with a gain
carried from the previous subframe. This reinforces periodicity that the
pulses alone would leave ragged, and it is one of the places where the
decoder's output depends on state rather than on the frame in hand.
