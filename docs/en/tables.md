# Reference tables

The codec's constant data: 62 tables, 6235 words in total, in `src/tables/`.
This document gives the numbers a reader most often wants to look up.

## Inventory

| group | tables | words | contents |
|---|---:|---:|---|
| `lsf` | 13 | 2262 | LSF codebooks, MA predictor, cos/acos conversions |
| `fcb` | 13 | 2215 | fixed-codebook structure, positions, shapes |
| `lpc` | 3 | 481 | analysis window, lag window, root-search grid |
| `bands` | 13 | 427 | critical-band geometry, suppression curves |
| `pitch` | 5 | 337 | interpolation filters, lag-to-slot maps |
| `ulaw` | 1 | 256 | μ-law expansion |
| `gain` | 5 | 118 | gain codebook, log/exp tables |
| `analysis` | 2 | 108 | taper, sine table |
| `filters` | 7 | 31 | fixed filter coefficients |

The two codebooks dominate: LSF and fixed-codebook data together are 72% of
the constant space, which is unsurprising — they *are* the coder's model of
speech.

## Critical bands

The 129 magnitude bins are folded into 20 bands. At a 256-sample window and
8 kHz, one bin is 31.25 Hz.

| band | bins | approx. Hz | width |
|---:|---|---|---:|
| 0 | 0 | 0 – 31 | 1 |
| 1 | 1 – 3 | 31 – 125 | 3 |
| 2 | 4 – 6 | 125 – 219 | 3 |
| 3 | 7 – 9 | 219 – 312 | 3 |
| 4 | 10 – 12 | 312 – 406 | 3 |
| 5 | 13 – 16 | 406 – 531 | 4 |
| 6 | 17 – 20 | 531 – 656 | 4 |
| 7 | 21 – 24 | 656 – 781 | 4 |
| 8 | 25 – 29 | 781 – 938 | 5 |
| 9 | 30 – 34 | 938 – 1094 | 5 |
| 10 | 35 – 40 | 1094 – 1281 | 6 |
| 11 | 41 – 47 | 1281 – 1500 | 7 |
| 12 | 48 – 55 | 1500 – 1750 | 8 |
| 13 | 56 – 64 | 1750 – 2031 | 9 |
| 14 | 65 – 74 | 2031 – 2344 | 10 |
| 15 | 75 – 86 | 2344 – 2719 | 12 |
| 16 | 87 – 100 | 2719 – 3156 | 14 |
| 17 | 101 – 118 | 3156 – 3719 | 18 |
| 18 | 119 – 127 | 3719 – 4000 | 9 |
| 19 | 128 | 4000 – | 1 |

The widths follow the ear's frequency resolution — roughly constant below
500 Hz, then growing — which is what makes 20 numbers a reasonable summary of
a speech spectrum. Bands 0 and 19 are single bins at the extremes and carry
little; band 0 is excluded from the spectral-flatness statistics entirely.

## Fixed-codebook classes

| class | index range | size | pulses | radices | shortlist | shapes | positions |
|---:|---|---:|---:|---|---|---|---:|
| 0 | 0x0000 – 0x3FFF | 16384 | 3 | 8 × 16 × 16 | 3, 6, 6 | 0, 1, 2 | 40 |
| 1 | 0x4000 – 0x7FFF | 16384 | 3 | 8 × 16 × 16 | 3, 6, 6 | 3, 4, 5 | 40 |
| 2 | 0x8000 – 0x98FF | 6400 | 2 | 40 × 40 | 10, 10 | 6, 7 | 80 |
| 3 | 0x9900 – 0xB1FF | 6400 | 2 | 40 × 40 | 10, 10 | 8, 9 | 80 |
| 4 | 0xB200 – 0xBF97 | 3480 | 2 | 30 × 29 | 10, 10 | 10, 11 | 59 |
| 5 | 0xBF98 – 0xFF97 | 16384 | 5+5 | 128 × 128 | — | — | — |

Each class is exactly full: multiply the radices by 2^pulses for the sign bits
and you get the size. For example class 4: 30 × 29 × 2² = 3480.

- **size** — indices the class occupies
- **radices** — the mixed-radix digits the class index decomposes into
- **shortlist** — candidate positions kept per track during the search
- **shapes** — which of the twelve 20-sample pulse shapes the class uses
- **positions** — distinct positions the class's tracks cover between them

Valid indices run to 0xFF97; 0xFF98 and above are clamped. Total: 65432 of the
65536 an unrestricted 16-bit field would allow.

## Quantiser dimensions

| quantiser | entries | dimension | bits |
|---|---:|---|---:|
| LSF stage 1 | 128 | 10 | 7 |
| LSF stage 2, low half | 64 | 5 | 6 |
| LSF stage 2, high half | 64 | 5 | 6 |
| LSF predictor mode | 2 | — | 1 |
| Gain stage 1 | 8 | 2 | 3 |
| Gain stage 2 | 16 | 2 | 4 |

The LSF second stage is applied twice — once to each half of the vector — from
the same 64-entry table, so 64 entries do the work of 4096.

## Filters

| filter | taps | where |
|---|---:|---|
| Input high-pass | 2 recursive + 3 direct | `preprocess` |
| Shaping high-pass | 2 recursive + 3 direct | `weighting` |
| Smoothing (symmetric) | 6 | `weighting` |
| Tilt | 5 | `weighting` |
| Pitch interpolation | 31, read at 3 phases | `pitch` |
| Closed-loop lag interpolation | 13, two views | `pitch_search` |
| Postfilter search | 28, 7 phases × 4 | `ltp` |
| Postfilter refinement | 112, 7 phases × 16 | `ltp` |
| Output reconstruction | 2 second-order sections | `postfilter` |
| Analysis window | 400 | `lpc` |
| Suppression window | 256 | `analysis` |

## Reading the values

Table values are `i16` and almost always Q15. Two conventions catch readers
out:

- **Halved coefficients.** The fractional multiply doubles every product, so a
  coefficient that should be 0.8 is stored as 0.4. Two band-width reciprocal
  tables exist for exactly this reason — one whole, one halved — used by the
  paths that do and do not get the doubling back.
- **Shared storage.** A few tables overlap in the original layout. The log2
  fraction table runs into the tables that follow it, and the sine table is
  read both directly and at a quarter-period offset by the rotator. Where the
  code relies on that, the tables here are merged and the offsets named, rather
  than duplicated.
