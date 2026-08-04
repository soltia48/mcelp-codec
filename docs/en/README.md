# M-CELP — technical documentation

The Mitsubishi CELP speech codec: a 3.6 kbit/s analysis-by-synthesis coder for
8 kHz telephony speech. This directory describes how it works, stage by stage,
at the level of detail needed to read or modify the implementation in `src/`.

## At a glance

| | |
|---|---|
| Sample rate | 8 kHz, mono |
| Audio format | μ-law, 8 bits per sample |
| Frame | 320 samples = 40 ms |
| Bit stream | 18 bytes per frame = **3.6 kbit/s** |
| Structure | 1 frame = 2 half-frames = 4 subframes of 80 samples |
| Prediction | 10th-order short-term, one long-term (pitch) predictor |
| Excitation | adaptive codebook + fixed codebook, both gain-scaled |

CELP means *code-excited linear prediction*. A short-term linear predictor
models the vocal tract; what it cannot predict — the excitation — is chosen
from two codebooks by synthesising every candidate and keeping the one whose
output sounds closest to the input. "Closest" is measured through a perceptual
weighting filter, so the coder spends its bits where the ear will notice.

Two things distinguish this codec from the ITU-T G.729 family it otherwise
resembles. Its frame is four times longer, which buys the low bit rate at the
cost of delay. And its encoder carries a full spectral noise-suppression front
end: the signal that reaches the analysis has already been cleaned.

## Reading order

Start with [architecture](architecture.md) for the shape of both pipelines,
then follow the signal:

| document | subject |
|---|---|
| [architecture.md](architecture.md) | frame structure, the two pipelines, what state survives a frame |
| [bitstream.md](bitstream.md) | the 18-byte frame and its fourteen fields |
| [noise-suppression.md](noise-suppression.md) | the encoder front end: spectrum, critical bands, voice activity, suppression |
| [lpc.md](lpc.md) | short-term prediction, line spectral frequencies, their quantiser |
| [pitch.md](pitch.md) | long-term prediction: the adaptive codebook and its two searches |
| [fixed-codebook.md](fixed-codebook.md) | the six innovation structures and the 16-bit index that selects one |
| [gain.md](gain.md) | the conjugate gain codebook and the energy predictor behind it |
| [synthesis.md](synthesis.md) | the decoder: synthesis, postfiltering, concealment |
| [fixed-point.md](fixed-point.md) | the arithmetic model everything is computed in |

## Bit budget

One frame carries 138 bits of parameters plus a suppression flag:

| parameter | per frame | bits |
|---|---|---|
| Line spectrum, predictor mode + first stage | 1 | 8 |
| Line spectrum, second stage (two halves) | 1 | 12 |
| Pitch lag, absolute (first subframe of each half) | 2 | 8 each |
| Pitch lag, relative (second subframe of each half) | 2 | 5 each |
| Fixed-codebook index | 4 | 16 each |
| Gain index | 4 | 7 each |
| Frame suppression flag | 1 | 1 |
| | | **139** |

The remaining 5 bits of the 144 transported are framing. The fixed codebook
takes 64 of the 138 parameter bits — nearly half the frame — which is where a
CELP coder's quality comes from at this rate.

## Conventions in these documents

Fixed-point quantities are written in Q notation: `Q15` means a signed 16-bit
value scaled by 2¹⁵, so 1.0 is not representable and 32767 stands for
0.99997. See [fixed-point.md](fixed-point.md).

*Frame* always means 320 samples, *half-frame* 160, *subframe* 80. Several
stages work per half-frame rather than per frame or per subframe, and the
distinction matters.
