# Mitsubishi CELP — Technical documentation

This directory contains a technical description of the Mitsubishi CELP
speech codec (三菱CELP方式音声コーデック) as implemented by this crate,
covering the bitstream, LPC analysis, excitation, gain quantization,
synthesis, and the persistent decoder state.

## English

- [Overview / reading order](en/README.md)
- [Decoder architecture and bitstream](en/architecture.md)
- [LPC: LSF → LSP → LPC](en/lpc.md)
- [Excitation: pitch and fixed codebooks](en/excitation.md)
- [Gain quantization](en/gain.md)
- [Synthesis, postfilter, and μ-law output](en/synthesis.md)
- [Decoder state](en/state.md)

## 日本語

- [概要 / 読む順番](ja/README.md)
- [デコーダのアーキテクチャとビットストリーム](ja/architecture.md)
- [LPC: LSF → LSP → LPC](ja/lpc.md)
- [励振: ピッチコードブックと固定コードブック](ja/excitation.md)
- [ゲイン量子化](ja/gain.md)
- [合成、ポストフィルタ、μ-law 出力](ja/synthesis.md)
- [デコーダ状態](ja/state.md)
