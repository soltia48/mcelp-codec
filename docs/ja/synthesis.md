# 合成、ポストフィルタ、μ-law 出力

励振生成とゲイン量子化の後に残るステップは以下の 4 つです:
(1) ピッチ強調と合成フィルタの両方で使われる synth-control 判定、
(2) LPC 合成フィルタそのもの、(3) ハーフフレームポストフィルタ、
(4) 線形 PCM から μ-law への変換。本ドキュメントでこれら 4 つすべてを扱います。

ソース: [src/synth.rs](../../src/synth.rs),
[src/postfilter.rs](../../src/postfilter.rs),
[src/ulaw.rs](../../src/ulaw.rs),
[src/tables/postfilter.rs](../../src/tables/postfilter.rs),
[src/tables/ulaw.rs](../../src/tables/ulaw.rs)。

## 1. <a name="1-synth-control-ピッチ強調のゲーティング"></a> Synth control (ピッチ強調のゲーティング)

`compute_pitch_enhance_gain` は、現サブフレームに対してピッチ周期性強調 IIR
([excitation.md](excitation.md#34-ピッチ周期性強調-pitch_enhance) 参照)
を有効化するかを決定します。入力は:

- `fcb_index` — このサブフレームの FCB 符号語
- `current_lag` — このサブフレームのピッチ lag
- `clamp_input` — このサブフレームの第 1 反射係数 $k_1$
  (`levinson_recursion` の出力。
  [lpc.md](lpc.md#35-反射係数-k_i-を-clamp_input-として用いる) を参照)
- `prev_lag` — 前ハーフフレームで選択された lag
- `prev_pitch_gain` — 前サブフレームの $g_p$

lag クラスでのディスパッチはハードコードされています。

| lag クラス | 動作                                                |
| ---------- | --------------------------------------------------- |
| 0, 3, 4    | 常に **有効** (gain = 16384 = Q15 0.5)              |
| 1, 2, 5    | `update_synth_control` の判定に従う                 |

`update_synth_control` 自体は 4 分岐のヒステリシスを実装します。まず
`prev_pitch_gain` を `[MIN_SYNTH_STATE, MAX_SYNTH_STATE] = [3277, 13017]` に
クランプして **前制御状態** として使い、その上で次のように判断します。

```
if clamp_input が大きすぎる (≈ |k_1| > 0.5):
    return DISABLE                    (gain = 0)
elif prev_control_state ≤ SYNTH_CONTROL_GATE = 13107:
    return prev_control_state         (持ち越し)
elif current_lag ≤ prev_lag · 24576/32768:
    return prev_control_state         (lag のドリフトが緩やか → 維持)
elif current_lag <  prev_lag · 24576/32768 · 2:
    return ENABLE                     (gain = 16384)
else:
    return prev_control_state         (持ち越し)
```

言葉にすると: LPC が非有声フレームを示す ($k_1$ が大きすぎる) ときは IIR を
**無効** にし、それ以外でピッチ lag が *素早く変化した* (前 lag を 24576/32768
でスケールした値の 1× から 2× の間でジャンプした) ときに有効化する — つまり、
フォーマットは話者が無声から有声に切り替わるとき、エンコーダがピッチに段階的
変化を駆動することを期待しています。lag クラス 0, 3, 4 — フォーマットが強く
有声な信号のために予約している FCB 範囲 — では、IIR は無条件に有効です。

このルーチンの定数は:

```
SYNTH_CONTROL_ENABLE  = 16384   (Q15 0.5 — 有効時のゲイン)
SYNTH_CONTROL_DISABLE = 0
SYNTH_CONTROL_GATE    = 13107   (≈ 0.4)
PREV_LAG_BLEND_GAIN   = 24576   (= 0.75)
MIN_SYNTH_STATE       = 3277
MAX_SYNTH_STATE       = 13017
```

## 2. 励振の混合

`mix_excitation(v, c, g_p, g_c)` は `n ∈ 0..80` の範囲で Q15 の励振
`e[n] ∈ i16` を生成します。

```
acc = g_c · c[n] · 8 + g_p · v[n] · 4 + 32768       (i64)
acc = sat32(acc)
e[n] = sat16(acc >> 16)
```

`· 8` および `· 4` の係数は、$g_c$ と $g_p$ の内部 Q フォーマットを合成
フィルタが期待する Q15 入力に揃える固定スケール補正です。`+ 32768` は
最終 `>> 16` 用の最近接丸めバイアスです。

## 3. LPC 合成フィルタ ($1/A(z)$)

[src/synth.rs](../../src/synth.rs) の `lpc_synthesis_filter(excitation, lpc_coeffs, history)` は
ストレートな全極 IIR です。

```
buf[0..10]  = history (前サブフレーム出力の最後 10 サンプル)
for n in 0..80:
    mac = Σ_{k=0..9} lpc_coeffs[k] · buf[9 + n - k]   (Q12 × Q15 = Q27)
    acc = (excitation[n] << 16) - mac · 16 + 32768   (i64)
    acc = sat32(acc)
    buf[10 + n] = sat16(acc >> 16)
output[0..80]    = buf[10..90]
new_history[0..10] = buf[80..90]
```

Q フォーマットの選択:

- `excitation[n]` は Q15。`excitation[n] << 16` で Q31 に持ち上げる。
- `lpc_coeffs[k]` は Q12。`mac` は Q27。`mac * 16` で Q31 に持ち上げる。
- アキュムレータは i32 飽和、その後丸め付きで右シフトし Q15 にする。

LPC 係数は `build_autocorrelation_with_state` のミラー出力から直接得られる
ため ([lpc.md](lpc.md) 参照)、合成パスでは **個別の LSP→LPC 変換は実行しません** —
`response_shaper.mirror_phase` が Q12 LPC を直接生成します。

`lpc_coeffs[0]` (暗黙の `a_0 = 1.0`) は合成フィルタに渡される 10 要素ベクトルに
**含まれません** — フィルタは $a_1 \dots a_{10}$ をループで処理し、
$a_0 \cdot e[n]$ の項は別途 `excitation[n] << 16` として加算します。

あるサブフレームから取り出された `history` は、次のサブフレームの `history`
となります。デコーダは `DecoderState` の `lpc_synth_history` としてそれを
保持します。

## 4. 過去励振バッファの管理

320 セルの `past_excitation` バッファは過去 2 ブロック分の励振サンプルを
保持します。

```
[0..160]   ← 前ブロックの e[n]
[160..240] ← 現ブロックの sub a (位相 0)  e[n]
[240..320] ← 現ブロックの sub b (位相 80) e[n]
```

各サブフレームは現ブロックスロットに 80 サンプルの励振を書き込み、
*次の* サブフレーム (または次のブロック) の `pitch_adaptive_codebook` が
`effective_lag` の後方オフセットでこのバッファから読み出します。

ブロックの後 — すなわち 2 番目のサブフレームの後 (`i % 2 == 1`) —
[src/lib.rs](../../src/lib.rs) は `[160..320]` を `[0..160]` にシフトダウンし、
`[160..320]` をクリアして、次ブロックの synth-control 判定が正しい参照を
持てるよう `prev_lag = lag_int` をコミットします。

## 5. ポストフィルタ

4 つのサブフレームすべてが合成された後、デコーダは 320 サンプルの Q15 バッファ
`synth_pcm` を保持します。ポストフィルタはその後、[src/postfilter.rs](../../src/postfilter.rs)
の `postfilter_apply` によって **ハーフフレーム** 単位 (160 サンプル = ブロック)
で適用されます。

ポストフィルタは **前向き + 後ろ向き IIR** で、フォルマント強調にチルト補正と
最終出力ゲイン調整を組み合わせています。サンプルごとに:

```
input_acc         = sample << 16
b1 = input_step(input_acc, fwd_hist, INPUT_COEFF)  (4-tap MAC + 3 セルのシフトレジスタ)
b2 = feedback_correction(b1, delay_fwd, FEEDBACK_COEFF)  (3 つの IIR フィードバックセクション)
shift_delay_line(delay_fwd, b2)
rounded = round_acc(b2)

b3 = input_step(rounded, rev_hist, INPUT_COEFF)
b4 = feedback_correction(b3, delay_rev, FEEDBACK_COEFF)
shift_delay_line(delay_rev, b4)

pf = output_sample(b4, OUTPUT_GAIN)
sample = (pf >> 16) & 0xffff
```

固定係数は:

```
INPUT_COEFF    = [7472, -22381,  22381, -7472]   (4-tap 入力 MAC)
FEEDBACK_COEFF = [-23033, 21666,  -6815]         (3-tap フィードバック IIR)
OUTPUT_GAIN    = 19661                           (出力トリム ≈ Q15 で 0.6)
```

ポストフィルタは 2 種類の状態を保持します。

- `postfilter_history[6]` — 3 セルの前向き + 3 セルの後ろ向き履歴シフト
  レジスタ
- `postfilter_delay[12]` — 6 セルの前向き + 6 セルの後ろ向き dword ペア
  ディレイライン (各 `feedback_correction` は 3 つの dword ペアを読む)

前向き段は各々の前半を、後ろ向き段は後半を使います。両半分とも `DecoderState`
にハーフフレーム間 (および同様にフレーム間) で永続化されます。

`round_acc` と `output_sample` は DSP ファームウェアの飽和順序
(`<< 1`, `+ 0x8000`, `sat32`) を再現する小さなヘルパで、各中間ステップを
オリジナル DSP が用いる 16 ビットの粒度に丸めます。

## 6. μ-law 変換

最終的な 320 個の Q15 サンプルは、[src/ulaw.rs](../../src/ulaw.rs) の
`linear_i16_to_ulaw` によって μ-law バイトに変換されます。これは標準の
G.711 μ-law エンコーダを実装します。

```
m0 = sample if sample >= 0 else !sample        (1 の補数の絶対値)
m  = min((m0 >> 2) + 33, 8191)                 (バイアス + 13 ビット飽和)

seg_bits  = m >> 6                              (セグメント 0..127)
seg_shift = 0 if seg_bits == 0 else floor(log2(seg_bits)) + 1
mantissa  = (m >> (seg_shift + 1)) & 0x0F
segment   = 7 - seg_shift

ulaw = (segment << 4) | (15 - mantissa)
if sample >= 0: ulaw |= 0x80
```

逆方向は 256 エントリの LUT (`ULAW_TO_LIN`) です。これは対称性のために
提供されていますが、**復号パスでは使われません** — このコーデックは μ-law
を出力するだけで、自身の操作の一部として線形 PCM に戻す復号は行いません。

細かな点: エンコーダは負ゼロを正ゼロに畳み込みます
(`linear_i16_to_ulaw(0) == 0xFF`)。これはテーブルラウンドトリップが厳密な
逆ではないことを意味します — バイト 127 と 255 はどちらも線形 0 にマップされ、
そこから 0xFF に戻ります。これは標準的な G.711 の負ゼロ畳み込みであり、
電話機テレフォニースタックの他部分がフォーマットを解釈する仕方と一致します。

## 7. なぜ μ-law 出力なのか?

このコーデックは G.711 / μ-law ベースのテレフォニーエンドポイントへの統合を
意図しており、そこではオーディオパスの入出力がすでに μ-law です。μ-law を
直接出力することで、デコーダは μ-law シンクへ送出するときの追加の線形化
ステップを 1 つ省きます。`Codec::decode_frame` API はこれを反映して、線形
PCM の `[i16; 320]` ではなく μ-law の `[u8; 320]` を返します。
