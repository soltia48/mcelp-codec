# LPC: LSF → LSP → LPC

三菱 CELP は、スペクトル包絡を 10 次の **線スペクトル周波数 (LSF)**
ベクトルとして送出し、これを split-VQ + フレーム間予測でエンコードしています。
デコーダ側では、フレームの LSF が LSP に変換され、ソースが
`response_shaper` と呼ぶ 3 パスの自己相関パイプラインを通して
サブフレームごとの LPC 係数に展開されます。本ドキュメントはこの一連の流れを
端から端まで説明します。

ソースファイル:

- [src/lsf.rs](../../src/lsf.rs) — LSF 復号と LSP 変換
- [src/lpc_analysis.rs](../../src/lpc_analysis.rs) — 自己相関 + Levinson +
  安定性チェック
- [src/tables/lsf.rs](../../src/tables/lsf.rs) — VQ コードブックと予測子
- [src/tables/lsp.rs](../../src/tables/lsp.rs) — cosine LUT

## 1. LSF 復号

フレーム最初の 2 つのビットストリームフィールド (F[0], F[1]) が LSF
インデックスを運びます。

```rust
struct LsfIndices {
    mode:  bool,  // F[0] bit 7        — 予測子 / 平滑化行を選択
    seed:  u8,    // F[0] bits 0..6    — stage-1 コードブックインデックス (7 ビット)
    upper: u8,    // F[1] bits 6..11   — k=0..4 用 stage-2 コードブックインデックス
    lower: u8,    // F[1] bits 0..5    — k=5..9 用 stage-2 コードブックインデックス
}
```

つまり 10 次元 LSF に対して **8 + 12 = 20 ビット** で、そのうち 1 ビットが
予測子モード、19 ビットが split-VQ のインデックスです。`upper` と `lower` は
LSP インデックスの範囲ではなく、**F[1] 内のビット位置** を指す名前です。

### 1.1 split-VQ コードブックの組み合わせ

`combine_codebook` ([`combine_codebook`][lsf] 参照) は 7 ビットの第 1 段
ベクトルと 2 つの 6 ビットの第 2 段ベクトルを足し合わせます。

```
scratch[0..5]  = STAGE1_CODEBOOK[seed][0..5]  + STAGE2_CODEBOOK[upper][0..5]
scratch[5..10] = STAGE1_CODEBOOK[seed][5..10] + STAGE2_CODEBOOK[lower][5..10]
```

テーブルは [src/tables/lsf.rs](../../src/tables/lsf.rs) に Q15 で格納されています
(`STAGE1_CODEBOOK` は 128×10、`STAGE2_CODEBOOK` は 64×10)。

### 1.2 最小ギャップの強制

組み合わせ後、`enforce_min_gap` が 2 回 — 1 回目は `min_gap = 10`、2 回目は
`min_gap = 5` で — 呼び出され、各ペア `(scratch[i], scratch[i+1])` を
`scratch[i+1] − scratch[i] ≥ −min_gap` となるように引き離します。
2 回 (10 → 5) に分けるのは、リファレンス DSP が用いる 2 段階の平滑化スケジュールを
再現したものです。

### 1.3 予測的組み合わせ

フレーム間予測は `predictive_combine` によって、過去 3 フレーム分の
**スクラッチ** ベクトルを保持する 3 スロットの LSF 履歴 (`LsfHistory`、
最古がスロット 2) に対して計算されます。各 `i ∈ 0..10` について:

```
out[i] = ( scratch[i]    · SMOOTH[mode][i]
         + history[0][i] · PREDICTOR[mode][i]
         + history[1][i] · PREDICTOR[mode][i + 10]
         + history[2][i] · PREDICTOR[mode][i + 20]
         ) · 2 / 65536
```

(`× 2 / 65536` は DSP の `frct_mul=2` の後に `>>16` を行う動作の再現。)

2 つの予測子モード (`mode ∈ {0, 1}`) は別の予測子 / 平滑化行を提供します。
選択ビットは F[0] の MSB です。一般に mode 0 は予測の度合いが弱い設定、
mode 1 は予測の度合いが強い設定に対応します。

予測的組み合わせの後、(予測出力ではなく) *生のスクラッチ* が
`LsfHistory::splice_new_vector` によって履歴に組み込まれます。

### 1.4 安定化

予測出力は `stabilize_and_finalize` によって最終化されます。

1. **Bubble-1 パス** が隣接する順序乱れの要素を入れ替え。
2. **下限クランプ** `lsf[0] ≥ 41` (Q15 で約 0.00125)。
3. **前向き最小ギャップ = 321** で LSP の分離可能性を保証。
4. **上限クランプ** `lsf[9] ≤ 25682` (Q15 で約 0.7838)。

このパスの後、LSF は単調増加かつ有界となります。後段の LSP 変換は
これらの不変条件に依存します。

## 2. LSF → LSP

`lsf_to_lsp` ([src/lsf.rs](../../src/lsf.rs) 参照) は、正規化された LSF
(Q15、範囲 `[0, 0.5)`) を対応する LSP (= `cos(π · lsf)`) に写す離散コサイン変換を
実装します。

```
scaled       = lsf · 20861 · 2     // i64
index        = clamp(scaled >> 24, 0, 63)
fraction_bits= (scaled >> 16) & 0xFF
lsp = COS_LUT_VALUE[index] + (COS_LUT_SLOPE[index] · fraction_bits) >> 12
```

定数 20861 は `64 · 2^15 / π ≈ 20861` で、`scaled >> 24` が 64 エントリの
コサイン LUT への整数インデックスとなり、残り 8 ビットの小数 (Q8) が
事前計算された傾きを使ってテーブルセル間を線形補間します。
`COS_LUT_VALUE` と `COS_LUT_SLOPE` は [src/tables/lsp.rs](../../src/tables/lsp.rs)
にあります。

デコーダはまた `prev_frame_lsp` を保持しており、`block_average_lsp`
(成分ごとの `(prev + curr) >> 1`) によって **block-0 の作業係数** —
すなわち前フレームの LSP と現フレームの LSP の平均 — を生成できるようにしています。
block-1 の作業係数は単純に現フレームの LSP そのものです。

## 3. 3 パスの `response_shaper` パイプライン

このコーデックの LPC 解析で興味深いのは、合成フィルタ係数を現フレームの LSP から
**直接** 計算しないという点です。代わりに `lib.rs` は、3 つの異なる「ウィンドウ」を
用いて `build_autocorrelation_with_state` をフレームあたり 3 回呼び出し、
得られた係数をサブフレームごとの LPC に分割します。各パスは共有された 35 セルの
スクラッチバッファ (`DecoderState` の `response_shaper_buffer`) を変更するため、
シーケンスは順序依存です。

### 3.1 パス A — 安定性チェック

入力ウィンドウ: `lsp_q15` の後ろに 4 つの 0 を追加 (`stability_seed_window =
[lsp_q15..., 0, 0, 0, 0]`)。

```
(stability_autocorr, post_buf) = build_autocorrelation_with_state(
    prev_frame_lsp,            // 「作業」係数
    stability_seed_window,     // 「シード」ウィンドウ
    response_shaper_buffer,
);
response_shaper_buffer ← post_buf
```

デコーダは `stability_autocorr` の **後半** (セル `12..22`) を仮の LPC
ベクトルとして取り出し、それに対して 10 ステップの Levinson-Durbin 再帰
(`levinson_recursion`) を実行して 10 個の反射係数 $k_1, \dots, k_{10}$ を得ます。
これらは `is_lpc_stable` に渡され、`true` が返るのは $|k_1| < 1$ かつ末端の指標 —
`d14b_metric_residual(k_i)` を `i ∈ 2..10` にわたって 32-bit dword 対として
累算したもの — も `STABILITY_REF` のしきい値をクリアした場合に限られます。

**block-0 LPC オーバライド** の判定規則は以下です。

```
override_block_0 = stable_now && prev_lpc_stability_flag == 0
```

すなわちオーバライドは安定性の *立ち上がりエッジ* でのみ発火します。
オーバライドが発火すると、block-0 の作業係数はブロック平均ではなく
生の `lsp_q15` で置換されます。`prev_lpc_stability_flag` はその後更新されます。
この立ち上がりエッジゲートは、安定性の境界付近でのジッタを防ぐためのヒステリシスです。

### 3.2 パス B — block-0 自己相関

入力ウィンドウ: `[block_0_coeffs_eff..., lsp_q15[0..4]]`。

```
(autocorr_block0, post_buf) = build_autocorrelation_with_state(
    prev_frame_lsp,
    block0_seed_window,
    response_shaper_buffer,
);
response_shaper_buffer ← post_buf
sub0_lpc = autocorr_block0[1..11]
sub1_lpc = autocorr_block0[12..22]
```

`build_autocorrelation_with_state` の内部では:

1. 14 セルの **first_autocorr** = `[(work + seed[0..10])>>1, seed[0..4]]` を構成。
2. `response_shaper(first_autocorr, …)` を呼んで `mirror_output[0..11]` を得る。
   これがこのパスの最初のサブフレーム用 11-tap LPC。
3. `response_shaper(seed_window, …)` を呼んで 2 番目のサブフレーム用の
   `mirror_output[0..11]` を得る。
4. これらを連結して 22 セルの結果にする。

`response_shaper` 単体も解析ステージです: 35 セルの dword バッファの前半・後半を
埋める 2 回の `levinson_step` 累算、対称ペアを加算し反対称ペアを減算する
**マージフェーズ** (`[10, 8, 6, 4, 2]` の和と `[22, 20, 18, 16, 14]` の差を 5 回反復)、
そして両半分の和と差を読み出して 11 セルの LPC を出力する **ミラーフェーズ**
(`mirror[0] = 4096 = Q12 1.0` は常に固定) から成ります。数学的な直観としては、
前半が *対称* (P 多項式風) の半分を、後半が *反対称* (Q 多項式風) の半分を
エンコードしており、ミラーフェーズは LSP→LPC の高速再構成です。
11 セルの出力は Q12 です。

加えて、`levinson_recursion` が各サブフレームの 10 セル LPC 末尾
(= `autocorr_block0` のセル `1..11` と `12..22`) に対して呼び出され、
そのサブフレームの `clamp_input` として $k_1$ が抽出されます — ゲインパイプラインは
これを synth-control 判定の `clamp_input` として使用します。

### 3.3 パス C — block-1 自己相関

入力ウィンドウ: `stability_seed_window` (パス A と同じ
`[lsp_q15..., 0, 0, 0, 0]`、ただし作業係数の入力には `block_0_coeffs_eff` を使用)。

```
(autocorr_block1, post_buf) = build_autocorrelation_with_state(
    block_0_coeffs_eff,
    stability_seed_window,
    response_shaper_buffer,
);
response_shaper_buffer ← post_buf
sub2_lpc = autocorr_block1[1..11]
sub3_lpc = autocorr_block1[12..22]
```

概念的には block 1 用の対称的な構成であり、「作業」(= block-0 係数) と
「シード」(= 現フレームの LSP) の役割は、block 1 の最初のサブフレームが
block-0 出力と現フレーム LSP の平均、block 1 の 2 番目のサブフレームが
現フレーム LSP そのものとなるように選ばれています。

### 3.4 なぜ 3 パスなのか?

各パスの役割:

| パス | ウィンドウの組み合わせ                       | 役割                                              |
| ---- | -------------------------------------------- | ------------------------------------------------- |
| A    | `prev_lsp ↔ [lsp                     0]`     | block-0 オーバライドを判定 (安定性チェック)       |
| B    | `prev_lsp ↔ [block_0_coeffs   lsp[0..4]]`    | sub 0 と sub 1 の LPC                             |
| C    | `block_0_coeffs ↔ [lsp                0]`    | sub 2 と sub 3 の LPC                             |

3 つのパスは 35 セルのスクラッチバッファを決定論的な前進順で共有します。
いずれかを単体の解析で置き換えるとビット一致が崩れます — これは dword
スクラッチがパス間で情報を運ぶためで、Rust 実装が忠実に再現している
オリジナル DSP ファームウェアの性質です。

### 3.5 反射係数 $k_i$ を `clamp_input` として用いる

各サブフレーム `i ∈ 0..4` について、値
`clamp_inputs[i] = levinson_recursion(sub_i_lpc)[0]` (すなわち 10 セルの
自己相関末尾に対して 10 ステップの Levinson 再帰を回して得られる第 1 反射係数
$k_1$) は synth-control 判定 ([synthesis.md](synthesis.md#1-synth-control-ピッチ強調のゲーティング) を参照) に渡されます。
$k_1$ は古典的な有声性の指標 (高度に有声なフレームでは `k_1 → -1`) なので、
デコーダはそれを使って当該サブフレームに対するピッチ周期性強調 IIR を
有効化するかどうかを判断します。

別途、Q12 LPC ベクトルに対して逆 Levinson (Schur 再帰) を実装する
`lpc_to_first_reflection` というユーティリティが [src/lsf.rs](../../src/lsf.rs)
に用意されています。これは現状ランタイムの復号経路では使われていませんが、
エンコーダ側の解析との対称性とユニットテスト用に保持されています。

## 4. Levinson-Durbin 再帰

`levinson_recursion` は古典的な Levinson-Durbin 再帰を、10 ステップに特化して、
コーデックの他部分と共有する i40 アキュムレータセマンティクスで実装したものです。
この再帰は、10 セルの自己相関末尾を入力として in-place で反射係数
$k_1, \dots, k_{10}$ を生成します。

各内部ステップは以下を計算します。

```
E_i      = E_{i-1} · (1 - k_i^2)
k_i      = -⟨a_{i-1}, R_i⟩ / E_{i-1}      (前向きステップ)
a_i[m]   = a_{i-1}[m] - k_i · a_{i-1}[i-m]  for m ∈ 1..i  (後ろ向きステップ)
```

特徴は次の 2 点です:

- **`q15_reciprocal`** は反復的な `subc` ベース除算によって `1 / E_{i-1}` を
  Q15 で計算します (DSP の `subc` 命令を `subc_step` として実装したものを 15 ステップ、
  続いて 32-bit マスクのスクラブと 40-bit 符号拡張)。これは DSP ファームウェアの
  逆数推定をビット単位で再現します。
- シフト量 `asm` は逆数の正規化指数の上位 5 ビットから導出され、5 ビットからの
  符号拡張を受けます。

`levinson_recursion` が返す Q15 反射係数は、安定性チェック (`is_lpc_stable`) と、
ゲインパイプラインのサブフレームごとの `clamp_input` 参照の両方に使用されます。

## 5. サブフレームごとの LPC まとめ

3 つのパスが完了すると、デコーダは **4 つの 10 セル Q12 LPC ベクトル** —
`sub0_lpc`, `sub1_lpc`, `sub2_lpc`, `sub3_lpc` — と **4 つの Q15 第 1 反射係数** —
`clamp_inputs[0..4]` — を保持します。この時点で LPC 解析は終了し、
合成パイプラインはこれらの係数をそれ以上変換しません。

[lsf]: ../../src/lsf.rs
