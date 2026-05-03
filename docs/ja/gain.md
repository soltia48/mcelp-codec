# ゲイン量子化

各サブフレームのピッチゲイン $g_p$ と FCB ゲイン $g_c$ は、ビットストリームの
F[4], F[7], F[10], F[13] が運ぶ **7 ビットの位相ワード** に共同で符号化されます。
復号側でそれらを復元するのは単純ではありません — 位相テーブル参照、いくつかの
固定テーブルを辿る予測リフレッシュ、予測状態の指数的履歴更新、そして
しきい値駆動のテールデケイリミッタを伴います。

本ドキュメントはそのパイプラインを順に解説します。
ソース: [src/gain.rs](../../src/gain.rs), [src/tables/gain.rs](../../src/tables/gain.rs)。
オーケストレータのエントリーポイントは `gain_orchestrate_codec` です。

## 1. 入出力

### 入力

```rust
struct GainOrchestrateCodecInput {
    phase_word:    i16,    // ビットストリームからの 7 ビットフィールド
    threshold_in:  i16,    // テールデケイのしきい値 (フレーム永続)
    counter_in:    i16,    // テールデケイのカウンタ (フレーム永続)
    candidate:     [i16; 80], // このサブフレームの FCB 出力 c[n]
    history:       [i16; 5],  // 過去ゲイン履歴 (フレーム永続)
}
```

### 出力

```rust
struct GainOrchestrateOutput {
    pitch_gain:     i16,   // g_p, mix_excitation の pitch_gain として渡す
    fcb_gain:       i16,   // g_c, mix_excitation の fcb_gain として渡す
    history_out:    [i16; 5],
    threshold_out:  i16,
    counter_out:    i16,
}
```

`phase_word` は両ゲインで共有されます — つまりエンコーダはサブフレームごとに
1 つの 7 ビットインデックスを送出し、それが共同量子化された 128 個の
`(g_p, g_c)` のペアの中から 1 つを選択します。他の入力は、サブフレーム単位で
更新されるフレームごとの状態 ([state.md](state.md) 参照) から得られます。

## 2. 位相セットアップ

`gain_phase_setup(phase_word, …)` は位相ワードを 2 つの 4 ビット半分に分解し、
それぞれを使って `PHASE_TABLE` の 2 セル ("upper" と "lower" のペア) を参照します。

```
phase_first_idx  = (phase_word >> 4) << 1     // 0..14
phase_second_idx = ((phase_word & 0xF) << 1) + 16   // 16..46

phase_first_upper = PHASE_TABLE[phase_first_idx]
phase_first_lower = PHASE_TABLE[phase_first_idx + 1]
phase_second_upper= PHASE_TABLE[phase_second_idx]
phase_second_lower= PHASE_TABLE[phase_second_idx + 1]
```

`PHASE_TABLE` は [src/tables/gain.rs](../../src/tables/gain.rs) にある
48 セルの量子化テーブルです。`phase_word` の **上位 3 ビット**
(すなわち `phase_word >> 4`、範囲 0..7) が `PHASE_TABLE[0..16]` のセル対を
選択し、**下位 4 ビット** (範囲 0..15) が `PHASE_TABLE[16..48]` のセル対を
選択します。概念的には 2 つの半分は $g_c$ と $g_p$ の量子化値をそれぞれ
運びますが、そのまま使われるのではなく、後段の予測リフレッシュを通します。

`gain_phase_setup` は続いて以下を計算します。

- `update_base = ((phase_first_upper + phase_first_lower) << 16) >> 16`
  — $g_p$ の初期候補。
- `saved_ah    = ((phase_second_upper + phase_second_lower) << 16) >> 16`
  — 履歴セットアップに渡す簿記値。

また、予測リフレッシュの入力となる位相ワードのビットでシードされた 2 つの
i40 アキュムレータ (`acc_main`, `acc_scratch`) を導出します。

## 3. 予測リフレッシュ

`gain_predictive_refresh` が支配的な計算です。これは
`history_shift_init` (整数のシフト数) と `history_sample_init` (Q15 サンプル) を
生成し、両者でゲイン量子化器が正規化対象とする **励振リファレンスレベル** を
定義します。関数の 10 フェーズはオリジナル DSP ルーチンを反映しますが、
高レベルの構造は次のとおりです。

```
1. FCB 候補のエネルギー:
       E = 10·2^17 + 2·Σ c[n]^2
   (10·2^17 は小エネルギーの下限)
2. PREDICT_NORMALIZE_TABLE に対する一次正規化
       — 正規化エネルギーから参照インデックスを選び、
         2 セル間を補間 (Q15 frac)。
3. 線形補正 (19488<<19 を引き、18432<<20 を加える)。
4. 過去ゲイン履歴との 4-tap 内積:
       acc_main = Σ_{k=0..3} COEFF_TABLE[k] · history[k]
   ここで COEFF_TABLE = [5571, 4096, 2703, 1311] は 4 極の AR 予測子。
5. 結合 + 5439 倍率での混合 → history_shift_init。
6. PREDICT_SECONDARY_TABLE に対する二次正規化
       (一次と同じ形状、ただし係数が異なる)。
7. history_sample_init = 二次後アキュムレータの上位 16 ビット。
```

2 つの正規化ステップは `gain_normalize_primary` と `gain_normalize_secondary`
で実装されています。両者ともに **区分線形** な (異なる) 非線形関数の近似で、
i40 アキュムレータの正規化された仮数部でインデックスする 33 セルの参照
テーブルと、セル間を補間する Q15 の小数重みを持ちます。構造は次のとおり。

```
acc          ← acc を 5 だけシフト (Q?? → Q??)
t_exp        ← exp_acc(acc)              (Q40 正規化)
acc          ← norm_acc_with_t(acc, t_exp)
lookup_index ← acc の上位ビット - 16384·1.0
interp_frac  ← acc の下位ビット
result       ← table[lookup_index]
             - interp_frac · (table[lookup_index] - table[lookup_index+1])
```

`PREDICT_NORMALIZE_TABLE` と `PREDICT_SECONDARY_TABLE` (各 33 セル) は
対数的な曲線をテーブル化したもの: 一次テーブルはおおむね 16 ビット数の
$\log_2$ を Q15 にスケールしたもので、二次テーブルは二次正規化のために結果を
整形するコンパニオンです。両者とも [src/tables/gain.rs](../../src/tables/gain.rs)
にあります。

`gain_orchestrate_codec` は対称的な `gain_orchestrate` に対して 1 つ *微妙な*
追加最適化を行います: 正規化関数を呼んでからその結果を伝播させるのではなく、
入力から参照インデックスを **予測** し、`PREDICT_NORMALIZE_TABLE[idx]` および
`PREDICT_NORMALIZE_TABLE[idx+1]` を正規化の *外側* で先取り取得します。
そして関数は事前計算されたこの 2 つのテーブルセルを渡された状態で呼び出され、
二次正規化に対しても予測リフレッシュのフェーズ 2–9 をもう一度再生して二次入力を
覗き見ることで同様のことを行います。このパターンは Rust コードでも、
リバースされた挙動とビット一致を保つために維持されています。

## 4. 履歴セットアップと更新

`gain_history_setup` はテールデケイ前に $g_c$ の候補となる **初期ゲイン** を
生成します。

```
restored_acc = saved_ah ∥ predict_acc_main_low_bits  (40 ビット復元)
shift_count  = history_shift_init - 7                (no-CMPT バックオフ)
acc          = sample · acc32_16(restored_acc) · 2
acc          = shift_acc40(acc, shift_count)
acc          = shift_acc40(acc, -1)
acc          = sat32(acc)
initial_gain = hi16(acc)
```

`gain_history_update` はその後 5 セルの過去ゲイン履歴を更新します。
これは別の `gain_normalize_primary` をもう 1 度実行し (今回は
`PREDICT_NORMALIZE_TABLE` 上の第 3 のルックアップ予測を使用)、結果を更新後
Q15 アキュムレータにシフトし、飽和済みアキュムレータの上位半分と下位半分の
定数倍率混合 (`0x6054 = 24660` を倍率に) として新履歴エントリを計算します。

```
new_gain = ((scratch_lo·24660·2) >> 16 + scratch_hi·24660·2) >> 16
history_out = [new_gain, history[0], history[1], history[2]]
```

リフレッシュ後、`history_out[0..4]` は `[new_gain, prev[0], prev[1], prev[2]]`
で、`history_out[4] = history[4]` は定数スロットとして保持されます — これは
フレーム起動時の定数を運びます (`DecoderState::new` で 0 に初期化)。

## 5. テールデケイ

`gain_tail_decay` は、抑制状態のしきい値が発火した後にゲインを抑える
ソフトリミッタです。サブフレームごとに以下を実行します。

```
threshold_in -= 8                  (TAIL_THRESHOLD_STEP)
if threshold_in >= 0:
    counter ← TAIL_COUNTER_RELOAD = 4
threshold_out = 0

if counter ≤ 0:
    return (デケイなし)              ← (g_p, g_c) を未スケールでパススルー

counter -= 1
t_value = -counter + 3             (TAIL_EXP_BIAS)
scale   = MIN_SYNTH_STATE · 2^t_value     (= 3277 << t_value)
update_base ·= scale          ← スケール後の g_p
initial_gain ·= scale          ← スケール後の g_c
return (デケイ適用済, counter_out)
```

効果としては、しきい値が発火するとゲインが `t_value` (= 末尾に行くほどデケイが強い)
に応じた 2 のべき乗で乗算的に減衰され、最大 `TAIL_COUNTER_RELOAD` = 4
サブフレームの間カウントダウンされます。定数 `MIN_SYNTH_STATE = 3277` は
コーデックの複数箇所で使われる最小状態下限です ([synthesis.md](synthesis.md) 参照)。

定常状態 (suppress なし) では、`threshold_in` は全サブフレームで 0 となり
`counter_out` も 0 のままなので、テールデケイブロックは素通りで
`(g_p, g_c) = (update_base, initial_gain)` となります。

## 6. 配線

サブフレームごとに出力される 4 つのゲイン結果は、[src/lib.rs](../../src/lib.rs)
で次のように接続されます。

```
pitch_gain   = tail.update_base_out   → mix_excitation
fcb_gain     = tail.initial_gain_out  → mix_excitation
threshold    = tail.threshold_out     → 次サブフレームへ
counter      = tail.counter_out       → 次サブフレームへ
history_out[5] → 次サブフレームへ (フレーム末でフレーム永続)
prev_pitch_gain = pitch_gain          → 次サブフレームの compute_pitch_enhance_gain
```

**コミット順序** に注目してください: `prev_pitch_gain` はゲインオーケストレーション
*の後*、しかし次サブフレームの `compute_pitch_enhance_gain` *の前* に更新されます。
これは必須で、ピッチ強調ゲインが synth-control ヒステリシスを通じて前ピッチゲインに
依存するためです ([synthesis.md](synthesis.md#1-synth-control-ピッチ強調のゲーティング) を参照)。

## 7. 初期状態

コーデック初期化時、`gain_history` は `[-17254, -17254, -17254, -17254, 0]` に
セットされ、`gain_threshold = gain_counter = 0` となります。非ゼロの履歴は AR
予測子に小さな負のバイアスを与え、静的なピッチゲイン包絡を近似します — これが
ないと、最初の数フレームは履歴が実値で埋まるまで予測が大きく外れます。
