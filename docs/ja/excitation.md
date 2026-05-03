# 励振: ピッチコードブックと固定コードブック

各 80 サンプルのサブフレームの励振信号 `e[n]` は、2 つの寄与の線形結合です。

```
e[n] = pitch_gain · v[n] + fcb_gain · c[n]      (合計は mix_excitation でスケール)
```

- `v[n]` は **適応コードブック** の出力 — 過去励振バッファから引き出される
  分数ピッチコピー (長期予測)。
- `c[n]` は **固定コードブック** の出力 — 疎な代数的符号語で、必要に応じて
  ピッチ周期性 IIR で強調されたもの。

本ドキュメントは `v` と `c` の生成を扱います。ゲイン量子化 (`pitch_gain` と
`fcb_gain` の計算) は [gain.md](gain.md) のテーマです。
ソース: [src/pitch.rs](../../src/pitch.rs),
[src/fcb.rs](../../src/fcb.rs), [src/tables/pitch.rs](../../src/tables/pitch.rs),
[src/tables/fcb_short.rs](../../src/tables/fcb_short.rs),
[src/tables/fcb_main.rs](../../src/tables/fcb_main.rs)。

## 1. ピッチ lag 復号

各ブロックは、最初のサブフレームに対しては絶対ピッチ lag を、2 番目の
サブフレームに対しては差分 lag を運びます。ビットストリームの配置は次のとおり。

| サブフレーム | フィールド | 幅    | 符号化              |
| ------------ | ---------- | ----: | ------------------- |
| sub 0        | F[2]       | 8     | 絶対値, block 0     |
| sub 1        | F[5]       | 5     | 差分, block 0       |
| sub 2        | F[8]       | 8     | 絶対値, block 1     |
| sub 3        | F[11]      | 5     | 差分, block 1       |

`src/pitch.rs` の `PitchLagState` 構造体は、ブロック内の `prev_lag` を保持します
(block 0 の sub 0 → sub 1, block 1 の sub 2 → sub 3)。block 0 の状態は
block 1 の状態とは独立です。

### 1.1 絶対値復号

`decode_lag_absolute` (8 ビット入力):

```
if phase < 197:
    lag     = (phase + 2) / 3 + 19            // 整数 lag, 範囲 20..85
    sub_lag = phase - lag · 3 + 58            // 分数オフセット, 符号付き
else:
    lag     = phase - 112                     // 85..143, 整数のみ
    sub_lag = 0
```

2 つの領域があります: 小さな lag (`< 85`) は 1/3 サンプルの分数分解能で
コード化され (したがって `(lag, sub_lag)` の組が `1/3` 刻みの分数ピッチを表す)、
大きな lag (`≥ 85`) は整数のみです。整数の範囲は 20..143 サンプル
(8 kHz で約 56 Hz から 400 Hz) をカバーします。

### 1.2 差分復号

`decode_lag_differential` (5 ビット入力 + 前サブフレームの lag):

```
anchor = max(prev_lag - 5, 20)
if anchor + 9 > 143: anchor = 134            // 上限クランプ
q  = (phase + 2) / 3
t1 = q - 1
lag     = anchor + t1                        // anchor-5..anchor+9
sub_lag = phase - 2 - t1 · 3
```

差分 lag は `[anchor-1, anchor+10]` の整数値にまたがり、これも 1/3 サンプルの
分数分解能を持ちます。アンカーは前 lag を 5 だけ下げたもの (もしくは上端を
超える場合は 134 に固定) — 一般的な CELP デルタピッチコーダと同じ形状です。

### 1.3 分数オフセットから補間タップへ

`decode_lag_fract(sub_lag)` は符号付き `sub_lag` を以下に変換します。

- **fract index** `∈ {0, 1, 2}` — 補間フィルタを選択
- **lag adjustment** `lag_adjust ∈ {0, 1}` — 整数 lag に加える補正

```
a = -sub_lag
lag_adjust = (a < 0) ? 1 : 0
a += lag_adjust · 3
return (a, lag_adjust)
```

最終的に `(integer_lag, sub_lag)` は半開区間の分数ピッチと、3 つの 10-tap
補間フィルタのうち 1 つを選択する位相インデックス `0..2` に対応付けられます。

## 2. 適応コードブック (分数ピッチ補間)

`pitch_adaptive_codebook` は、過去励振サンプル間を 10-tap ポリフェーズ
フィルタで補間して `v[n]` の 80 サンプルを合成します。補間フィルタは
[src/tables/pitch.rs](../../src/tables/pitch.rs) にあります。

```
INTERP_FILTER[3][10] (Q15):
    fract=0:  長く緩やかに減衰するローブ (フル分解能の中心フィルタ)
    fract=1:  +1/3 サンプルシフトのフィルタ
    fract=2:  +2/3 サンプルシフトのフィルタ
```

前向きフィルタ `h` は `INTERP_FILTER[fract]` をそのまま使用します。
後ろ向きフィルタ `h_rev` は `fract == 0` の場合は `[h[1..10], 0]`、
そうでない場合は *他方* の分数フィルタ (`fract=1 → INTERP_FILTER[2]`,
`fract=2 → INTERP_FILTER[1]`) を使用します。これにより整数 lag タップを
中心とする対称な両側畳み込みが得られます。

各出力サンプル `n ∈ 0..80` について:

```
acc = Σ_{k=0..9} past[base + n - k] · h[k]
    + Σ_{k=0..9} past[base + n + 1 + k] · h_rev[k]
v[n] = clamp_i16( (acc + (1 << 14)) >> 15 )

past[write_offset + n] = v[n]   ← 同じ呼び出し内での自己フィードバック
```

`base` は `write_offset − effective_lag` で、
`effective_lag = lag + lag_adjust` です。重要なのはループ内の
**自己フィードバック**: 各出力サンプルが `write_offset + n` の位置で過去励振
バッファに書き戻されるので、同じサブフレームの後続サンプルがそれを参照できる、
という点です。これにより非常に短いピッチ lag (`lag < 80`) は再帰的なピッチ
予測子として振る舞い、lag がサブフレーム長より短くても周期構造を生成します。

`write_offset` は `state::PAST_EXCITATION_OUTPUT_BASE + sub_in_block`
(すなわちブロックの最初のサブフレームでは 160、2 番目では 240) です。
過去励振バッファのレイアウトは [state.md](state.md) で説明します。

呼び出し周辺の境界チェックは `base ≥ 9` および `base + 80 + 10` が範囲内に
収まることを保証します。そうでなければ `v` はゼロのまま残ります
(ピッチ状態がまだ完全に充填されていないときの優雅な退化)。

## 3. 固定コードブック

FCB は **代数的** です: 各 16 ビットの FCB インデックスは小さな単位振幅パルス
位置の集合に復号され、固定パルス振幅にスケーリングされ、必要に応じて IIR で
周期化されます。インデックスの大きさに応じて 2 つの経路が選択されます。

```
fcb_index in {short path | main path 0 | 1 | 2 | 3 | 4}
            ↑
            └ DISPATCH_THRESHOLDS テーブルが経路 (= "lag class") を選択
```

### 3.1 ディスパッチと lag クラス

`fcb_dispatch_lag_class(fcb_index)` ([src/fcb.rs](../../src/fcb.rs)) は
16 ビットインデックスを 6 エントリの `DISPATCH_THRESHOLDS` で分類します。

```
DISPATCH_THRESHOLDS = [65432, 49048, 45568, 39168, 32768, 16384]
```

対応関係:

| 範囲 (符号なし)            | lag_class | 経路        |
| -------------------------- | :-------: | ----------- |
| `[0, 16384)`               | 0         | main, k=0   |
| `[16384, 32768)`           | 1         | main, k=1   |
| `[32768, 39168)`           | 2         | main, k=2   |
| `[39168, 45568)`           | 3         | main, k=3   |
| `[45568, 49048)`           | 4         | main, k=4   |
| `[49048, 65432)`           | 5         | **short**   |
| `[65432, ∞)` (クランプ済) | 5         | short       |

つまりフォーマットは *6 つの lag クラス* を予約しています — 段階的に狭くなる
バイアスオフセットを持つメイン代数的コードブック用に 5 つ、よりタイトな
"short" コードブック用に 1 つ。`clamped_fcb_index` は `THRESHOLDS[0] - 1` の
上限にクランプされた後のインデックスです。

### 3.2 short 経路 (κ = 5)

lag クラス 5 では、`fcb_short_path` が符号語を 5 パルストラック 2 本
(Track A と Track B) 用の 14 ビットサブインデックス対に復号します。

```
delta = (clamped_index - CODEWORD_BIAS) mod 2^16
track_a = delta >> 7     // 7 ビット
track_b = delta & 0x7F   // 7 ビット
```

各トラック行は、`(sign, |pos| - 1)` としてエンコードされた *符号付き* パルス
位置の 5 要素 `i16` 配列です。

```
TRACK_A[track_a] = [b0, b1, b2, b3, b4]    // 各 |b| - 1 は 0..79 の位置
TRACK_B[track_b] = [b0, b1, b2, b3, b4]
```

`fcb_short_pulse_synthesis` は各トラックの位置に振幅 **8192**
(Q15 ≈ 0.25) の単位パルスを置きます。2 トラックが同じ位置に来ると、
そこでパルスは重畳します。short コードブック全体は `128 · 128 = 16,384`
符号語にまたがります。

整数 lag が完全なサブフレーム長より短い場合 (`lag < 80`)、`pitch_enhance` が
適用されます (§3.4 参照)。

### 3.3 main 経路 (κ ∈ 0..5)

main lag クラス 0..4 では、符号語は *位置コードブック* から選択された位置に
配置される 2 個または 3 個のパルス列に復号されます。位置コードブックは小さな
「再帰出力」のシーケンスでインデックスされます。復号パイプラインは:

```
fcb_main_pulse_decode → (bit_count, bit_decomposition[3], recursion_outputs[3])
fcb_main_pulse_template → ブロックあたり 80 サンプルのパルステンプレート (振幅)
fcb_main_codebook_synth → 最終的な 80 サンプルの c[n]
```

クラスごとのパラメータ (`BIT_COUNT`, `BIAS`, `T_MPYA_TABLE`, `T_MPY_TABLE`,
`PULSE_SEEDS`, `PULSE_DATA`, `PULSE_POSITION_CODEBOOK`) は
[src/tables/fcb_main.rs](../../src/tables/fcb_main.rs) にあります。
復号は次のように動作します。

1. **実効インデックス** — クラスごとの `BIAS` を引く (各 lag クラスの窓内で
   インデックス 0 がそのクラスのサブコードブックの 0 番から始まるように)。
2. **ビット分解** — 実効インデックスの下位 `bit_count` ビットを個々の符号
   ビット `bit_decomposition[k] ∈ {0, 1}` に分割。
3. **再帰** — インデックスの上位ビットが小さな再帰を駆動し、各ステップ `k` で
   以下を計算:

   ```
   m1                = (T_MPYA[k] · a_hi) >> 17
   recursion_outputs[k] = a_hi - T_MPY[k] · m1
   a_hi              = m1
   ```

   これは上位ビットを `T_MPY[k]` を底とする表現に展開するもの — 1 パルスにつき
   `T_MPY[k]` 通りの位置から 1 つを選ぶ固定小数点の商 / 剰余カスケードです。
4. **パルステンプレート** — 5 行のテンプレート (各 20 サンプルの後に 60 個の
   ゼロが続いて 80 サンプルブロックを埋める) が `PULSE_DATA` から、
   `PULSE_SEEDS[lag_class] · 20` のオフセットを起点に読み出されます。
   ブロック数は `bit_count + 1`。テンプレートはこの lag クラスの **パルスの振幅**
   を与えます (これらは単位パルスではなく重み付きで減衰する包絡。
   [src/tables/fcb_main.rs](../../src/tables/fcb_main.rs) を参照)。
5. **コードブック合成** — `bit_count` 個のパルスそれぞれについて、
   `PULSE_POSITION_CODEBOOK[lag_class · 120 + k · 40 + position]` が 80 サンプル
   サブフレーム内の *開始位置* を与えます。ここで `position` は
   `recursion_outputs[bit_count - 1 - k]` です。`pulse_offset = k · 80` における
   パルステンプレートが、`bit_decomposition[k] == 1` ならそこから加算、
   `== 0` なら減算で `c[]` に足し込まれます。サンプルごとに飽和が適用されます。

各 main lag クラスのコードブックサイズは、隣接する `DISPATCH_THRESHOLDS` の
ギャップ (= そのクラスにルーティングされる固有の `clamped_fcb_index` の数) です:

| `lag_class` | `BIT_COUNT` | パルスクラスタ | コードブックサイズ |
| :---------: | :---------: | :------------: | -----------------: |
| 0           | 3           | 3              | 16384              |
| 1           | 3           | 3              | 16384              |
| 2           | 2           | 2              |  6400              |
| 3           | 2           | 2              |  6400              |
| 4           | 2           | 2              |  3480              |

コードブック合成の後、`lag < 80` であれば再び `pitch_enhance` が適用されます。

### 3.4 ピッチ周期性強調 (`pitch_enhance`)

short と main いずれの経路でも、整数ピッチ lag が完全なサブフレームより小さい
(`lag < 80`) 場合、FCB 出力 `c[n]` は `pitch_enhance` で後処理されます。
役割は **ピッチの周期を FCB 励振に注入する** ことであり、これにより有声フレームでも
FCB 符号語自体は非周期的でありながら、周期的な FCB 成分が保たれます。

実装は単一タップの再帰フィルタで、適応コードブックと同じ 10-tap ピッチ補間
フィルタを再利用します。

```
lag_internal = lag - 10 + lag_adjust
buf[0..lag_internal] をスクラッチライン (オフセット 21) にコピー
for n in lag_internal..80:
    acc = Σ_k scratch[centre - k] · fwd[k]
        + Σ_k scratch[centre + 1 + k] · rev[k]
    new = sat16( (buf[n] << 16 + (acc · 2 · gain) >> 14 + 0x8000) >> 16 )
    buf[n]    = new
    scratch[…] = new       ← 新値はスクラッチラインにフィードバック
```

`gain` は `compute_pitch_enhance_gain` が計算する **ピッチ強調ゲイン** で、
lag クラスでディスパッチされます。

- lag クラス `{0, 3, 4}`: gain は常に **有効** (`= 16384 = Q15 0.5`)。
- lag クラス `{1, 2, 5}`: gain は `update_synth_control` の出力
  ([synthesis.md](synthesis.md#1-synth-control-ピッチ強調のゲーティング) を参照)。
  これは現ピッチ lag, 前ピッチ lag, LPC 解析からの第 1 反射係数 $k_1$,
  および前ピッチゲインを調べて有効化するか判断します。

`gain = 0` のとき IIR は素通り (再帰の寄与がなくなる) となり、FCB 出力は
代数的符号語のまま残ります。

## 4. 励振の混合

`v` と `c` の両方が形成された後、`src/synth.rs` の [`mix_excitation`][synth] は
以下を生成します。

```
e[n] = sat16(((fcb_gain · c[n]) << 3 + (pitch_gain · v[n]) << 2 + 32768) >> 16)
```

— 32 ビットアキュムレータが `8 · fcb_gain · c[n] + 4 · pitch_gain · v[n]`
に Q15 の丸めバイアスを加えたものを保持し、その後 i32 飽和とともに右シフトで
Q15 にします。8 と 4 の係数はゲインパイプラインの Q フォーマット選択を
補正する固定スケール係数です。

混合された励振 `e[n]` は `write_offset + n` の位置で過去励振バッファに
書き戻されます (一時的な `v[n]` 自己フィードバックコピーを上書き)。
これは同時に LPC 合成フィルタにも渡されます ([synthesis.md](synthesis.md) を参照)。

[synth]: ../../src/synth.rs
