# デコーダ状態

このデコーダは、ストリーミングでフレーム単位の状態を持つマシンです: 各 18
バイトのフレームは、前フレームから引き継がれるおよそ 1 フレーム分の状態に
依存します。本ドキュメントはその状態を列挙し、各要素がフレーム・ブロック・
サブフレーム間でどのように流れるかを説明します。
ソース: [src/state.rs](../../src/state.rs), [src/lib.rs](../../src/lib.rs)。

## 1. `DecoderState` 構造体

```rust
pub struct DecoderState {
    pub lsf_history:               LsfHistory,           // 3 スロットの LSF ローリング履歴
    pub prev_frame_lsp:            [i16; 10],            // 前フレームの LSP (Q15)
    pub prev_block_coeffs:         [i16; 10],            // 前ブロックの LPC 作業係数
    pub past_excitation:           [i16; 320],           // スライディング励振バッファ
    pub lpc_synth_history:         [i16; 10],            // 合成出力の最後 10 サンプル
    pub pitch_state_block0:        PitchLagState,        // sub 0/1 の prev_lag
    pub pitch_state_block1:        PitchLagState,        // sub 2/3 の prev_lag
    pub gain_history:              [i16; 5],             // 5 セルの過去ゲインリング
    pub gain_threshold:            i16,                  // テールデケイのしきい値
    pub gain_counter:              i16,                  // テールデケイのカウントダウン
    pub postfilter_history:        [i16; 6],             // 前向き+後ろ向き 4-tap シフトレジスタ
    pub postfilter_delay:          [i16; 12],            // 前向き+後ろ向き 3 セクションのディレイ
    pub prev_lag:                  i16,                  // 直近ブロックで選択された lag (synth control)
    pub prev_pitch_gain:           i16,                  // 直近サブフレームの g_p (synth control)
    pub response_shaper_buffer:    [i16; 35],            // パス間で持ち越す dword スクラッチ
    pub prev_lpc_stability_flag:   i16,                  // 0/1 — block-0 オーバライドのエッジ判定用
}
```

総計の状態は約 800 バイト — 安価にコピーできるほど小さく、コーデックを厳密に
*ステートフル* にするほど大きい量です。

## 2. 初期値

```rust
DecoderState::new() = {
    lsf_history          = LsfHistory::new()       // 3 スロット全てが INIT_LSF_TEMPLATE
    prev_frame_lsp       = INIT_LSP_TEMPLATE
    prev_block_coeffs    = [0; 10]
    past_excitation      = [0; 320]
    lpc_synth_history    = [0; 10]
    pitch_state_block0   = PitchLagState::default()  // prev_lag = 0
    pitch_state_block1   = PitchLagState::default()
    gain_history         = [-17254, -17254, -17254, -17254, 0]
    gain_threshold       = 0
    gain_counter         = 0
    postfilter_history   = [0; 6]
    postfilter_delay     = [0; 12]
    prev_lag             = 60
    prev_pitch_gain      = 3277
    response_shaper_buffer = [0; 35]
    prev_lpc_stability_flag = 0
}
```

非ゼロの種値のうち、特に注目すべきものを挙げます。

- `INIT_LSF_TEMPLATE` は均等間隔の 10-LSF ベクトル
  `[2340, 4679, ..., 23396]` で、すなわち Q15 の `[0, 0.7]` を等分割したもの。
  デコーダ起動時、LSF 履歴はこのテンプレートを 3 つコピーで保持するので、
  最初のフレームで `predictive_combine` がゴミと積算することがありません。
- `INIT_LSP_TEMPLATE` はその等分割の cos (オフラインで計算しテーブル化)。
  最初のフレームの `prev_frame_lsp` はこのテンプレートです。
- `gain_history` は `-17254` を 4 回繰り返してゼロパディングした 5 番目の
  スロットでシードされます。これは平均的な過去ゲイン予測状態を近似する
  小さな負のバイアスで、これがないと最初の数フレームの予測が大幅に外れます。
- `prev_lag = 60` と `prev_pitch_gain = 3277` は、起動時に synth-control
  ヒステリシスを安定した中間状態にするための静的な値です。

## 3. フレーム内の状態の流れ

フレーム内では、状態は次の順で更新されます。

1. **LSF 履歴** (`lsf_history.splice_new_vector`) — `decode_lsf_direct_mode`
   *の最中*、安定化の前に、生のスクラッチベクトル (予測出力ではない) で
   コミットされる。
2. **前フレームの LSP** (`prev_frame_lsp`) — ステージ 1 の最後 (`lsf_to_lsp`
   の後) でコミットされ、*現*フレームの LSP が *次* フレームの
   `prev_frame_lsp` になる。
3. **`response_shaper_buffer`** — 3 回 (パスごとに 1 回) 書き込まれ、毎回
   新しい値がコミットされる。フレーム内の 3 パスの順序は重要。
4. **`prev_lpc_stability_flag`** — 安定性チェックの後にコミット。
5. **サブフレームごと** (i = 0..3):
   - `past_excitation[write_offset + n]` ← `e[n]`
   - `lpc_synth_history` ← `s_hat` の最後 10 サンプル
   - `prev_pitch_gain` ← `gain_out.pitch_gain`
6. **ブロックごと** (sub 1 と sub 3 の後):
   - `past_excitation[160..320]` を `[0..160]` にシフトし、`[160..320]` を
     クリア
   - `prev_lag` ← このブロックの最終サブフレーム lag
7. **フレームごと** (最終サブフレームの後):
   - `gain_history`, `gain_threshold`, `gain_counter` ← ゲインオーケストレータの
     出力から
   - `postfilter_history`, `postfilter_delay` ← `postfilter_apply` から
     (in-place)

ピッチ状態 `pitch_state_block0` と `pitch_state_block1` は `decode_subframe_lag`
の呼び出し中に変更され (sub-1 / sub-3 の差分復号用にブロックごとの `prev_lag`
を保持)、明示的にはコミットされません。

## 4. リセットフレーム

コーデックはインバンドのリセットフレーム
(`(frame[17] & 0x0F) == 0x08` — [architecture.md](architecture.md#33-リセットフレーム) を参照)
を認識し、`DecoderState::reset()` を呼び出すことで反応します。これは
`*self = DecoderState::new()` です。リセット後はデコーダはちょうどインスタンス化
されたばかりのように振る舞います。リセットフレーム自体は音声を生成しません —
`Codec::decode_frame` は `None` を返します。

## 5. バッファの所有権

スクラッチ *のように見える* がコール間で永続化される状態フィールドがいくつかあります。

- **`response_shaper_buffer`**: この 35 セルの `[i16]` は LPC 解析の dword
  スクラッチバッファ。フレームごとに 3 回 (安定性パス、block-0 パス、block-1
  パス) 書き込まれ、得られたバッファは *意図的に* 次フレームへ持ち越され、
  そこで次の安定性パスの初期バッファになります。これをフレームごとのゼロ
  初期化に置き換えると、ビット一致が崩れます。
- **`prev_block_coeffs`**: 構造体に宣言されていますが、現在のデコーダはより
  単純なブロック平均規則に従うため、このフィールドを読み取りません。
  リバースエンジニアリングの表面が露出していて、しかしランタイムの復号経路
  では使われない代替的なサブフレーム LSP 補間経路 (`subframe_lsp_interp`) との
  互換性のために保持されています。

## 6. サブフレーム / ブロック / フレームの関係まとめ

| 量                          | 周期         | 更新タイミング                                  |
| --------------------------- | ------------ | ----------------------------------------------- |
| `lsf_history`               | フレーム毎    | LSF 復号中 (安定化前)                            |
| `prev_frame_lsp`            | フレーム毎    | `lsf_to_lsp` の後                                |
| `response_shaper_buffer`    | パス毎        | 各 `build_autocorrelation_with_state` の後       |
| `prev_lpc_stability_flag`   | フレーム毎    | 安定性チェックの後                                |
| `past_excitation`           | サブ毎        | `pitch_adaptive_codebook` 内 + 混合の後          |
| `past_excitation` (シフト)  | ブロック毎    | サブフレーム 1 とサブフレーム 3 の後              |
| `lpc_synth_history`         | サブ毎        | `lpc_synthesis_filter` の後                      |
| `prev_pitch_gain`           | サブ毎        | `gain_orchestrate_codec` の後                    |
| `prev_lag`                  | ブロック毎    | サブフレーム 1 とサブフレーム 3 の後              |
| `gain_history`              | サブ毎        | `gain_history_update` 内、最終コミットはフレーム毎 |
| `gain_threshold/counter`    | サブ毎        | `gain_tail_decay` 内、最終コミットはフレーム毎    |
| `postfilter_history/delay`  | ハーフフレーム毎 | `postfilter_apply` 内                           |
| `pitch_state_block{0,1}`    | サブ毎        | `decode_subframe_lag` 内                         |
