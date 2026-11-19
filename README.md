# Airfoil Simulation

2次元翼周りの非圧縮性粘性流体シミュレーションコード（有限差分法・スタガード格子）を Julia で実装・高速化するプロジェクトです．

元の Fortran 実装および Julia 移植初期版から，**Julia 標準配列（1-based Array）への移行** および **Polyester によるマルチスレッド並列化** を行い，計算結果の数値的再現性を担保しつつ高速化を行っています．

---

## 主な機能・特徴

- **有限差分法・スタガード格子** による2次元ナビエ・ストークス方程式の数値解析
- **時間積分スキーム**:
  - 初回ステップ: オイラー陽解法 (`fs0!`)
  - 2ステップ目以降: アダムス・バシュフォース法 (`fs1!`)
- **圧力ポアソン方程式ソルバー**: 赤黒SOR法 (`sor_row`, `fs2!`)
- **マルチスレッド並列化**: `Polyester.jl` の `@batch` による軽量スレッド並列化と `@fastmath` 最適化
- **Julia パッケージおよびスクリプト両対応**:
  - パッケージモジュール `AirfoilSimulation` として呼び出し可能
  - CLI から `julia src/AirfoilSimulation.jl` として単体実行可能

---

## ディレクトリ構成

- `src/AirfoilSimulation.jl`: シミュレーション本体（高速化・標準配列版）
- `test/runtests.jl`: 回帰テストスイート（`Pkg.test()` で実行）
- `setting.dat`: 計算条件設定ファイル（格子数，レイノルズ数，迎え角，時間刻み幅等）
- `clark_y.dat`: 翼形状データファイル（Clark Y）
- `project/`:
  - `01_Porting_code_from_fortran/`: Fortran から移植した初期の改善版コード
  - `03_original_benchmark/`: 各ステージの実行時間・メモリ使用量を計測・比較するベンチマーク環境
  - `04_Adapting_syntax_for_Julia/`: Juliaらしい構文へ適応したリファレンス実装（OffsetArrays版）

---

## 環境構築

- **推奨環境**: Julia 1.13.0
- 依存関係はルートの `Project.toml` および `Manifest.toml` で管理されています．

```sh
# パッケージの依存関係をセットアップ
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

---

## 実行方法

### 1. スクリプトとして直接実行

カレントディレクトリに `setting.dat` および翼形状データ（例: `clark_y.dat`）を配置した状態で実行します．

```sh
julia --project=. src/AirfoilSimulation.jl
```

マルチスレッドを指定して実行する場合（例: 16スレッド）:
```sh
julia -t 16 --project=. src/AirfoilSimulation.jl
```

### 2. Julia REPL / 外部スクリプトからの呼び出し

```julia
using AirfoilSimulation

# setting.dat を読み込んでシミュレーションを実行
AirfoilSimulation.main()

# 設定ファイルを指定して実行
AirfoilSimulation.run_simulation("setting.dat")
```

### 出力ファイル
実行するとカレントディレクトリに以下のファイルが生成されます．
- `TEXTDATA/time-series.dat`: 各ステップの計算状況（発散値 DIV, CFL数, 平均流速, 揚力係数 $C_L$, 抗力係数 $C_D$, ポアソン反復回数 NCOUNT）
- `TEXTDATA/output.dat`: 格子情報および翼表面圧力分布サマリー
- `TEXTDATA/puv_field0000.dat` など: 格子点ごとの流速・圧力場データ
- `seed_u.dat`, `seed_p.dat`, `seed_old.dat`: 計算再開（リスタート）用倍精度データ

---

## テストの実行

回帰テストスイートを実行し，設定読み込み，配列確保，メッシュ生成，境界条件計算，時間積分ソルバーの健全性を検証します．

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```
または
```sh
julia --project=. test/runtests.jl
```

---

## ベンチマーク・プロファイル

各関数の実行時間計測やベースラインとの比較は `project/03_original_benchmark` で行えます．詳細は [`project/03_original_benchmark/README.md`](project/03_original_benchmark/README.md) を参照してください．

```sh
# src/AirfoilSimulation.jl の性能を計測・既存基準と比較
julia --project=project/03_original_benchmark project/03_original_benchmark/benchmark.jl \
  --source=src/AirfoilSimulation.jl \
  --name="StandardArray_Polyester"
```
