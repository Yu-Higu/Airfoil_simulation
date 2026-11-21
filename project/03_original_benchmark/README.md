# Airfoil simulation ベンチマーク

このディレクトリでは、シミュレーションコードの各処理時間と全体の実行時間を計測し、保存済みの基準データおよび過去の計測結果と比較します。

- 既定の比較対象（候補コード）: [`src/3ME_experiment_improved.jl`](./src/3ME_experiment_improved.jl)
- 既定の基準コード: `../01_Porting_code_from_fortran/3ME_experiment_improved.jl`
- 基準レポートの保存先: `results/baseline/benchmark.txt`

候補コードは絶対パスで指定できるため、ベンチマーク用ディレクトリの外にあるコードも計測できます。計測時の依存パッケージはこのディレクトリの Julia プロジェクトを使用します。

## 1. 初回セットアップ

リポジトリのルートから、ベンチマーク用プロジェクトに移動して依存パッケージを準備します。

```sh
cd project/03_original_benchmark
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## 2. 基準データを一度だけ取得する

最初に既定の基準コードを計測し、レポートを保存します。

```sh
julia --project=. benchmark.jl --record-baseline
```

この実行では、既定の基準コードをベンチマーク用プロジェクトで計測し、`results/baseline/benchmark.txt` に保存します。以後の通常計測ではこのレポートを読み込むため、基準コードは再計測されません。

基準ファイルがすでにある状態で `--record-baseline` を実行すると、上書きせずエラーになります。基準を意図的に更新するときは `--replace-baseline` を指定してください。このオプションは基準の再計測と置き換えを行います。

```sh
julia --project=. benchmark.jl --replace-baseline
```

既存の基準もあとから履歴比較できるよう残したい場合は、更新前に別名でコピーしてください。`benchmark.txt` というファイル名で `results/` 以下に置くと、次回以降の履歴比較対象になります。

```sh
mkdir -p results/baseline-before-update
cp results/baseline/benchmark.txt \
   results/baseline-before-update/benchmark.txt
julia --project=. benchmark.jl --replace-baseline
```

別のコードを基準にする場合は、基準を記録するときにパスを指定します。

```sh
julia --project=. benchmark.jl --record-baseline \
  --baseline-source /path/to/reference.jl
```

`--baseline` も `--baseline-source` と同じ指定として使用できます。基準レポートの保存先は常に `results/baseline/benchmark.txt` です。

## 3. 候補コードを名前付きで計測する
Thread 数を指定する場合，

```sh
JULIA_NUM_THREADS=16 taskset -c 0-15
```

を文頭に追加

基準レポートを作成した後、候補コードを指定して計測します。

```sh
julia --project=. benchmark.jl \
  --source Airfoil_simulation/src/AirfoilSimulation.jl \
  --name cache-layout-v2
```

`--source` を省略した場合は `src/3ME_experiment_improved.jl` が候補になります。`--name` には英数字、ピリオド、`_`、`-` を使用できます。名前は `baseline` にできず、同名の実験結果がすでにある場合は上書きせずエラーになります。

この例では、次の2つのファイルが作成されます。

```text
results/
├── baseline/
│   └── benchmark.txt
└── cache-layout-v2/
    ├── benchmark.txt
    └── comparison.txt
```

- `benchmark.txt`: 候補コードの生の計測結果、環境情報、設定、機械可読の時間データ
- `comparison.txt`: 保存済み基準と過去のベンチマーク結果との比較

`--name` を省略すると `run-yyyymmdd-HHMMSS` 形式の名前が自動で付けられます。実験内容が分かる名前を指定すると、結果を後から探しやすくなります。

## 4. ベンチマーク用ディレクトリ以外から実行する

`--project` にベンチマーク用プロジェクトを、続けて `benchmark.jl` にその絶対パスを指定します。候補コードも絶対パスで指定すれば、現在の作業ディレクトリに依存せず実行できます。

```sh
cd /home/yuya_higuchi/dev/Airfoil_simulation

JULIA_NUM_THREADS=16 BENCHMARK_SECONDS=5 taskset -c 0-15 \
  /home/yuya_higuchi/.juliaup/bin/julia \
  --project=project/03_original_benchmark \
  project/03_original_benchmark/benchmark.jl \
  --source src/AirfoilSimulation.jl \
  --name HOGEHOGE
```

```sh
julia --project=/path/to/Airfoil_simulation/project/03_original_benchmark \
  /path/to/Airfoil_simulation/project/03_original_benchmark/benchmark.jl \
  --source /path/to/my-experiments/3ME_experiment_improved.jl \
  --name my-experiment
```

相対パスで指定した `--source` や `--baseline-source` は、コマンドを実行した時点の作業ディレクトリから解決されます。パスに空白が含まれる場合は、シェルで引用符を付けてください。

## 5. 比較レポートの見方

比較レポートは候補コードの計測後に作成され、以下を比較します。

1. 候補コードと `results/baseline/benchmark.txt`
2. 候補コードと `results/` 以下に保存されている過去のベンチマークレポート（基準レポート自体は除く）

新しい名前付き実験の `benchmark.txt` に加えて、従来の `benchmark-*.txt` 形式のレポートも履歴として読み込みます。比較可能な測定項目について、候補時間、比較対象時間、差分率を表示します。差分率は `(候補時間 / 比較対象時間 - 1) × 100` で、**プラスは候補が遅い、マイナスは候補が速い**ことを示します。

計測時間は Julia のバージョン、CPU、スレッド数、設定値に影響されます。比較レポートには環境情報も記録されるため、履歴と条件が異なる場合はその点を考慮してください。既定の入力データは `../01_Porting_code_from_fortran/setting.dat` と `clark_y.dat` で、実行時に一時ディレクトリへコピーされます。候補コードの置き場所を変えても、これらの入力データの場所は変わりません。

### 同じ条件で計測する

各ベンチマーク対象の計測時間は `BENCHMARK_SECONDS` で変更できます（既定値は5秒）。

```sh
BENCHMARK_SECONDS=15 julia --project=. benchmark.jl \
  --source /path/to/3ME_experiment_improved.jl \
  --name cache-layout-v2-15s
```

この値はベンチマーク対象ごとの目安であり、実行全体の上限時間ではありません。完全なシミュレーションや逐次計測も実行するため、処理には時間がかかる場合があります。基準と候補を比較する際は Julia のバージョン、CPU、スレッド数、入力設定、`BENCHMARK_SECONDS` を揃え、他の重い処理を避けて実行してください。

## 6. 計測内容

`BenchmarkTools` によって、設定ファイル読み込み、配列確保、メッシュ生成、初期・後続ステップの速度予測と圧力計算、境界条件、出力処理、再開用データの読み書き、完全な `main()` 実行を計測します。変更状態を使う計測では、各サンプルの前に同じ入力状態を複製します。出力処理は一時ディレクトリを使い、ディレクトリ準備・後片付けは計測時間に含めません。

加えて、完全なシミュレーションを一度逐次実行し、実際に変化する計算状態で各処理の合計時間と全体に占める割合を計測します。SOR反復回数もレポートに記録します。

## 7. プロファイル

CPU、アロケーション、型安定性の調査には `profile.jl` を使用します。

```sh
julia --project=. profile.jl
```

レポートには主要関数の `@code_warntype`、代表的な1ステップおよび完全な `main()` のCPUプロファイル、1ステップのアロケーションプロファイルが含まれます。完全実行の反復数は `PROFILE_REPEATS`（既定値1）、アロケーションのサンプリング率は `ALLOC_SAMPLE_RATE`（既定値0.001）で変更できます。

```sh
PROFILE_REPEATS=2 ALLOC_SAMPLE_RATE=0.005 \
  julia --project=. profile.jl
```

プロファイルには計測オーバーヘッドがあるため、最適化箇所の特定に使用し、性能差の判定には通常のベンチマークを使用してください。`seedin!` は再開用ユーティリティとして調査対象ですが、既定の `main()` からは呼び出されません。
