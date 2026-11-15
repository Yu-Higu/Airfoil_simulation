# Airfoil simulation

3MEの実験流体シミュレーションコードを再現可能な環境へ移し，数値結果の正しさを保ちながら高速化するプロジェクトです．

## 構成

- `src/`：Juliaパッケージ本体
- `test/`：パッケージのテスト
- `project/01_Porting_code_from_fortran/`：パッケージ本体へ未統合のシミュレーションコードと入力データ。既存の回帰テストをルートのテストから実行
- `run.jl`：現在のパッケージ実行入口

旧シミュレーションはパッケージ本体にはまだ統合されていません．一方、既存の回帰テストはルートの `test/runtests.jl` から実行され、GitHub Actions の `Pkg.test()` でも検証されます．旧コードの入力データ `setting.dat` と `clark_y.dat` は保持しています．

## 環境

- Julia 1.13.0
- 依存関係は `Project.toml` と `Manifest.toml` で管理します．

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. run.jl
julia --project=. -e 'using Pkg; Pkg.test()'
```

GitHub Actionsでも同じJulia 1.13系と `Pkg.test()` を使います．

## 実行とテスト

ルートの `run.jl` は現時点ではモデル・ソルバーが未設定であることを表示します．`Pkg.test()` はパッケージの読み込みと旧ソルバーの回帰テストを実行します．回帰テストだけを単独実行する場合：

```sh
julia --project=. project/01_Porting_code_from_fortran/test_3ME_experment_improved.jl
```

旧コードを実行する場合：

```sh
julia --project=. project/01_Porting_code_from_fortran/run.jl
```

旧コードの生成結果はGit管理対象外です．物理モデルの正しさや高速化効果の検証は今後の課題です．

## コードの出典と変更記録

移行するコードごとに，参照元リポジトリ，元課題，参照元commit，ファイル，ライセンス，変更内容を記録します．移行元と対象コードが決まった時点で，この節と関連するテストを更新します．
