# AirfoilCFD.jl

Fortran 77版 `simulation_airfoil.f` のJulia移植プロジェクトです。

## 実行

Julia 1.9以降を想定しています。

```powershell
cd AirfoilCFD_julia
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. run.jl
```

対話入力の代わりに引数を指定することもできます。

```powershell
julia --project=. run.jl 100000 5.0
```

第3引数は計算ステップ数NC、第4引数は時系列出力間隔NGです。省略時は原文と同じ `NC=10000`, `NG=5` です。

出力は `TEXTDATA` に保存されます。計算終了時に `seed_u.dat`, `seed_p.dat`, `seed_old.dat` を作成します。

## 移植方針

元Fortranの計算式、配列添字、ループ順序、境界条件、SOR反復条件を維持し、Julia化に必要なI/O・配列範囲の修正を行います。翼座標探索で52点配列の53番目を参照する原文の範囲外アクセスは、隣接点探索を1:51に制限しています。

## 注意

本プロジェクトは移植作業中です。Julia実行環境のない作業環境では実行検証できていません。初回実行後は、元Fortran版との数値比較を行ってください。

配列添字の初期確認は次で実行できます。

```powershell
julia --project=. test/runtests.jl
```
