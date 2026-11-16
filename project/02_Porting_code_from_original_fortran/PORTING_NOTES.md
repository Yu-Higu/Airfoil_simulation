# 移植作業メモ（2026-11-13）

## 方針
- Fortran の計算式、反復順序、境界条件、出力頻度を維持する。
- Julia 化に必要な添字・初期化・I/O の修正だけを行う。
- 数値一致・実行確認を行っていないものは未検証として扱う。

## 今回確認した修正点

### 1. MESH 配列の範囲
`CALBETA` は `I=-2:LX+5` および `J=-2:LY+5` で配列を設定し、`DX(-2)=X(-2)-X(-3)` のように `-3` の要素も参照する。元の COMMON 宣言の `-1:LX+4` では不足する。Julia 側の MESH 配列は `-3:LX+5` / `-3:LY+5` とする。

`src/types.jl` の X/XX/DX/DXX と Y/YY/DY/DYY をこの範囲に修正した。

### 2. CALBETA の翼座標探索
翼座標配列は `DIMENSION WNG_*(IGS)` で 52 要素だが、元コードの `DO II=1,IGS` 内で `WNG_*(II+1)` を参照している。最終反復で 53 番目を読むため範囲外となる。隣接点を使う探索ループは `II=1:IGS-1` に修正する必要がある。

### 3. FS2 のSORループ順序
元コードのSOR反復は `DO I=1,IG` が外側、`DO J=1,JG` が内側。これは Gauss-Seidel 型の更新順序に影響するため、Juliaでも `for i in 1:IG; for j in 1:JG` の順序を維持する。以前提示したFS2案のSOR部分はループ順序が逆だったため、そのまま採用しない。

### 4. SEEDIN のドメイン判定
元コードの `IF ((XLD.NE.XL).OR.(XLD.NE.XL))` は YLD を比較していない。意図された判定は `XLD != XL || YLD != YL` と判断し、Julia側ではこの誤記を修正する。

## 未完了
- FS0 / FS1 / FS2 / UBC のJulia実装と原文の最終照合
- CALBETA / OUTPUT / SEEDIN / SEEDOUT のJulia実装
- MAIN相当の実行制御、TEXTDATA出力ディレクトリ作成
- Julia実行環境での構文・実行・数値検証

## 検証上の制約
この作業環境では Julia 実行ファイルを確認できなかったため、現時点でJuliaによる実行確認はしていない。

## Julia実装を追加したサブルーチン
- `FS0!`, `FS1!`, `FS2!`: `src/solver.jl`
- `UBC!`: `src/boundary.jl`
- `CALBETA!`: `src/mesh.jl`
- `OUTPUT!`, `SEEDIN!`, `SEEDOUT!`: `src/io.jl`
- `run_simulation`: `src/AirfoilCFD.jl`

## 追加で見つかった原文の配列範囲問題
`FS0` は `U(0:IG+1,0:JG+1,*)` の初期化で `BETA(I,J,*)` を同じ添字で参照します。また `CALBETA` の `IDY` 計算は `BETA(I+1,J,2)` を参照しますが、COMMON宣言のBETAは `0:IG,0:JG` です。Julia側ではBETAの領域を `0:IG+1,0:JG+1` に拡張し、追加したゴーストセルを自由流値 `1.0` としています。これは複数箇所の範囲外参照を避けるための補正で、ゴーストセルの物理的な意図を自由流と推定した箇所です。

## 現時点の検証状況
全サブルーチンのJuliaコードを作成しましたが、この実行環境にJulia実行ファイルがないため、Juliaの構文検査・実行検証は未実施です。特に以下は実行後に重点確認が必要です。
- `CFDData` のOffsetArray軸とコンストラクタの整合
- Juliaの固定幅出力とFortran書式の完全一致
- 翼座標の範囲外探索修正によるBETA分布
- `FS0` ゴーストセル初期値の妥当性
- 1ステップおよび複数ステップのFortran/Jデータ比較

## 原文のビルド確認
作業環境のgfortranで元Fortranソースのコンパイルは成功しました（`-O0 -ffixed-line-length-none`）。ただしコンパイラは、元ソースにMESH配列、BETA配列の複数の範囲外参照があると警告しています。Julia側では次のように範囲を拡張・補正しました。
- MESH: `-3:LX+5` / `-3:LY+5`
- BETA: `0:IG+1` / `0:JG+1` とし、拡張したゴーストセルを1.0で初期化
- 翼表面探索: `II=1:IGS-1`

このgfortranコンパイルは元Fortranのコンパイル確認であり、Julia版の実行確認ではありません。
