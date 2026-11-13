using Test
using Printf
using OffsetArrays

include("3ME_experiment_improved.jl")

const TEST_SETTING = "test_setting.dat"
const TEST_WING = "test_wing.dat"

mktempdir() do test_dir
    cd(test_dir) do
        println("Starting Airfoil Simulation Regression Tests...")

        t_start = time()

        @testset "Airfoil Simulation Regression Test Suite" begin

            # テスト実行環境の準備とダミーデータの作成
            @testset "0. Environment Setup" begin
                mkpath("TEXTDATA")
                @test isdir("TEXTDATA")

                # テストを高速化するため格子数やステップ数を小さくした setting.dat を生成する
                open(TEST_SETTING, "w") do io
                    println(io, "32      ! IG (X grid count)")
                    println(io, "16      ! JG (Y grid count)")
                    println(io, "0       ! IGS")
                    println(io, "100.0   ! RE")
                    println(io, "5.0     ! DEG (AOA)")
                    println(io, "0.001   ! DT")
                    println(io, "2       ! NC (Total steps reduced for fast regression)")
                    println(io, "1       ! NG (Output interval)")
                    println(io, "5.0     ! XL")
                    println(io, "5.0     ! YL")
                    println(io, "'test_wing.dat' ! NWING")
                end
                @test isfile(TEST_SETTING)

                # テスト用のシンプルな翼形状データを生成する
                open(TEST_WING, "w") do io
                    println(io, "Test Wing Header")
                    println(io, "  0.00000  0.00000  0.00000")
                    println(io, "  0.25000  0.02000 -0.02000")
                    println(io, "  0.50000  0.03000 -0.03000")
                    println(io, "  0.75000  0.01000 -0.01000")
                    println(io, "  1.00000  0.00000  0.00000")
                end
                @test isfile(TEST_WING)
            end

            # 設定ファイルおよび翼形状ファイルの読み込み関数のテスト
            @testset "1. Setting and Wing Reader Tests" begin
                p = read_setting(TEST_SETTING)
                @test p.IG == 32
                @test p.JG == 16
                @test p.RE == 100.0
                @test p.DEG == 5.0
                @test p.DT == 0.001
                @test p.NC == 2
                @test p.NWING == "test_wing.dat"

                xt, yt, yb = read_wing(TEST_WING)
                @test length(xt) == 5
                @test xt[1] == 0.0
                @test xt[5] == 1.0
            end

            # 配列のメモリ確保とOffsetArraysのインデックス範囲のテスト
            @testset "2. Allocation & OffsetArray Tests" begin
                p = read_setting(TEST_SETTING)
                S = allocate_sim(p)
                @test S.IG == 32
                @test S.JG == 16
                @test axes(S.U, 1) == 0:(p.IG+1)
                @test axes(S.P, 1) == -4:(p.IG+4)
                @test all(S.BETA .== 1.0)
            end

            # 座標変換、メッシュ生成、および翼面の内外判定（BETA）のテスト
            @testset "3. Mesh Generation & Beta Calculation Tests" begin
                p = read_setting(TEST_SETTING)
                S = allocate_sim(p)
                err = calbeta!(S)
                @test err == 0
                @test S.DX[1] > 0.0
                @test S.DY[1] > 0.0
                @test size(S.PDPHI) == (p.IG + 1, p.JG + 1)
            end

            # 境界条件の適用と、発散値・揚抗力計算（UBC）の数値的健全性テスト
            @testset "4. Boundary Condition & UBC Tests" begin
                p = read_setting(TEST_SETTING)
                S = allocate_sim(p)
                calbeta!(S)

                DIV, CFL, USUM, CL, CD = ubc!(S, 0)
                @test isfinite(DIV)
                @test isfinite(CFL)
                @test isfinite(USUM)
                @test isfinite(CL)
                @test isfinite(CD)
            end

            # 時間積分スキームとポアソン方程式ソルバー（赤黒SOR法）の結合テスト
            @testset "5. Time-Integration & Poisson Solver (Smoke Test)" begin
                p = read_setting(TEST_SETTING)
                S = allocate_sim(p)
                calbeta!(S)

                # 初回ステップの実行（オイラー陽解法）
                fs0!(S)
                ncount_1 = fs2!(S)
                @test ncount_1 > 0
                @test ncount_1 <= 50001

                # 2ステップ目以降の実行（アダムス・バシュフォース法）
                fs1!(S)
                ncount_2 = fs2!(S)
                @test ncount_2 > 0
                @test ncount_2 <= 50001
            end

        end

        t_elapsed = time() - t_start
        println("--------------------------------------------------------")
    end
end