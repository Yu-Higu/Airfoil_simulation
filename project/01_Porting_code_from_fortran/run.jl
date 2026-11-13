# =======================================================================
# run.jl - 2次元翼まわり流れシミュレーション 本番実行スクリプト
# =======================================================================

cd(@__DIR__)

# コアとなるシミュレーションコード（関数定義）を読み込む
include("3ME_experiment_improved.jl")

function run_simulation()
    println("=== Starting Airfoil Simulation (Production) ==/n")
    # 本体側の main() を呼び出す（setting.dat を読み込んで実行）
    main()
    println("/n=== Simulation Completed Successfully ===")
end

run_simulation()