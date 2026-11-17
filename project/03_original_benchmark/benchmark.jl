using BenchmarkTools
using Dates

if get(ENV, "BENCHMARK_CHILD", "") != "1"
    include("compare.jl")
    run_comparison()
    exit()
end

include("common.jl")

seconds = parse(Float64, get(ENV, "BENCHMARK_SECONDS", "5"))
seconds > 0 || throw(ArgumentError("BENCHMARK_SECONDS must be greater than zero"))

source_dir = make_workdir()
seed_dir = mktempdir()
try
    p = read_setting(joinpath(source_dir, "setting.dat"))
    p.NC >= 2 || throw(ArgumentError("setting.dat must specify NC >= 2"))
    p.NG >= 1 || throw(ArgumentError("setting.dat must specify NG >= 1"))
    input_state = allocate_sim(p)
    mesh_state = prepare_mesh_state(source_dir, p)
    fs0_input_state = prepare_first_step_input(mesh_state)
    first_state = deepcopy(fs0_input_state)
    fs0!(first_state)
    post_first_state = prepare_post_first_step_state(first_state)

    initial_ubc_state = deepcopy(mesh_state)
    run_states = prepare_run_states(source_dir, p)
    final_state = run_states.final_state
    snapshots = run_states.snapshots
    first_ncount = run_states.first_ncount
    first_snapshot = only(filter(s -> s.num == 2, snapshots))
    middle_snapshot = only(filter(s -> s.num == cld(p.NC, 2), snapshots))
    final_snapshot = only(filter(s -> s.num == p.NC, snapshots))
    cd(seed_dir) do
        seedout(final_state, p.NC, p.NC * p.DT)
    end

    suite = BenchmarkGroup()
    setting_path = joinpath(source_dir, "setting.dat")
    suite["01 Read setting.dat"] =
        @benchmarkable bench_read_setting($setting_path) evals=1
    suite["02 Allocate simulation arrays"] =
        @benchmarkable bench_alloc_sim($p) evals=1
    suite["03 Mesh and airfoil: calbeta!"] =
        @benchmarkable bench_calbeta!(S, $source_dir) setup=(S=deepcopy($input_state)) evals=1
    suite["04 Initial boundaries and diagnostics: ubc! core"] =
        state_benchmark(bench_ubc_initial!, initial_ubc_state)
    suite["05 First velocity predictor: fs0!"] = state_benchmark(fs0!, fs0_input_state)
    suite["06 First pressure solve: fs2! (SOR iterations=$(first_ncount))"] =
        state_benchmark(fs2!, first_state)
    suite["07a Later velocity predictor: fs1! at NUM=2"] =
        state_benchmark(fs1!, first_snapshot.predictor_input)
    suite["07b Later velocity predictor: fs1! at NUM=$(middle_snapshot.num)"] =
        state_benchmark(fs1!, middle_snapshot.predictor_input)
    suite["07c Later velocity predictor: fs1! at NUM=$(final_snapshot.num)"] =
        state_benchmark(fs1!, final_snapshot.predictor_input)
    suite["08a Later pressure solve: fs2! at NUM=2 (SOR iterations=$(first_snapshot.ncount))"] =
        state_benchmark(fs2!, first_snapshot.pressure_input)
    suite["08b Later pressure solve: fs2! at NUM=$(middle_snapshot.num) (SOR iterations=$(middle_snapshot.ncount))"] =
        state_benchmark(fs2!, middle_snapshot.pressure_input)
    suite["08c Later pressure solve: fs2! at NUM=$(final_snapshot.num) (SOR iterations=$(final_snapshot.ncount))"] =
        state_benchmark(fs2!, final_snapshot.pressure_input)
    suite["09a Later boundaries and diagnostics: ubc! at NUM=2"] =
        state_benchmark(bench_ubc_step!, first_snapshot.boundary_input)
    suite["09b Later boundaries and diagnostics: ubc! at NUM=$(middle_snapshot.num)"] =
        state_benchmark(bench_ubc_step!, middle_snapshot.boundary_input)
    suite["09c Later boundaries and diagnostics: ubc! at NUM=$(final_snapshot.num)"] =
        state_benchmark(bench_ubc_step!, final_snapshot.boundary_input)
    suite["10 Periodic boundaries plus field output"] =
        io_benchmark(bench_ubc_with_output!, final_snapshot.boundary_input)
    suite["11 Field output: write_puv"] =
        @benchmarkable bench_write_puv(S, dir, $p.NC) setup=(S=deepcopy($final_state); dir=mktempdir(); mkpath(joinpath(dir, "TEXTDATA"))) teardown=(rm(dir; recursive=true, force=true)) evals=1
    output_num = p.NC % 5 == 0 ? p.NC - 1 : p.NC
    suite["12 Summary output: output.dat"] =
        @benchmarkable bench_output(S, dir, $output_num, $(p.NC * p.DT)) setup=(S=deepcopy($final_state); dir=mktempdir(); mkpath(joinpath(dir, "TEXTDATA"))) teardown=(rm(dir; recursive=true, force=true)) evals=1
    suite["13 Restart output: seedout"] =
        @benchmarkable bench_seedout(S, dir, $(p.NC), $(p.NC * p.DT)) setup=(S=deepcopy($final_state); dir=mktempdir()) teardown=(rm(dir; recursive=true, force=true)) evals=1
    suite["14 Restart input: seedin! (not called by main)"] =
        @benchmarkable bench_seedin!(S, $seed_dir) setup=(S=deepcopy($input_state)) evals=1
    suite["15 Complete run: main()"] =
        @benchmarkable bench_full_run(dir) setup=(dir=make_workdir()) teardown=(rm(dir; recursive=true, force=true)) evals=1

    result_dir = joinpath(@__DIR__, "results")
    mkpath(result_dir)
    result_path = joinpath(result_dir,
                           "benchmark-$(Dates.format(now(), "yyyymmdd-HHMMSS")).txt")
    timed_dir = make_workdir()
    println("Settings: ", p.IG, " x ", p.JG, ", RE=", p.RE, ", DEG=", p.DEG,
            ", NC=", p.NC, ", NG=", p.NG, ", DT=", p.DT)
    println("Benchmarking every stage for approximately ", seconds, " seconds each...")
    results = run(suite; seconds=seconds, evals=1)
    sequential_timing = timed_main(timed_dir)

    open(result_path, "w") do io
        println(io, "Airfoil simulation benchmark using project/01 setting.dat")
        println(io, "Date: ", now())
        println(io, "Source: ", BENCHMARK_SOURCE)
        println(io, "Julia: ", VERSION, "; CPU: ", Sys.CPU_NAME,
                "; threads: ", Threads.nthreads())
        println(io, "Conditions: IG=", p.IG, ", JG=", p.JG, ", RE=", p.RE,
                ", DEG=", p.DEG, ", DT=", p.DT, ", NC=", p.NC, ", NG=", p.NG)
        println(io, "BenchmarkTools seconds per target: ", seconds)
        println(io, "Each mutating sample uses a fresh deepcopy of its fixed stage input.")
        println(io, "All output targets use temporary directories; setup/cleanup is excluded.")
        println(io, "\n=== BenchmarkTools results ===")
        for (name, trial) in results
            println(io, "\n", name)
            show(io, MIME("text/plain"), trial)
            println(io)
        end
        println(io, "\n=== Machine-readable timing metrics ===")
        for (name, trial) in results
            println(io, "METRIC\tbenchmark\t", name, "\t",
                    BenchmarkTools.median(trial).time)
        end
        println(io, "\nMeasured SOR iteration counts in the full simulation:")
        println(io, "First step: ", first_ncount)
        for snapshot in snapshots
            println(io, "Step ", snapshot.num, ": ", snapshot.ncount)
        end
        println(io, "\n=== Per-stage time from one complete, sequential run ===")
        println(io, "These timers follow main() order and the exact configured state trajectory.")
        println(io, "Complete run: ", round(sequential_timing.elapsed / 1e9; digits=6), " s")
        for (name, ns) in sort(collect(sequential_timing.timings); by=last, rev=true)
            println(io, rpad(name, 44), " ",
                    round(ns / 1e9; digits=6), " s  ",
                    round(100 * ns / sequential_timing.elapsed; digits=1), "%")
        end
        println(io, rpad("Time-series, logging, and uninstrumented work", 44), " ",
                round(sequential_timing.other / 1e9; digits=6), " s  ",
                round(100 * sequential_timing.other / sequential_timing.elapsed; digits=1), "%")
        println(io, "\n=== SOR iteration counts from the exact configured run ===")
        println(io, "First step: ", first_ncount)
        for snapshot in snapshots
            println(io, "Step ", snapshot.num, ": ", snapshot.ncount)
        end
        println(io, "Sequential timing run counts: min=", minimum(sequential_timing.ncounts),
                ", max=", maximum(sequential_timing.ncounts))
        println(io, "\n=== Machine-readable sequential metrics ===")
        println(io, "METRIC\tsequential\tComplete run\t", sequential_timing.elapsed)
        for (name, ns) in sequential_timing.timings
            println(io, "METRIC\tsequential\t", name, "\t", ns)
        end
        println(io, "METRIC\tsequential\tOther\t", sequential_timing.other)
    end

    println("\n=== Median stage time / memory (sorted by time) ===")
    ordered = sort(collect(results); by=entry -> BenchmarkTools.median(last(entry)).time,
                   rev=true)
    for (name, trial) in ordered
        estimate = BenchmarkTools.median(trial)
        println(rpad(name, 52), " ",
                lpad(string(round(estimate.time / 1e6; digits=3)), 10), " ms  ",
                estimate.memory, " bytes / ", estimate.allocs, " allocs")
    end
    println("\n=== Sequential stage-time shares ===")
    for (name, ns) in sort(collect(sequential_timing.timings); by=last, rev=true)
        println(rpad(name, 44), " ",
                lpad(string(round(ns / 1e9; digits=4)), 9), " s  ",
                round(100 * ns / sequential_timing.elapsed; digits=1), "%")
    end
    println(rpad("Time-series, logging, and uninstrumented work", 44), " ",
            lpad(string(round(sequential_timing.other / 1e9; digits=4)), 9), " s  ",
            round(100 * sequential_timing.other / sequential_timing.elapsed; digits=1), "%")
    println("Sequential full run: ", round(sequential_timing.elapsed / 1e9; digits=4), " s")
    println("\n=== SOR iteration counts from the exact configured run ===")
    println("First step: ", first_ncount)
    for snapshot in snapshots
        println("Step ", snapshot.num, ": ", snapshot.ncount)
    end
    println("Saved results to: ", result_path)
finally
    rm(source_dir; recursive=true, force=true)
    isdefined(@__MODULE__, :timed_dir) && rm(timed_dir; recursive=true, force=true)
    rm(seed_dir; recursive=true, force=true)
end
