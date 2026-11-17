const BENCHMARK_SOURCE = abspath(get(ENV, "BENCHMARK_SOURCE",
                                     joinpath(@__DIR__, "src",
                                              "3ME_experiment_improved.jl")))
isfile(BENCHMARK_SOURCE) ||
    throw(ArgumentError("benchmark source file does not exist: $BENCHMARK_SOURCE"))
include(BENCHMARK_SOURCE)

const SOURCE_DIR = joinpath(@__DIR__, "..", "01_Porting_code_from_fortran")
const SETTING_SOURCE = joinpath(SOURCE_DIR, "setting.dat")
const WING_SOURCE = joinpath(SOURCE_DIR, "clark_y.dat")

function make_workdir()
    dir = mktempdir()
    cp(SETTING_SOURCE, joinpath(dir, "setting.dat"))
    cp(WING_SOURCE, joinpath(dir, "clark_y.dat"))
    return dir
end

function quiet(f)
    return redirect_stdout(devnull) do
        f()
    end
end

function prepare_mesh_state(dir, p)
    S = allocate_sim(p)
    cd(dir) do
        quiet() do
            calbeta!(S)
        end
    end
    return S
end

function prepare_first_step_input(mesh_state)
    S = deepcopy(mesh_state)
    ubc!(S, 1)
    return S
end

function prepare_first_step_state(mesh_state)
    S = prepare_first_step_input(mesh_state)
    fs0!(S)
    return S
end

function prepare_post_first_step_state(first_state)
    S = deepcopy(first_state)
    fs2!(S)
    ubc!(S, 1)
    return S
end

function prepare_run_states(dir, p)
    S = allocate_sim(p)
    checkpoints = unique((2, cld(p.NC, 2), p.NC))
    snapshots = NamedTuple[]
    first_ncount = 0
    cd(dir) do
        quiet() do
            mkpath("TEXTDATA")
            calbeta!(S)
            ubc!(S, 0)
            for NUM in 1:p.NC
                capture = NUM in checkpoints
                capture && (predictor_input = deepcopy(S))
                NUM == 1 ? fs0!(S) : fs1!(S)
                capture && (pressure_input = deepcopy(S))
                ncount = fs2!(S)
                capture && (boundary_input = deepcopy(S))
                ubc!(S, NUM)
                if capture
                    push!(snapshots, (num=NUM, predictor_input=predictor_input,
                                      pressure_input=pressure_input,
                                      boundary_input=boundary_input,
                                      ncount=ncount))
                end
                NUM == 1 && (first_ncount = ncount)
            end
        end
    end
    return (final_state=S, snapshots=snapshots, first_ncount=first_ncount)
end

function benchmark_main(dir)
    return cd(dir) do
        quiet() do
            main()
        end
    end
end

function timed_stage!(f, timings::Dict{String,UInt64}, name::String)
    start = time_ns()
    value = f()
    elapsed = UInt64(time_ns() - start)
    timings[name] = get(timings, name, UInt64(0)) + elapsed
    return value
end

function timed_main(dir)
    return cd(dir) do
        quiet() do
            timings = Dict{String,UInt64}()
            ncounts = Int[]
            run_start = time_ns()

            p = timed_stage!(timings, "Read setting.dat") do
                read_setting("setting.dat")
            end
            S = timed_stage!(timings, "Allocate arrays") do
                allocate_sim(p)
            end
            DT, NC, NG = S.DT, S.NC, S.NG

            NUM = 0
            TIME = 0.0
            AVETIM = 0.0
            ierr = timed_stage!(timings, "Mesh generation: calbeta!") do
                calbeta!(S)
            end
            ierr == 0 || error("calbeta! returned $ierr")

            timed_stage!(timings, "Create TEXTDATA directory") do
                mkpath("TEXTDATA")
            end
            div, cfl, usum, cl, cd = timed_stage!(timings, "Initial UBC + field output") do
                ubc!(S, NUM)
            end

            NUM += 1
            TIME += DT
            AVETIM += DT
            if NUM == 1
                timed_stage!(timings, "First predictor: fs0!") do
                    fs0!(S)
                end
            else
                timed_stage!(timings, "Later predictor: fs1!") do
                    fs1!(S)
                end
            end
            ncount = timed_stage!(timings, "Pressure solve: fs2!") do
                fs2!(S)
            end
            push!(ncounts, ncount)
            div, cfl, usum, cl, cd = timed_stage!(timings, "Step UBC + periodic field output") do
                ubc!(S, NUM)
            end

            open("TEXTDATA/time-series.dat", "w") do ts
                println(ts, TS_HEADER)
                println(TS_HEADER)
                line = ts_line(NUM, TIME, div, cfl, usum, cl, cd, ncount)
                println(ts, line)
                println(line)

                for _ in 2:NC
                    NUM += 1
                    TIME += DT
                    AVETIM += DT

                    timed_stage!(timings, "Later predictor: fs1!") do
                        fs1!(S)
                    end
                    ncount = timed_stage!(timings, "Pressure solve: fs2!") do
                        fs2!(S)
                    end
                    push!(ncounts, ncount)
                    div, cfl, usum, cl, cd =
                        timed_stage!(timings, "Step UBC + periodic field output") do
                            ubc!(S, NUM)
                        end
                    if NUM % NG == 0
                        line = ts_line(NUM, TIME, div, cfl, usum, cl, cd, ncount)
                        println(ts, line)
                        println(line)
                    end
                end
            end

            timed_stage!(timings, "Restart output: seedout") do
                seedout(S, NUM, TIME)
            end
            timed_stage!(timings, "Summary output: output") do
                output(S, NUM, TIME)
            end

            elapsed = UInt64(time_ns() - run_start)
            measured_stages = sum(values(timings); init=UInt64(0))
            return (timings=timings, elapsed=elapsed,
                    other=elapsed - measured_stages, ncounts=ncounts)
        end
    end
end

function state_benchmark(f, template)
    return @benchmarkable $f(S) setup=(S=deepcopy($template)) evals=1
end

function io_benchmark(f, template)
    return @benchmarkable $f(S, dir) setup=(S=deepcopy($template); dir=mktempdir(); mkpath(joinpath(dir, "TEXTDATA"))) teardown=(rm(dir; recursive=true, force=true)) evals=1
end

bench_read_setting(path) = read_setting(path)
bench_alloc_sim(p) = allocate_sim(p)
bench_calbeta!(S, dir) = cd(dir) do
    quiet() do
        calbeta!(S)
    end
end
bench_ubc_initial!(S) = ubc!(S, 1)
bench_ubc_step!(S) = ubc!(S, 1)
bench_ubc_with_output!(S, dir) = cd(dir) do
    ubc!(S, 5)
end
bench_write_puv(S, dir, num) = cd(dir) do
    write_puv(S, num)
end
bench_output(S, dir, num, time) = cd(dir) do
    quiet() do
        output(S, num, time)
    end
end
bench_seedout(S, dir, num, time) = cd(dir) do
    seedout(S, num, time)
end
bench_seedin!(S, dir) = cd(dir) do
    quiet() do
        seedin!(S)
    end
end
bench_full_run(dir) = benchmark_main(dir)

function later_step!(S)
    fs1!(S)
    fs2!(S)
    return ubc!(S, 1)
end
