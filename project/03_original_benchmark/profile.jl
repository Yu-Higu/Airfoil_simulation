using BenchmarkTools
using InteractiveUtils
using Profile
using Dates
include("common.jl")

repeats = parse(Int, get(ENV, "PROFILE_REPEATS", "1"))
repeats >= 1 || throw(ArgumentError("PROFILE_REPEATS must be at least one"))
allocation_sample_rate = parse(Float64, get(ENV, "ALLOC_SAMPLE_RATE", "0.001"))
0 < allocation_sample_rate <= 1 ||
    throw(ArgumentError("ALLOC_SAMPLE_RATE must be in the interval (0, 1]"))

source_dir = make_workdir()
profile_dir = make_workdir()
try
    p = read_setting(joinpath(source_dir, "setting.dat"))
    mesh_state = prepare_mesh_state(source_dir, p)
    first_state = prepare_first_step_state(mesh_state)
    post_first_state = prepare_post_first_step_state(first_state)
    final_state = post_first_state

    result_dir = joinpath(@__DIR__, "results")
    mkpath(result_dir)
    result_path = joinpath(result_dir,
                           "profile-$(Dates.format(now(), "yyyymmdd-HHMMSS")).txt")

    open(result_path, "w") do io
        println(io, "Airfoil simulation profile using project/01 setting.dat")
        println(io, "Date: ", now())
        println(io, "Julia: ", VERSION, "; CPU: ", Sys.CPU_NAME,
                "; threads: ", Threads.nthreads())
        println(io, "Conditions: IG=", p.IG, ", JG=", p.JG, ", RE=", p.RE,
                ", DEG=", p.DEG, ", DT=", p.DT, ", NC=", p.NC, ", NG=", p.NG)
        println(io, "Complete-run repetitions: ", repeats)
        println(io, "Allocation sample rate: ", allocation_sample_rate)

        println(io, "\n=== @code_warntype ===")
        for f in (read_setting, allocate_sim, calbeta!, fs0!, fs1!, fs2!, ubc!,
                  write_puv, output, seedout, seedin!)
            println(io, "\n--- ", f, " ---")
            argtypes = f === read_setting ? (String,) :
                       f === allocate_sim ? (typeof(p),) :
                       f === write_puv ? (typeof(final_state), Int) :
                       f === output || f === seedout ? (typeof(final_state), Int, Float64) :
                       f === seedin! ? (typeof(final_state),) :
                       f === ubc! ? (typeof(final_state), Int) :
                       (typeof(final_state),)
            code_warntype(io, f, argtypes)
        end

        println(io, "\n=== One representative FS1! + FS2! + UBC! CPU profile ===")
        step_state = deepcopy(post_first_state)
        Profile.clear()
        @profile later_step!(step_state)
        Profile.print(io; format=:flat, sortedby=:count, maxdepth=25)

        println(io, "\n=== Complete main() CPU profile ===")
        benchmark_main(profile_dir)
        Profile.clear()
        @profile for _ in 1:repeats
            benchmark_main(profile_dir)
        end
        Profile.print(io; format=:flat, sortedby=:count, maxdepth=25)

        println(io, "\n=== One representative step allocation profile ===")
        Profile.Allocs.clear()
        allocation_state = deepcopy(post_first_state)
        Profile.Allocs.@profile sample_rate=allocation_sample_rate begin
            later_step!(allocation_state)
        end
        allocations = Profile.Allocs.fetch()
        Profile.Allocs.print(io, allocations; format=:flat, sortedby=:count, maxdepth=25)
        Profile.Allocs.clear()
    end

    println("Saved CPU, allocation, and type-stability reports to: ", result_path)
finally
    rm(source_dir; recursive=true, force=true)
    rm(profile_dir; recursive=true, force=true)
end
