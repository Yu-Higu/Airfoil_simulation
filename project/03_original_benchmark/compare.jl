const BENCHMARK_DIR = @__DIR__
const BENCHMARK_SCRIPT = joinpath(BENCHMARK_DIR, "benchmark.jl")
const DEFAULT_CANDIDATE =
    joinpath(BENCHMARK_DIR, "src", "3ME_experiment_improved.jl")
const DEFAULT_BASELINE =
    joinpath(BENCHMARK_DIR, "..", "01_Porting_code_from_fortran",
        "3ME_experiment_improved.jl")
const DEFAULT_BASELINE_REPORT =
    joinpath(BENCHMARK_DIR, "results", "baseline", "benchmark.txt")
const NS_PER_UNIT = Dict("ns" => 1.0, "μs" => 1e3, "ms" => 1e6, "s" => 1e9)

function parse_source_args(args)
    candidate = get(ENV, "BENCHMARK_SOURCE", DEFAULT_CANDIDATE)
    baseline = get(ENV, "BENCHMARK_BASELINE", DEFAULT_BASELINE)
    name = nothing
    record_baseline = false
    replace_baseline = false
    index = 1
    while index <= length(args)
        arg = args[index]
        if arg in ("--source", "--baseline", "--baseline-source", "--name")
            index < length(args) ||
                throw(ArgumentError("$arg requires a value"))
            value = args[index+1]
            if arg == "--source"
                candidate = value
            elseif arg == "--name"
                name = value
            else
                baseline = value
            end
            index += 2
        elseif startswith(arg, "--source=")
            candidate = split(arg, "="; limit=2)[2]
            index += 1
        elseif startswith(arg, "--baseline=") ||
            startswith(arg, "--baseline-source=")
            baseline = split(arg, "="; limit=2)[2]
            index += 1
        elseif startswith(arg, "--name=")
            name = split(arg, "="; limit=2)[2]
            index += 1
        elseif arg == "--record-baseline"
            record_baseline = true
            index += 1
        elseif arg == "--replace-baseline"
            record_baseline = true
            replace_baseline = true
            index += 1
        else
            throw(ArgumentError("unknown benchmark option: $arg"))
        end
    end
    if name !== nothing
        occursin(r"^[A-Za-z0-9][A-Za-z0-9._-]*$", name) ||
            throw(ArgumentError("--name must contain only letters, numbers, '.', '_' or '-'"))
        name == "baseline" &&
            throw(ArgumentError("'baseline' is reserved for the reference report"))
    end
    record_baseline && name !== nothing &&
        throw(ArgumentError("--name cannot be used when recording the baseline"))
    return (candidate=abspath(candidate), baseline=abspath(baseline),
        name=name, record_baseline=record_baseline,
        replace_baseline=replace_baseline)
end

function run_child(source, report_path)
    cmd = `$(Base.julia_cmd()) --project=$(BENCHMARK_DIR) $(BENCHMARK_SCRIPT)`
    child_env = copy(ENV)
    child_env["BENCHMARK_CHILD"] = "1"
    child_env["BENCHMARK_SOURCE"] = source
    child_env["BENCHMARK_RESULT_PATH"] = report_path
    output = read(setenv(cmd, child_env), String)
    report_match = match(r"Saved results to: (.+)", output)
    report_match === nothing &&
        error("benchmark process did not report its results file for $source")
    println("Saved benchmark report: ", strip(report_match.captures[1]))
    saved_path = abspath(strip(report_match.captures[1]))
    isfile(saved_path) || error("benchmark report was not created: $saved_path")
    return saved_path
end

function add_metric!(metrics, category, name, value)
    metrics[(category, name)] = value
end

function read_metrics(path)
    metrics = Dict{Tuple{String,String},Float64}()
    lines = readlines(path)
    active_trial = nothing
    pending_name = ""
    for line in lines
        if startswith(line, "METRIC\t")
            fields = split(line, '\t'; limit=4)
            length(fields) == 4 && add_metric!(metrics, fields[2], fields[3],
                parse(Float64, fields[4]))
            continue
        end
        if startswith(line, "BenchmarkTools.Trial:")
            active_trial = pending_name
            pending_name = ""
        elseif active_trial !== nothing
            median_match = match(r"Time\s+\(median\):\s*([0-9.]+)\s*(ns|μs|ms|s)",
                line)
            single_match = match(r"Single result which took ([0-9.]+) (ns|μs|ms|s)",
                line)
            time_match = median_match === nothing ? single_match : median_match
            if time_match !== nothing
                add_metric!(metrics, "benchmark", active_trial,
                    parse(Float64, time_match.captures[1]) *
                        NS_PER_UNIT[time_match.captures[2]])
                active_trial = nothing
            end
        end
        if !isempty(strip(line))
            pending_name = strip(line)
        end
        sequential_match = match(r"Complete run:\s*([0-9.]+)\s*s", line)
        if sequential_match !== nothing
            add_metric!(metrics, "sequential", "Complete run",
                parse(Float64, sequential_match.captures[1]) * 1e9)
        end
    end
    return metrics
end

function report_metadata(path)
    lines = readlines(path)
    julia = findfirst(line -> startswith(line, "Julia: "), lines)
    conditions = findfirst(line -> startswith(line, "Conditions: "), lines)
    return (julia=julia === nothing ? "not recorded" : lines[julia],
        conditions=conditions === nothing ? "not recorded" : lines[conditions])
end

function old_benchmark_reports(result_dir, baseline_report)
    paths = String[]
    if isdir(result_dir)
        for (dir, _, files) in walkdir(result_dir)
            for file in files
                is_report = file == "benchmark.txt" ||
                    (startswith(file, "benchmark-") && endswith(file, ".txt"))
                if is_report
                    path = abspath(joinpath(dir, file))
                    path == baseline_report || push!(paths, path)
                end
            end
        end
    end
    sort!(paths; by=path -> (stat(path).mtime, path))
    return paths
end

function metric_label(metric)
    category, name = metric
    return category == "benchmark" ? name : "$category / $name"
end

function print_comparison(io, title, source_report, candidate_metrics)
    metrics = read_metrics(source_report)
    println(io, "\n=== ", title, " ===")
    println(io, "Report: ", source_report)
    metadata = report_metadata(source_report)
    println(io, metadata.julia)
    println(io, metadata.conditions)
    println(io, rpad("Measurement", 54), lpad("Candidate", 14),
        lpad("Compared", 14), lpad("Change", 12))

    shared = sort(collect(intersect(keys(candidate_metrics), keys(metrics)));
        by=metric_label)
    if isempty(shared)
        println(io, "No shared timing metrics are available.")
        return
    end
    for key in shared
        candidate = candidate_metrics[key]
        compared = metrics[key]
        delta = 100 * (candidate / compared - 1)
        println(io, rpad(metric_label(key), 54),
            lpad(string(round(candidate / 1e6; digits=3), " ms"), 14),
            lpad(string(round(compared / 1e6; digits=3), " ms"), 14),
            lpad(string(round(delta; digits=1), "%"), 12))
    end
end

function run_comparison()
    paths = parse_source_args(ARGS)
    result_dir = joinpath(BENCHMARK_DIR, "results")
    baseline_report = abspath(DEFAULT_BASELINE_REPORT)
    mkpath(result_dir)

    if paths.record_baseline
        isfile(paths.baseline) ||
            throw(ArgumentError("baseline source file does not exist: $(paths.baseline)"))
        isfile(baseline_report) && !paths.replace_baseline &&
            throw(ArgumentError("baseline already exists at $baseline_report; use --replace-baseline to refresh it"))
        println("Recording baseline from: ", paths.baseline)
        mkpath(dirname(baseline_report))
        temporary_dir = mktempdir(result_dir)
        try
            temporary_report = joinpath(temporary_dir, "benchmark.txt")
            saved_path = run_child(paths.baseline, temporary_report)
            mv(saved_path, baseline_report; force=paths.replace_baseline)
        finally
            rm(temporary_dir; recursive=true, force=true)
        end
        println("Baseline saved to: ", baseline_report)
        return
    end

    isfile(paths.candidate) ||
        throw(ArgumentError("candidate source file does not exist: $(paths.candidate)"))
    isfile(baseline_report) ||
        throw(ArgumentError("no saved baseline found at $baseline_report; run benchmark.jl --record-baseline first"))
    run_name = paths.name === nothing ?
               "run-$(Dates.format(now(), "yyyymmdd-HHMMSS"))" : paths.name
    run_dir = joinpath(result_dir, run_name)
    ispath(run_dir) &&
        throw(ArgumentError("benchmark run already exists at $run_dir; choose another --name"))
    history = old_benchmark_reports(result_dir, baseline_report)

    println("Candidate source: ", paths.candidate)
    println("Using saved baseline: ", baseline_report)
    println("Run name: ", run_name)
    println("Running candidate benchmark...")
    flush(stdout)

    temporary_dir = mktempdir(result_dir)
    try
        candidate_report = run_child(paths.candidate,
            joinpath(temporary_dir, "benchmark.txt"))
        candidate_metrics = read_metrics(candidate_report)
        mkpath(run_dir)
        final_report = joinpath(run_dir, "benchmark.txt")
        mv(candidate_report, final_report)
        comparison_path = joinpath(run_dir, "comparison.txt")
        open(comparison_path, "w") do io
            println(io, "Airfoil simulation benchmark comparison")
            println(io, "Date: ", now())
            println(io, "Run name: ", run_name)
            println(io, "Candidate source: ", paths.candidate)
            println(io, "Saved baseline report: ", baseline_report)
            println(io, "Candidate report: ", final_report)
            println(io, "\nPositive change percentages indicate the candidate is slower.")
            print_comparison(io, "Candidate versus saved baseline", baseline_report,
                candidate_metrics)
            for report in history
                print_comparison(io, "Candidate versus historical run: " *
                    relpath(report, result_dir), report,
                    candidate_metrics)
            end
        end
        println("Saved comparison report: ", comparison_path)
    finally
        isdir(temporary_dir) &&
            rm(temporary_dir; recursive=true, force=true)
    end
end