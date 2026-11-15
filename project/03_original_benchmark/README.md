# Original Airfoil simulation benchmark

The benchmarked code is the isolated working copy at
[`src/3ME_experiment_improved.jl`](./src/3ME_experiment_improved.jl). It was
initially copied byte-for-byte from
`../01_Porting_code_from_fortran/3ME_experiment_improved.jl`; the original is
left untouched. Edit the working copy for optimization experiments.

The scripts continue to use `setting.dat` and `clark_y.dat` from
`../01_Porting_code_from_fortran/`. The settings are copied unchanged into
temporary run directories, so current runs use `IG=384`, `JG=144`, `RE=500`,
`DEG=20`, `DT=0.0025`, `NC=50`, and `NG=5`. Change that `setting.dat` only
when intentionally changing benchmark conditions.

## Optimization loop

1. Edit `src/3ME_experiment_improved.jl`; do not edit the original under
   `../01_Porting_code_from_fortran/`.
2. Run the benchmark and profiler below with unchanged settings, Julia version,
   and thread count.
3. Compare the newly generated files in `results/` with the previous baseline.
4. Keep each experiment in version control (for example, make a commit before
   starting the next optimization) so it is easy to compare or revert one
   change without losing the original.

## Setup

```sh
cd project/03_original_benchmark
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Measure each stage and identify bottlenecks

```sh
julia --project=. benchmark.jl
```

`BenchmarkTools` reports median time, memory, and allocation count for:
setting-file read, array allocation, mesh/airfoil setup, initial predictor and
pressure solve, later predictor/pressure solve/boundary stages sampled at
early/middle/final steps, periodic field output, summary output, restart
input/output, and the complete `main()` run.
Mutating calculations are measured from a fresh copy of the same prepared
stage state on each sample. Output measurements use temporary directories;
directory setup and cleanup are outside the timed expression. The full run
uses the original `main()` and includes all computation and output.

The script also performs one sequential pass with the exact settings, timing
each phase on the live evolving state and reporting each phase's total and
percentage of wall time. This captures changes in solver cost across the
simulation. The report records SOR iteration counts for the first, early,
middle, and final solver states. BenchmarkTools' complete `main()` measurement
remains the uninstrumented end-to-end reference.

Set `BENCHMARK_SECONDS` to change the duration per target (default 5 seconds):

```sh
BENCHMARK_SECONDS=15 julia --project=. benchmark.jl
```

Results are saved under `results/`. Compare runs only with the same Julia
version, machine, thread count, and unchanged `setting.dat`. Run on an
otherwise idle machine; warm-up and state-copy setup are excluded by
BenchmarkTools.

## CPU, allocation, and type-stability profiles

```sh
julia --project=. profile.jl
```

The report includes `@code_warntype` output for the main numerical and I/O
functions, a CPU profile of one representative solver step, a CPU profile of
the complete `main()` workflow, and a sampled allocation profile of one solver
step. The complete-run CPU profile uses `PROFILE_REPEATS` runs (default 1).
`ALLOC_SAMPLE_RATE` controls allocation sampling (default 0.001) to keep the
large SOR allocation trace manageable:

```sh
PROFILE_REPEATS=2 ALLOC_SAMPLE_RATE=0.005 julia --project=. profile.jl
```

Use the CPU profile's hottest functions/lines to locate work, check the
allocation profile for avoidable allocations on those paths, then inspect
`@code_warntype` for unexpected `Any` or unstable values. Profiling adds
overhead; use it to locate bottlenecks, and use the unprofiled benchmark to
compare performance changes. `seedin!` is measured/analyzed as a restart
utility but is not called by the original `main()` workflow.
