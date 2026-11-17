# Julia syntax adaptation

This directory contains an independent adaptation of
`../01_Porting_code_from_fortran/3ME_experiment_improved.jl`. The numerical
scheme, existing array bounds, output formats, and uppercase simulation-state
fields are preserved for compatibility. Local variables use conventional
lowercase Julia names.

The first optimization step adds persistent `ROWMAX` and `RHS` work arrays to
the simulation state. `fs2!` computes the pressure-equation right-hand side
once per step and reuses it during SOR; the convergence criterion and update
order are unchanged.

Run the regression tests from the repository root:

```sh
julia --project=. project/04_Adapting_syntax_for_Julia/test_3ME_experiment_julia.jl
```

The original implementation and the benchmark working copy remain unchanged.
