#!/usr/bin/env julia
include(joinpath(@__DIR__, "src", "AirfoilCFD.jl"))
using .AirfoilCFD

function read_number(prompt, T=Float64)
    print(prompt)
    flush(stdout)
    line = readline()
    return parse(T, strip(line))
end

re = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : read_number("REYNOLDS NUMBER, RE = ")
deg = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : read_number("ANGLE OF ATTACK, Degree = ")
nc = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 10000
ng = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 5

root = @__DIR__
airfoil = joinpath(root, "data", "clark_y.dat")
run_simulation(re, deg; nc=nc, ng=ng, airfoil_file=airfoil,
               output_dir=joinpath(root, "TEXTDATA"))
