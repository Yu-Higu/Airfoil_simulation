using Test
include(joinpath(@__DIR__, "..", "src", "AirfoilCFD.jl"))
using .AirfoilCFD

d = CFDData()
@test axes(d.U, 1) == 0:385
@test axes(d.U, 2) == 0:145
@test axes(d.BETA, 1) == 0:385
@test axes(d.BETA, 2) == 0:145
@test axes(d.X, 1) == -3:390
@test axes(d.Y, 1) == -3:150
@test axes(d.P, 1) == 0:385
@test axes(d.P, 2) == 0:145

println("CFDData axis tests passed.")
