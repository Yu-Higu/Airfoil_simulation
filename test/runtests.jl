using Test
using AirfoilSimulation

@testset "Package load smoke test" begin
    @test AirfoilSimulation isa Module
end

include(joinpath(@__DIR__, "..", "project", "01_Porting_code_from_fortran",
                 "test_3ME_experment_improved.jl"))
