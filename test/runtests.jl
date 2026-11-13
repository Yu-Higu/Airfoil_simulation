using Test
using AirfoilSimulation

@testset "Package load smoke test" begin
    @test AirfoilSimulation isa Module
end
