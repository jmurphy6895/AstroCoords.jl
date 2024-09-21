using AstroCoords
using Aqua
using Test

@testset "AstroCoords.jl" begin
    include("test_coordinate_changes.jl")
end

@testset "Aqua.jl" begin
    Aqua.test_all(AstroForceModels; ambiguities=(recursive = false))
end
