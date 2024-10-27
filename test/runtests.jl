using AstroCoords
using Aqua
using Test

using DifferentiationInterface
using Enzyme
using FiniteDiff
using FiniteDifferences
using ForwardDiff
using PolyesterForwardDiff
using Zygote
#TODO: GET THESE WORKING
#import FastDifferentiation, Mooncake, ReverseDiff, Symbolics, Tracker

using AllocCheck

const _COORDINATE_SETS = [
    Cartesian,
    Delaunay,
    Keplerian,
    Milankovich,
    ModEq,
    Cylindrical,
    Spherical,
    USM7,
    USM6,
    USMEM,
]

@testset "AstroCoords.jl" begin
    include("test_coordinate_changes.jl")
    include("test_anomalies.jl")
    include("test_quantities.jl")
end

@testset "Differentiation" begin
    include("test_differentiability.jl")
end

@testset "Code Performance" begin
    include("test_JET.jl")
    include("test_allocs.jl")
end

@testset "Aqua.jl" begin
    Aqua.test_all(AstroCoords; ambiguities=(recursive = false))
end
