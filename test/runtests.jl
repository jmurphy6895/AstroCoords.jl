using AstroCoords
using Aqua
using Test

using DifferentiationInterface
using ForwardDiff: ForwardDiff
using Enzyme: Enzyme
using FiniteDifferences: FiniteDifferences
using PolyesterForwardDiff: PolyesterForwardDiff
using Zygote: Zygote
#TODO: GET THESE WORKING
#import FastDifferentiation, Mooncake, ReverseDiff, Symbolics, Tracker

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
    include("differentiability.jl")
end

@testset "Aqua.jl" begin
    Aqua.test_all(AstroCoords; ambiguities=(recursive = false))
end
