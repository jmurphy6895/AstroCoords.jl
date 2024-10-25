########################################################################################
#
# Tests for Differentiation Across the Coordinate Sets
#
########################################################################################
# Currently Supported:
# Enzyme, ForwardDiff, FiniteDiff, FiniteDifferences, PolyesterForwardDiff
#
# Partially Supported:
# Zygote
#
# Not Yet Supported:
# FastDifferentiation, Mooncake, ReverseDiff, Symbolics, Tracker
########################################################################################
using AstroCoords
using Test

const _COORDINATE_SETS2 = [
    #Cartesian, #TODO: Skipped for DifferentiationInterface -- Zygote
    Delaunay,
    Keplerian,
    Milankovich,
    ModEq,
    #Cylindrical, #TODO: Skipped for DifferentiationInterface -- Zygote
    #Spherical, #TODO: Skipped for DifferentiationInterface -- Zygote
    USM7,
    USM6,
    USMEM,
]

const _BACKENDS = (
    ("ForwardDiff", AutoForwardDiff()),
    #("Diffractor", AutoDiffractor()), #! Error with iszero()?
    ("Enzyme", AutoEnzyme()),
    #("FastDifferentiation", AutoFastDifferentiation()), #! Doesn't Yet Support if Statements
    ("FiniteDifferences", AutoFiniteDifferences(; fdm=FiniteDifferences.central_fdm(5, 1))),
    #("Mooncake", AutoMooncake(;config=nothing)), #! Problem with StaticVector
    ("PolyesterForwardDiff", AutoPolyesterForwardDiff()),
    #("ReverseDiff", AutoReverseDiff()), #! Problem with normalize()
    #("Symbolics", AutoSymbolics()), #! Problem with normalize()
    #("Tracker", AutoTracker()), #! Problem with AstroCoord Construction, tries to convert TrackedReal to Float
    ("Zygote", AutoZygote()), #TODO: DiffInterface Zygote Doesn't Handle Constant Functions Yet (Stand-alone Zygote should work)
)

@testset "Coordinate Transformation Differentiation" begin
    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ]

    μ = 3.986004415e5

    for set in _COORDINATE_SETS
        f_fd, df_fd = value_and_jacobian(
            (x) -> set(Cartesian(x), μ), AutoFiniteDiff(), state
        )

        f_fd2, df_fd2 = value_and_derivative(
            (x) -> set(Cartesian(state), x), AutoFiniteDiff(), μ
        )

        for backend in _BACKENDS3
            @eval @testset $("Coordinate Set $set " * string(backend[1])) begin
                f_ad, df_ad = value_and_jacobian(
                    (x) -> $set(Cartesian(x), μ), $backend[2], $state
                )

                @test $f_fd == f_ad
                @test $df_fd ≈ df_ad rtol = 1e-4

                f_ad2, df_ad2 = value_and_derivative(
                    (x) -> $set(Cartesian($state), x), $backend[2], $μ
                )

                @test $f_fd2 == f_ad2
                @test $df_fd2 ≈ df_ad2 rtol = 1e-4
            end
        end
    end
end
