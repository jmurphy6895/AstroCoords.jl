using AstroCoords
using Test, JET

@testset "JET Testing" begin
    rep = JET.test_package(AstroCoords; toplevel_logger=nothing)
end
