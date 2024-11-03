@testset "Test Constructor" begin
    u0 = SVector{6}(randn(6)...)

    for set in _COORDINATE_SETS
        if set ∈ [Milankovich, USM7]
            u0 = SVector{7}(randn(7)...)
        else
            u0 = SVector{6}(randn(6)...)
        end

        coord = set(u0)
        @test params(coord) == u0
    end
end