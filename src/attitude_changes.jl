export EP2MRP, MRP2EP

function EP2MRP(β::AbstractVector{T}) where {T<:Number}
    β0, β1, β2, β3 = β

    σ1 = β1 / (1.0 + β0)
    σ2 = β2 / (1.0 + β0)
    σ3 = β3 / (1.0 + β0)

    σ = SVector{3}(σ1, σ2, σ3)

    σ_sq = sum(abs2.(σ))

    if √(σ_sq) > 1.0
        σ = -σ / (σ_sq)
    end

    return σ
end

function MRP2EP(σ::AbstractVector{T}) where {T<:Number}
    σ_sq = sum(abs2.(σ))

    β0 = (1.0 - σ_sq) / (1.0 + σ_sq)
    β = (2.0 * (σ)) ./ (1.0 + σ_sq)

    return SVector{4}(β0, β[1], β[2], β[3])
end
