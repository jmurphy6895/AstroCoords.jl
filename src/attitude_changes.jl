export EP2MRP, MRP2EP

function EP2MRP(β::AbstractArray{T,1}) where {T<:Number}
    β0, β1, β2, β3 = β

    σ1 = β1 / (1.0 + β0)
    σ2 = β2 / (1.0 + β0)
    σ3 = β3 / (1.0 + β0)

    σ = [σ1; σ2; σ3]

    if norm(σ) > 1.0
        σ = -σ / (dot(σ, σ))
    end

    return σ
end

function MRP2EP(σ::AbstractArray{T,1}) where {T<:Number}
    β0 = (1.0 - norm(σ)^2) / (1.0 + norm(σ)^2)
    β = (2.0 * (σ)) ./ (1.0 + norm(σ)^2)

    return [β0; β]
end
