module AstroCoordsZygoteExt

using AstroCoords

if isdefined(Base, :get_extension)
    using Zygote
    using ChainRulesCore: ChainRulesCore
    import ChainRulesCore: Tangent, NoTangent, ProjectTo
else
    using ..Zygote
    import ..Zygote.ChainRulesCore
    import ..Zygote.ChainRulesCore: Tangent, NoTangent, ProjectTo
end

function ChainRulesCore.rrule(
    new_coord::Type{<:AstroCoords.AstroCoord}, coord::AbstractArray
)
    # Forward evaluation (Keplerian transformation)
    y = new_coord(coord)

    function AstroCoord_pullback(Δ::AbstractVector)
        # Define the pullback (how gradients propagate backwards)
        Δcoords = typeof(coord)(Δ)
        return (NoTangent(), Δcoords)
    end

    function AstroCoord_pullback(Δ::Tangent)
        return (NoTangent(), collect(values(Δ)))
    end

    return y, AstroCoord_pullback
end

end
