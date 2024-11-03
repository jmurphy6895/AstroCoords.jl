export ModEq
"""
    ModEq{T} <: AstroCoord

Modified Equinoctial Orbital Elements. 6D parameterziation of the orbit.
p - semi-parameter
f - eccetricity projection onto longitude of perigee
g - eccetricity projection onto ⟂ longitude of perigee
h - projection of half inclination onto RAAN
k - projection of half inclination onto ⟂ RAAN
L - true longitude

Constructors
ModEq(p, f, g, h, k, l)
ModEq(X::AbstractArray)
ModEq(X::AstroCoord, μ::Number)

"""
struct ModEq{T} <: AstroCoord{6,T}
    p::T
    f::T
    g::T
    h::T
    k::T
    L::T
    @inline ModEq{T}(p, f, g, h, k, L) where {T} = new{T}(p, f, g, h, k, L)
    @inline ModEq{T}(X::ModEq{T}) where {T} = new{T}(X.p, X.f, X.g, X.h, X.k, X.L)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
ModEq(X::AbstractVector{T}) where {T} = ModEq{T}(X[1], X[2], X[3], X[4], X[5], X[6])
function ModEq(p::P, f::F, g::G, h::H, k::K, L::LT) where {P,F,G,H,K,LT}
    return ModEq{promote_type(P, F, G, H, K, LT)}(p, f, g, h, k, L)
end

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::ModEq{T}) where {T<:Number} = SVector{6,T}(g.p, g.f, g.g, g.h, g.k, g.L)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
function Base.one(::Type{M}; T::DataType=Float64) where {M<:ModEq}
    return M{T}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::ModEq{T}, i::Int) where {T<:Number}
    if i == 1
        return p.p
    elseif i == 2
        return p.f
    elseif i == 3
        return p.g
    elseif i == 4
        return p.h
    elseif i == 5
        return p.k
    elseif i == 6
        return p.L
    else
        throw(BoundsError(p, i))
    end
end
