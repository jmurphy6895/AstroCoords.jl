export ModEq
"""
    ModEq{T} <: AstroCoord

Modified Equinoctial Orbital Elements. 6D parameterziation of the orbit.
p -
f -
g -
h -
k -
L -

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
    @inline ModEq{T}(X::ModEq) where {T} = new{T}(X.p, X.f, X.g, X.h, X.k, X.L)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
ModEq(X::AbstractArray{T,1}) where {T} = ModEq{T}(X...)
function ModEq(p::P, f::F, g::G, h::H, k::K, L::LT) where {P,F,G,H,K,LT}
    return USM7{promote_type(P, F, G, H, K, LT)}(p, f, g, h, k, L)
end
(::Type{M})(g::StaticVector) where {M<:ModEq} = ModEq(g[1], g[2], g[3], g[4], g[5], g[6])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(M::ModEq) = SVector{6}(M.p, M.f, M.g, M.h, M.k, M.L)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
Base.one(::Type{M}) where {M<:ModEq} = M(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::ModEq, i::Int)
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
        throw(BoundsError(r, i))
    end
end
