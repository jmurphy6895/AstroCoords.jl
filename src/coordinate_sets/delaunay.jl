export Delaunay
"""
    Delaunay{T} <: AstroCoord

Delaunay Orbital Elements. 6D parameterziation of the orbit.
L - Canonical Keplerian Energy
G - Canonical Total Angular Momentum
H - Canonical Normal Angular Momentum (Relative to Equator)
M - Mean Anomaly
ω - Arugment of Periapsis
Ω - Right Ascention of the Ascending Node

Constructors
Delaunay(L, G, H, M, ω, Ω)
Delaunay(X::AbstractVector{<:Number})
Delaunay(X::AstroCoord, μ::Number)

"""
struct Delaunay{T} <: AstroCoord{6,T}
    L::T
    G::T
    H::T
    M::T
    ω::T
    Ω::T
    @inline Delaunay{T}(L, G, H, M, ω, Ω) where {T} = new{T}(L, G, H, M, ω, Ω)
    @inline Delaunay{T}(p::Delaunay) where {T} = new{T}(p.L, p.G, p.H, p.M, p.ω, p.Ω)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
Delaunay(X::AbstractVector{T}) where {T} = Delaunay{T}(X...)
function Delaunay(L::LT, G::GT, H::HT, M::MT, ω::PT, Ω::OmT) where {LT,GT,HT,MT,PT,OmT}
    return Delaunay{promote_type(LT, GT, HT, MT, PT, Omt)}(L, G, H, M, ω, Ω)
end
(::Type{D})(g::StaticVector) where {D<:Delaunay} = D(g[1], g[2], g[3], g[4], g[5], g[6])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::Delaunay) = SVector{6}(g.L, g.G, g.H, g.M, g.ω, g.Ω)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
Base.one(::Type{D}) where {D<:Delaunay} = D(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::Delaunay, i::Int)
    if i == 1
        return p.L
    elseif i == 2
        return p.G
    elseif i == 3
        return p.H
    elseif i == 4
        return p.M
    elseif i == 5
        return p.ω
    elseif i == 6
        return p.Ω
    else
        throw(BoundsError(r, i))
    end
end
