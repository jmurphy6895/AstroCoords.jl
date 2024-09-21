export Keplerian
"""
    Keplerian{T} <: AstroCoord

Keplerian Orbital Elements. 6D parameterziation of the orbit.
a - semi-major axis
e - eccetricity 
i - inclination 
Ω - Right Ascension of Ascending Node 
ω - Argument of Perigee
f - True Anomaly

Constructors
Keplerian(a, e, i, Ω, ω, f)
Keplerian(X::AbstractArray)
Keplerian(X::AstroCoord, μ::Number)

"""
struct Keplerian{T} <: AstroCoord{6,T}
    a::T
    e::T
    i::T
    Ω::T
    ω::T
    f::T
    @inline Keplerian{T}(a, e, i, Ω, ω, f) where {T} = new{T}(a, e, i, Ω, ω, f)
    @inline Keplerian{T}(p::Keplerian) where {T} = new{T}(p.a, p.e, p.i, p.Ω, p.ω, p.f)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
Keplerian(X::AbstractArray{T,1}) where {T} = Keplerian{T}(X...)
function Keplerian(a::A, e::E, i::I, Ω::O, ω::W, f::F) where {A,E,I,O,W,F}
    return Keplerian{promote_type(A, E, I, O, W, F)}(a, e, i, Ω, ω, f)
end
(::Type{K})(g::StaticVector) where {K<:Keplerian} = K(g[1], g[2], g[3], g[4], g[5], g[6])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::Keplerian) = SVector{6}(g.a, g.e, g.i, g.Ω, g.ω, g.f)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
Base.one(::Type{K}) where {K<:Keplerian} = K(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::Keplerian, i::Int)
    if i == 1
        return p.a
    elseif i == 2
        return p.e
    elseif i == 3
        return p.i
    elseif i == 4
        return p.Ω
    elseif i == 5
        return p.ω
    elseif i == 6
        return p.f
    else
        throw(BoundsError(r, i))
    end
end

#TODO: Return other anomalies
