export Spherical
"""
    Spherical{T} <: AstroCoord

Spherical Orbital Elements. 6D parameterziation of the orbit
r - radius
θ - in-plane angle
ϕ - out-of-plane angle
ṙ - instantaneous rate of change of radius
θdot - instantaneous rate of change of θ
ϕdot - instantaneous rate of change of ϕ

Constructors
Spherical(r, θ, ϕ, ṙ, ω, Ω)
Spherical(X::AbstractArray)
Spherical(X::AstroCoord, μ::Number)

"""
struct Spherical{T} <: AstroCoord{6,T}
    r::T
    θ::T
    ϕ::T
    ṙ::T
    θdot::T
    ϕdot::T
    @inline Spherical{T}(r, θ, ϕ, ṙ, θdot, ϕdot) where {T} =
        new{T}(r, θ, ϕ, ṙ, θdot, ϕdot)
    @inline Spherical{T}(X::Spherical) where {T} =
        new{T}(X.r, X.θ, X.ϕ, X.ṙ, X.θdot, X.ϕdot)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
Spherical(X::AbstractArray{T,1}) where {T} = Spherical{T}(X...)
function Spherical(r::R, θ::T, ϕ::P, ṙ::RD, θdot::TD, ϕdot::PD) where {R,T,P,RD,TD,PD}
    return Spherical{promote_type(R, T, P, RD, TD, PD)}(r, θ, ϕ, ṙ, θdot, ϕdot)
end
function (::Type{S})(g::StaticVector) where {S<:Spherical}
    return Spherical(g[1], g[2], g[3], g[4], g[5], g[6])
end

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(S::Spherical) = SVector{6}(S.r, S.θ, S.ϕ, S.ṙ, S.θdot, S.ϕdot)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
Base.one(::Type{S}) where {S<:Spherical} = S(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::Spherical, i::Int)
    if i == 1
        return p.r
    elseif i == 2
        return p.θ
    elseif i == 3
        return p.ϕ
    elseif i == 4
        return p.ṙ
    elseif i == 5
        return p.θdot
    elseif i == 6
        return p.ϕdot
    else
        throw(BoundsError(r, i))
    end
end

#TODO: Spherical ODE

export Cylindrical
"""
    Cylindrical{T} <: AstroCoord

Cylindrical Orbital Elements. 6D parameterziation of the orbit
r - radius
θ - in-plane angle
z - out-of-plane distance
ṙ - instantaneous rate of change of radius
θdot - instantaneous rate of change of θ
ż - instantaneous rate of change of z

Constructors
Cylindrical(r, θ, z, ṙ, θdot, ż)
Cylindrical(X::AbstractArray)
Cylindrical(X::AstroCoord, μ::Number)

"""
struct Cylindrical{T} <: AstroCoord{6,T}
    r::T
    θ::T
    z::T
    ṙ::T
    θdot::T
    ż::T
    @inline Cylindrical{T}(r, θ, z, ṙ, θdot, ż) where {T} = new{T}(r, θ, z, ṙ, θdot, ż)
    @inline Cylindrical{T}(X::Cylindrical) where {T} =
        new{T}(X.r, X.θ, X.z, X.ṙ, X.θdot, X.ż)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
Cylindrical(X::AbstractArray{T,1}) where {T} = Cylindrical{T}(X...)
function Cylindrical(r::R, θ::T, z::Z, ṙ::RD, θdot::TD, ż::ZD) where {R,T,Z,RD,TD,ZD}
    return Cylindrical{promote_type(R, T, P, RD, TD, ZD)}(r, θ, z, ṙ, θdot, ż)
end
function (::Type{C})(g::StaticVector) where {C<:Cylindrical}
    return Cylindrical(g[1], g[2], g[3], g[4], g[5], g[6])
end

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(C::Cylindrical) = SVector{6}(C.r, C.θ, C.z, C.ṙ, C.θdot, C.ż)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
Base.one(::Type{C}) where {C<:Cylindrical} = C(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::Cylindrical, i::Int)
    if i == 1
        return p.r
    elseif i == 2
        return p.θ
    elseif i == 3
        return p.z
    elseif i == 4
        return p.ṙ
    elseif i == 5
        return p.θdot
    elseif i == 6
        return p.ż
    else
        throw(BoundsError(r, i))
    end
end
