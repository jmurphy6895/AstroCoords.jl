export Cartesian
"""
    Cartesian{T} <: AstroCoord

Cartesian Orbital Elements. 6D parameterziation of the orbit.
x - X-position
y - Y-position
z - Z-position
ẋ - X-velocity
ẏ - Y-velocity
ż - Z-velocity

Constructors
Cartesian(x, y, z, ẋ, ẏ, ż)
Cartesian(X::AbstractArray)
Cartesian(X::AstroCoord, μ::Number)

"""
struct Cartesian{T} <: AstroCoord{6,T}
    x::T
    y::T
    z::T
    ẋ::T
    ẏ::T
    ż::T
    @inline Cartesian{T}(x, y, z, ẋ, ẏ, ż) where {T} = new{T}(x, y, z, ẋ, ẏ, ż)
    @inline Cartesian{T}(p::Cartesian) where {T} = new{T}(p.x, p.y, p.z, p.ẋ, p.ẏ, p.ż)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
Cartesian(X::AbstractArray{T,1}) where {T} = Cartesian{T}(X...)
function Cartesian(x::X, y::Y, z::Z, ẋ::XV, ẏ::YV, ż::ZV) where {X,Y,Z,XV,YV,ZV}
    return Cartesian{promote_type(X, Y, Z, XV, YV, ZV)}(x, y, z, ẋ, ẏ, ż)
end
(::Type{C})(g::StaticVector) where {C<:Cartesian} = C(g[1], g[2], g[3], g[4], g[5], g[6])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::Cartesian) = SVector{6}(g.x, g.y, g.z, g.ẋ, g.ẏ, g.ż)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
Base.one(::Type{C}) where {C<:Cartesian} = C(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::Cartesian, i::Int)
    if i == 1
        return p.x
    elseif i == 2
        return p.y
    elseif i == 3
        return p.z
    elseif i == 4
        return p.ẋ
    elseif i == 5
        return p.ẏ
    elseif i == 6
        return p.ż
    else
        throw(BoundsError(r, i))
    end
end
