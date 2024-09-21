export Milankovich
"""
    Milankovich{T} <: AstroCoord

Milankovich Orbital Elements. 7D parameterziation of the orbit.
hx - X-component of Angular Momentum Vector
hy - Y-component of Angular Momentum Vector
hz - Z-component of Angular Momentum Vector
ex - X-component of Eccentricity Vector
ey - Y-component of Eccentricity Vector
ez - Z-component of Eccentricity Vector
L - Mean Anomaly

Constructors
Milankovich(hx, hy, hz, ex, ey, ez, L)
Milankovich(X::AbstractArray)
Milankovich(X::AstroCoord, μ::Number)

"""
struct Milankovich{T} <: AstroCoord{7,T}
    hx::T
    hy::T
    hz::T
    ex::T
    ey::T
    ez::T
    L::T
    @inline Milankovich{T}(hx, hy, hz, ex, ey, ez, L) where {T} =
        new{T}(hx, hy, hz, ex, ey, ez, L)
    @inline Milankovich{T}(X::Milankovich) where {T} =
        new{T}(X.hx, X.hy, X.hz, X.ex, X.ey, X.ez, X.L)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
Milankovich(X::AbstractArray{T,1}) where {T} = Milankovich{T}(X...)
function Milankovich(
    hx::HX, hy::HY, hz::HZ, ex::EX, ey::EY, ez::EZ, L::LT
) where {HX,HY,HZ,EX,EY,EZ,LT}
    return Milankovich{promote_type(HX, HY, HZ, EX, EY, EZ, LT)}(hx, hy, hz, ex, ey, ez, L)
end
function (::Type{M})(g::StaticVector) where {M<:Milankovich}
    return Milankovich(g[1], g[2], g[3], g[4], g[5], g[6], g[7])
end

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(M::Milankovich) = SVector{7}(M.hx, M.hy, M.hz, M.ex, M.ey, M.ez, M.L)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
Base.one(::Type{M}) where {M<:Milankovich} = M(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::Milankovich, i::Int)
    if i == 1
        return p.hx
    elseif i == 2
        return p.hy
    elseif i == 3
        return p.hz
    elseif i == 4
        return p.ex
    elseif i == 5
        return p.ey
    elseif i == 6
        return p.ez
    elseif i == 7
        return p.L
    else
        throw(BoundsError(r, i))
    end
end
