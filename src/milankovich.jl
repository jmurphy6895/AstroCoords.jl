export Milankovich
"""
    Milankovich{T} <: AstroCoord

Milankovich Orbital Elements. 7D parameterziation of the orbit
hx - x-component of Angular Momentum Vector
hy - y-component of Angular Momentum Vector
hz - z-component of Angular Momentum Vector
ex - x-component of Eccentricity Vector
ey - y-component of Eccentricity Vector
ez - z-component of Eccentricity Vector
L - 

Constructors
Milankovich(hx, hy, hz, ex, ey, ez, L)
Milankovich(X::AbstractArray)
Milankovich(X::AstroCoord, μ::Number)

"""
struct Milankovich{T} <: AstroCoord{7, T}
    hx::T
    hy::T
    hz::T
    ex::T
    ey::T
    ez::T
    L::T
    @inline Milankovich{T}(hx, hy, hz, ex, ey, ez, L) where T = new{T}(hx, hy, hz, ex, ey, ez, L)
    @inline Milankovich{T}(X::Milankovich) where T = new{T}(X.hx, X.hy, X.hz, X.ex, X.ey, X.ez, X.L)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
Milankovich(X::AbstractArray{T, 1}) where T = Milankovich{T}(X...)
Milankovich(hx::HX, hy::HY, hz::HZ, ex::EX, ey::EY, ez::EZ, L::LT) where{HX,HY,HZ,EX,EY,EZ,LT} = Milankovich{promote_type(HX,HY,HZ,EX,EY,EZ,LT)}(hx, hy, hz, ex, ey, ez, L)
(::Type{M})(g::StaticVector) where M<:Milankovich = Milankovich(g[1], g[2], g[3], g[4], g[5], g[6], g[7])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(M::Milankovich) = SVector{7}(M.hx, M.hy, M.hz, M.ex, M.ey, M.ez, M.L)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
Base.one(::Type{M}) where M <: Milankovich = M(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::Milankovich, i::Int)
    if i == 1 return p.hx
    elseif i == 2 return p.hy
    elseif i == 3 return p.hz
    elseif i == 4 return p.ex
    elseif i == 5 return p.ey
    elseif i == 6 return p.ez
    elseif i == 7 return p.L
    else throw(BoundsError(r,i))
    end
end



export Milankovich_ODE
"""
Milankovich ODE System

Arguments:
-'u::AbstractArray': USM State Vector [H; e; L]
-'p::AbstractArray': Parameter Vector
-'t::AbstractFloat': Time
-'accel::Function': Function computing perturbing acceleration in the inertial frame 

Returns:
-'du::AbstractArray{AbstractFloat, 1}': USM ODE Change in State Vector [dH; de; dL]
"""
function Milankovich_ODE(u::AbstractVector{<:Number}, p::AbstractVector{<:Number}, t::Number, accel::Function)
    Hx, Hy, Hz, ex, ey, ez, L = u

    x, y, z, vx, vy, vz = Mil2cart(u, p)

    r = [x; y; z]
    v = [vx; vy; vz]

    μ = p[1]

    ad = accel(u, p, t)

    H = [Hx; Hy; Hz]
    ẑ = [0.; 0.; 1.]

    du[1:3] = skew_sym(r) * ad
    du[4:6] = 1/μ * (skew_sym(v) * skew_sym(r) - skew_sym(H)) * ad
    du[7] = (dot(ẑ, r)/(norm(H)*(norm(H) + dot(ẑ,H))))*dot(H,ad) + norm(H)/(norm(r)^2)
end