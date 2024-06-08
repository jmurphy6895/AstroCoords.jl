export Keplerian
"""
    Keplerian{T} <: AstroCoord

Keplerian Orbital Elements. 6D parameterziation of the orbit
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
struct Keplerian{T} <: AstroCoord{6, T}
    a::T
    e::T
    i::T
    Ω::T
    ω::T
    f::T
    @inline Keplerian{T}(a, e, i, Ω, ω, f) where T = new{T}(a, e, i, Ω, ω, f)
    @inline Keplerian{T}(p::Keplerian) where T = new{T}(p.a, p.e, p.i, p.Ω, p.ω, p.f)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
Keplerian(X::AbstractArray{T, 1}) where T = Keplerian{T}(X...)
Keplerian(a::A, e::E, i::I, Ω::O, ω::W, f::F) where{A,E,I,O,W,F} = Keplerian{promote_type(A,E,I,O,W,F)}(a, e, i, Ω, ω, f)
(::Type{K})(g::StaticVector) where K<:Keplerian = K(g[1], g[2], g[3], g[4], g[5], g[6])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::Keplerian) = SVector{6}(g.a, g.e, g.i, g.Ω, g.ω, g.f)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
Base.one(::Type{K}) where K <: Keplerian = K(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# ~~~~~~~~~~~~~~~~ Coordinate Changes ~~~~~~~~~~~~~~~~~~ #
@inline function (::Type{C})(p::Keplerian, μ::T) where C <: Cartesian where T

    return C(koe2cart(p, μ))

end

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::Keplerian, i::Int)
    if i == 1 return p.a
    elseif i == 2 return p.e
    elseif i == 3 return p.i
    elseif i == 4 return p.Ω
    elseif i == 5 return p.ω
    elseif i == 6 return p.f
    else throw(BoundsError(r,i))
    end
end

export GaussVE
"""
Gaussian Variational Equations for the Keplerian Orbital Elements

Arguments:
-'u::AbstractArray': Keplerian State Vector [a; e; i; Ω; ω; f]
-'p::AbstractArray': Parameter Vector
-'t::AbstractFloat': Time
-'accel::Function': Function computing perturbing acceleration in the inertial frame 

Returns:
-'du::AbstractArray{AbstractFloat, 1}': Keplerian ODE Change in State Vector [da; de; di; dΩ; dω; df]
"""
function GaussVE(u::AbstractVector{<:Number}, p::AbstractVector{<:Number}, t::Number, accel::Function)

    a, e, i, Ω, ω, f = u
    μ = p[1]

    acc = accel(u, p, t)

    p = a * (1.0 - e^2)
    r = p/(1.0 + e*cos(f))
    θ = ω + f
    h = √(μ*p)

    ar = acc[1]
    aθ = acc[2] * sin(f)
    ah = acc[3]

    da = (2.0 * a^2)/h * (e*sin(f)*ar + p/r*aθ)
    de = 1.0/h * (p*sin(f)*ar + ((p+r)*cos(f) + r*e)*aθ)
    di = (r*cos(θ)/h) * ah
    dΩ = (r*sin(θ))/(h*sin(i)) * ah
    dω = 1.0/(h*e) * (-p*cos(f)*ar + (p+r)*sin(f)*aθ) - (r*sin(θ)*cos(i))/(h*sin(i))*ah
    df = h/(r^2) + 1.0/(h*e) * (p*cos(f)*ar - (p+r)*sin(f)*aθ)

    return [da; de; di; dΩ; dω; df]

end

#TODO: Change To Same Class as Keplerian?
export KeplerianM
"""
    KeplerianM{T} <: AstroCoord

Keplerian Orbital Elements. 6D parameterziation of the orbit
a - semi-major axis
e - eccetricity 
i - inclination 
Ω - Right Ascension of Ascending Node 
ω - Argument of Perigee
f - True Anomaly

Constructors
Keplerian(a, e, i, Ω, ω, f)
Keplerian(X::AbstractArray)

"""
struct KeplerianM{T} <: AstroCoord{6, T}
    a::T
    e::T
    i::T
    Ω::T
    ω::T
    M::T
    @inline KeplerianM{T}(a, e, i, Ω, ω, M) where T = new{T}(a, e, i, Ω, ω, M)
    @inline KeplerianM{T}(p::KeplerianM) where T = new{T}(p.a, p.e, p.i, p.Ω, p.ω, p.M)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
KeplerianM(X::AbstractArray{T, 1}) where T = KeplerianM{T}(X...)
KeplerianM(a::A, e::E, i::I, Ω::O, ω::W, M::MT) where{A,E,I,O,W,MT} = Keplerian{promote_type(A,E,I,O,W,MT)}(a, e, i, Ω, ω, M)
(::Type{K})(g::StaticVector) where K<:KeplerianM = K(g[1], g[2], g[3], g[4], g[5], g[6])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::KeplerianM) = SVector{6}(g.a, g.e, g.i, g.Ω, g.ω, g.f)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
Base.one(::Type{K}) where K <: KeplerianM = K(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# ~~~~~~~~~~~~~~~~ Coordinate Changes ~~~~~~~~~~~~~~~~~~ #
@inline function (::Type{C})(p::KeplerianM, μ::T) where C <: Cartesian where T

    return C(koeM2cart(p, μ))

end

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::KeplerianM, i::Int)
    if i == 1 return p.a
    elseif i == 2 return p.e
    elseif i == 3 return p.i
    elseif i == 4 return p.Ω
    elseif i == 5 return p.ω
    elseif i == 6 return p.M
    else throw(BoundsError(r,i))
    end
end


#TODO: Gauss Variational for M
export GaussVE
"""
Gaussian Variational Equations for the Keplerian Orbital Elements

Arguments:
-'u::AbstractArray': Keplerian State Vector [a; e; i; Ω; ω; M]
-'p::AbstractArray': Parameter Vector
-'t::AbstractFloat': Time
-'accel::Function': Function computing perturbing acceleration in the inertial frame 

Returns:
-'du::AbstractArray{AbstractFloat, 1}': Keplerian ODE Change in State Vector [da; de; di; dΩ; dω; dM]
"""
function GaussMVE(u::AbstractVector{<:Number}, p::AbstractVector{<:Number}, t::Number, accel::Function)

    a, e, i, Ω, ω, f = u
    μ = p[1]

    acc = accel(u, p, t)

    p = a * (1.0 - e^2)
    r = p/(1.0 + e*cos(f))
    θ = ω + f
    h = √(μ*p)

    ar = acc[1]
    aθ = acc[2] * sin(f)
    ah = acc[3]

    da = (2.0 * a^2)/h * (e*sin(f)*ar + p/r*aθ)
    de = 1.0/h * (p*sin(f)*ar + ((p+r)*cos(f) + r*e)*aθ)
    di = (r*cos(θ)/h) * ah
    dΩ = (r*sin(θ))/(h*sin(i)) * ah
    dω = 1.0/(h*e) * (-p*cos(f)*ar + (p+r)*sin(f)*aθ) - (r*sin(θ)*cos(i))/(h*sin(i))*ah
    dM = h/(r^2) + 1.0/(h*e) * (p*cos(f)*ar - (p+r)*sin(f)*aθ)

    return [da; de; di; dΩ; dω; df]

end

#TODO: Keplerian Propagator
#export KeplerianPropagator
"""
Gaussian Variational Equations for the Keplerian Orbital Elements

Arguments:
-'u::AbstractArray': Keplerian State Vector [a; e; i; Ω; ω; M]
-'p::AbstractArray': Parameter Vector
-'t::AbstractFloat': Time
-'accel::Function': Function computing perturbing acceleration in the inertial frame 

Returns:
-'du::AbstractArray{AbstractFloat, 1}': Keplerian ODE Change in State Vector [da; de; di; dΩ; dω; dM]
"""