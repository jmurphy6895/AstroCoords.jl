using SALAMANDER

export cart2koe
"""
Computes the Keplerian Orbital Elements from a Cartesian Set

Arguments:
-'u::AbstractVector{<:Number}': Cartesian State Vector [x; y; z; ẋ; ẏ; ż]
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_koe::Vector{<:Number}': Keplerian Orbital Element Vector [a; e; i; Ω(RAAN); ω(AOP); f(True Anomaly)] 
*All Angles are in Radians
"""
function cart2koe(u::AbstractVector{T}, μ::Number) where {T<:Number}
    x, y, z, ẋ, ẏ, ż = u

    r = [x, y, z]
    v = [ẋ, ẏ, ż]

    rmag = norm(r)
    vmag = norm(v)

    #* Angular Momentum
    h = cross(r, v)
    hmag = norm(h)

    k̂ = [0., 0., 1.]

    #* Inclination
    i = AngleBetweenVectors.angle(h, k̂)

    #* Line of Nodes
    n = cross(k̂, h)
    Ω = atan(n[2], n[1])

    #* Semi-Major Axis
    a = 1.0/(2.0/rmag - vmag^2/μ)

    #* Eccentricity
    e = cross(v, h)/μ - r/rmag
    emag = norm(e)

    #* RAAN, AOP, True Anomaly
    if abs(i) < 1E-15 && abs(emag) < 1E-15
        Ω = 0.0
        ω = 0.0
        f = rem2pi(atan(y, x), RoundDown)
    elseif abs(i) < 1E-15 && abs(emag) > 1E-15
        Ω = 0.0
        ω = rem2pi(atan(e[2], e[1]), RoundDown)
        f = rem2pi(atan(dot(h, cross(e, r))/norm(h), dot(r, e)), RoundDown)
    elseif abs(i) > 1E-15 && abs(emag) < 1E-15
        Ω = rem2pi(atan(n[2], n[1]), RoundDown)
        ω = 0.0
        f = rem2pi(atan(dot(r, cross(h. n))/norm(h), dot(r, n)), RoundDown)
    else
        #* Eccentric & True Anomaly
        if a > 0.0
            eSE = dot(r, v) / √(μ * a)
            eCE = (rmag * vmag^2)/μ - 1.0
            emag = norm([eSE; eCE])
            E = atan(eSE, eCE)
            f = eccentricAnomaly2TrueAnomaly(E, emag)
        else
            eSH = dot(r, v) / √(-μ * a)
            eCH = (rmag * vmag^2)/μ - 1.0
            emag = √(1.0 - hmag^2/(μ*a))
            F = log((eCH + eSH)/(eCH - eSH))/2.0
            f = 2.0*atan(√(1.0+e)*sinh(F/2.0), √(e-1.0)*cosh(F/2.0))
        end

        #* AOP
        px = dot(r, n)
        py = dot(r, cross(h, n)/hmag)
        ω = atan(py, px) - f
    end

    return [a; emag; i; Ω; ω; f]

end

export koe2cart
"""
Computes the Cartesian Orbital Elements from a Keplerian Set

Arguments:
-'u::AbstractVector{<:Number}': Keplerian State Vector [a; e; i; Ω(RAAN); ω(AOP); f(True Anomaly)]  
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_cart::Vector{<:Number}': Keplerian Orbital Element Vector [x; y; z; ẋ; ẏ; ż]
*All Angles are in Radians
"""
function koe2cart(u::AbstractVector{<:Number}, μ::Number)
    a, e, i, Ω, ω, f = u

    rmag = (a * (1.0 - e^2)/(1.0 + e*cos(f)))

    θ = ω + f

    x = rmag * (cos(Ω)*cos(θ) - sin(Ω)*sin(θ)*cos(i))
    y = rmag * (sin(Ω)*cos(θ) + cos(Ω)*sin(θ)*cos(i))
    z = rmag * (sin(θ)*sin(i))

    h = √(μ*a*(1.0 - e^2))

    ẋ = -μ/h * (cos(Ω)*(sin(θ) + e*sin(ω)) + sin(Ω)*(cos(θ) + e*cos(ω))*cos(i))
    ẏ = -μ/h * (sin(Ω)*(sin(θ) + e*sin(ω)) - cos(Ω)*(cos(θ) + e*cos(ω))*cos(i))
    ż = μ/h * (cos(θ) + e*cos(ω)) * sin(i)

    return [x; y; z; ẋ; ẏ; ż]
end

export koe2USM7
"""
Converts Keplerian Orbital Elements into the Unified State Model Set
Van den Broeck, Michael. "An Approach to Generalizing Taylor Series Integration for Low-Thrust Trajectories." (2017).
https://repository.tudelft.nl/islandora/object/uuid%3A2567c152-ab56-4323-bcfa-b076343664f9

Arugments:
-'u:AbstractVector{<:Number}': Keplerian State Vector [a; e; i; Ω(RAAN); ω(AOP); f(True Anomaly)]
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_USM::Vector{<:Number}': Unified State Model Vector [C; Rf1; Rf2; ϵO1; ϵO2; ϵO3; η0]
*All Angles are in Radians
"""
function koe2USM7(u::AbstractVector{<:Number}, μ::Number)
    a, e, i, Ω, ω, f = u

    ##TODO: NEED TO ADD PARABOLIC CASE
    C = √(μ / (a *(1.0 - e^2)))
        
    R = e*C
    Rf1 = -R*sin(Ω + ω)
    Rf2 = R*cos(Ω + ω)

    u = ω + f

    ϵO1 = sin(i/2) * cos((Ω - u)/2)
    ϵO2 = sin(i/2) * sin((Ω - u)/2)
    ϵO3 = cos(i/2) * sin((Ω + u)/2)
    η0 = cos(i/2) * cos((Ω + u)/2)

    return [C; Rf1; Rf2; ϵO1; ϵO2; ϵO3; η0]

end

export USM72koe
"""
Converts Unified State Model elements into the Keplerian Orbital Element Set
Van den Broeck, Michael. "An Approach to Generalizing Taylor Series Integration for Low-Thrust Trajectories." (2017).
https://repository.tudelft.nl/islandora/object/uuid%3A2567c152-ab56-4323-bcfa-b076343664f9

Arugments:
-'u::AbstractVector{<:Number}': Unified State Model Vector [C; Rf1; Rf2; ϵO1; ϵO2; ϵO3; η0]
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_koe:Vector{<:Number}': Keplerian State Vector [a; e; i; Ω(RAAN); ω(AOP); f(True Anomaly)]
*All Angles are in Radians
"""
function USM72koe(u::AbstractVector{<:Number}, μ::Number)
    C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0 = u

    sinλ = (2.0*ϵO3*η0) / (ϵO3^2 + η0^2)
    cosλ = (η0^2 - ϵO3^2) / (ϵO3^2 + η0^2)
    λ = atan(sinλ, cosλ)

    ve1 = Rf1*cosλ + Rf2*sinλ
    ve2 = C - Rf1*sinλ + Rf2*cosλ

    R = √(Rf1^2 + Rf2^2)
    e = R/C

    a = μ / (2*C*ve2 - (ve1^2 + ve2^2))
    i = acos(1.0 - 2.0*(ϵO1^2 + ϵO2^2))

    if (ϵO1 == 0.0 && ϵO2 == 0.0) || (ϵO3 == 0.0 && η0 == 0.0)
        Ω = 0.0
    else
        Ω = atan((ϵO1*ϵO3 + ϵO2*η0)/√((ϵO1^2 + ϵO2^2)*(η0^2 + ϵO3^2)), (ϵO1*η0 - ϵO2*ϵO3)/√((ϵO1^2 + ϵO2^2)*(η0^2 + ϵO3^2)))
        while Ω < 0.0
            Ω = Ω + 2*π
        end
    end

    if R == 0.0
        ω = 0.0
        f = λ - Ω
        while f < 0.0
            f = f + 2*π
        end
    else
        f = atan(ve1/R, (ve2 - C)/R)
        while f < 0.0
            f = f + 2*π
        end
        ω = λ - Ω - f
        while ω < 0.0
            ω = ω + 2*π
        end
    end

    return [a; e; i; Ω; ω; f]

end

export USM72USM6
"""
Converts USM with Quaternions to USM with Modified Rodrigue Parameters

Arugments:
-'u::AbstractVector{<:Number}': Unified State Model Vector [C; Rf1; Rf2; ϵO1; ϵO2; ϵO3; η0]
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_USM6::Vector{Number}': USM6 State Vector [C; Rf1; Rf2; σ1; σ2; σ3]
"""
function USM72USM6(u::AbstractVector{<:Number}, μ::Number)

    C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0 = u

    σ = EP2MRP([η0; ϵO1; ϵO2; ϵO3])
    
    return [C; Rf1; Rf2; σ]
end

export USM62USM7
"""
Converts USM with Modified Rodrigue Parameters to USM with Quaternions

Arugments:
-'u::AbstractVector{<:Number}': USM6 Vector [C; Rf1; Rf2; σ1; σ2; σ3]
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_USM::Vector{Number}': Unified State Model Vector [C; Rf1; Rf2; ϵO1; ϵO2; ϵO3; η0]
"""
function USM62USM7(u::AbstractArray{<:Number}, μ::Number)

    C, Rf1, Rf2, σ1, σ2, σ3 = u

    η0, ϵO1, ϵO2, ϵO3 = EP2MRP([σ1; σ2; σ3])
    
    return [C; Rf1; Rf2; ϵO1; ϵO2; ϵO3; η0]
end


export USM72USMEM
"""
Converts USM with Quaternions to USM with Exponential Mapping

Arugments:
-'u::AbstractVector{<:Number}': Unified State Model Vector [C; Rf1; Rf2; ϵO1; ϵO2; ϵO3; η0]
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_USMEM::Vector{Number}': USMEM State Vector [C; Rf1; Rf2; a1; a2; a3, Φ]
"""
function USM72USMEM(u::AbstractVector{<:Number}, μ::Number)

    C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0 = u

    Φ = 2.0 * acos(η0)
    a = [ϵO1; ϵO2; ϵO3]
    
    return [C; Rf1; Rf2; Φ.*a]
end

export USMEM2USM7
"""
Converts USM with Exponential Mapping to USM with Quaternions

Arugments:
-'u::AbstractVector{<:Number}': USMEM Vector [C; Rf1; Rf2; a1; a2; a3, Φ]
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_USM::Vector{<:Number}': Unified State Model Vector [C; Rf1; Rf2; ϵO1; ϵO2; ϵO3; η0]
"""
function USMEM2USM7(u::AbstractVector{<:Number}, μ::Number)

    C, Rf1, Rf2, a1, a2, a3 = u

    a = [a1; a2; a3]
    Φ = norm(a)

    ϵ = sin(Φ/2.0)/Φ * a
    η0 = cos(Φ/2.0)
    
    return [C; Rf1; Rf2; ϵ; η0]
end

export koe2ModEq
"""
Converts Keplerian Elements into the Modified Equinoctial Elements

Arugments:
-'u:AbstractVector{<:Number}': Keplerian State Vector [a; e; i; Ω(RAAN); ω(AOP); ν(True Anomaly)]
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_ModEq::Vector{<:Number}': Modified Equinoctial State Vector [p; f; g; h; k; l] 
*All Angles are in Radians
"""
function koe2ModEq(u::AbstractVector{<:Number}, μ::Number)
    a, e, i, Ω, ω, ν = u

    p = a*(1 - e^2)
    f = e*cos(ω + Ω)
    g = e*sin(ω + Ω)
    h = tan(i/2) * cos(Ω)
    k = tan(i/2) * sin(Ω)
    L = Ω + ω + ν

    return [p; f; g; h; k; L]

end

export ModEq2koe
"""
Converts Modified Equinoctial Elements into the Keplerian Elements

Arugments:
-'u:AbstractVector{<:Number}': Modified Equinoctial State Vector [p; f; g; h; k; l] 
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_koe::Vector{<:Number}': Keplerian State Vector [a; e; i; Ω(RAAN); ω(AOP); ν(True Anomaly)]
*All Angles are in Radians
"""
function ModEq2koe(u::AbstractVector{<:Number}, μ::Number)
    p, f, g, h, k, L = u

    a = p/(1 - f^2 - g^2)
    e = √(f^2 + g^2)
    i = atan(2*√(h^2 + k^2), 1-h^2-k^2)
    Ω = atan(k, h)
    ω = atan(g*h - f*k, f*h + g*k)
    ν = L - atan(g, f)

    return [a; e; i; Ω; ω; ν]

end

export koe2ModEqN
"""
Converts Keplerian Elements into the Modified Equinoctial Elements with Mean Motion

Arugments:
-'u:AbstractVector{<:Number}': Keplerian State Vector [a; e; i; Ω(RAAN); ω(AOP); ν(True Anomaly)]
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_ModEq::Vector{<:Number}': Modified Equinoctial State Vector [n; f; g; h; k; l] 
*All Angles are in Radians
"""
function koe2ModEqN(u::AbstractVector{<:Number}, μ::Number)
    a, e, i, Ω, ω, ν = u

    n = √(μ/(a^3))
    f = e*cos(ω + Ω)
    g = e*sin(ω + Ω)
    h = tan(i/2) * cos(Ω)
    k = tan(i/2) * sin(Ω)
    L = Ω + ω + ν

    return [n; f; g; h; k; L]

end

export ModEqN2koe
"""
Converts Modified Equinoctial Elements with Mean Motion into the Keplerian Elements

Arugments:
-'u:AbstractVector{<:Number}': Modified Equinoctial State Vector [p; f; g; h; k; l] 
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_koe::Vector{<:Number}': Keplerian State Vector [a; e; i; Ω(RAAN); ω(AOP); ν(True Anomaly)]
*All Angles are in Radians
"""
function ModEqN2koe(u::AbstractVector{<:Number}, μ::Number)
    n, f, g, h, k, L = u

    a = ∛(μ/(n^2))
    e = √(f^2 + g^2)
    i = atan(2*√(h^2 + k^2), 1-h^2-k^2)
    Ω = atan(k, h)
    ω = atan(g*h - f*k, f*h + g*k)
    ν = L - atan(g, f)

    return [a; e; i; Ω; ω; ν]

end

export cart2Mil
"""
Converts Cartesian State Vector into the Milankovich State

Arugments:
-'u::AbstractVector{<:Number}': Cartesian State Vector [x; y; z; ẋ; ẏ; ż] 
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_Mil::Vector{Number}': Milankovich State Vector [H; e; L] 
"""
function cart2Mil(u::AbstractVector{<:Number}, μ::Number)
    r = u[1:3]
    v = u[4:6]

    H = cross(r, v)
    e = cross(v/μ, H) - r/norm(r)

    x̂ = [1.; 0.; 0.]
    ŷ = [0.; 1.; 0.]
    ẑ = [0.; 0.; 1.]

    if abs(H[1]) > eps(T) || abs(H[2]) > eps(T) 
        nΩ = cross(ẑ, H)/norm(cross(ẑ, H))
        Ω = atan(dot(ŷ, nΩ), dot(x̂, nΩ))
        nΩperp = cross(H, nΩ/norm(H))

        if norm(e) > eps(T)
            ω = atan(dot(e, nΩperp), dot(e, nΩ))

            eperp = cross(H, e)/(norm(H)*norm(e))
            f = atan(dot(r, eperp), dot(r, e/norm(e)))

            L = Ω + ω + f
        else
            L = Ω + atan(dot(r, nΩperp), dot(r, nΩ))
        end
    else
        if norm(e) > eps(T)

            ω_hat = atan(dot(e, ŷ), dot(e, x̂))
            eperp = cross(H, e)/(norm(H)*norm(e))

            f = atan(dot(r, eperp), dot(r, e/norm(e)))

            L = ω_hat + f
        else
            L = atan(dot(r, ŷ), dot(r, x̂))
        end
    end

    return [H; e; L]
end


export Mil2cart
"""
Converts Milankovich State Vector into the Cartesian State

Arugments:
-'u::AbstractVector{<:Number}': Milankovich State Vector [H; e; L] 
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_cart::Vector{<:Number}': Cartesian State Vector [x; y; z; ẋ; ẏ; ż] 
"""
function Mil2cart(u::AbstractVector{<:Number}, μ::Number)

    H = u[1:3]
    e = u[4:6]
    L = u[7]

    x̂ = [1.; 0.; 0.]
    ŷ = [0.; 1.; 0.]
    ẑ = [0.; 0.; 1.]

    if abs(H[1]) > eps(T) || abs(H[2]) > eps(T) 
        nΩ = cross(ẑ, H)/norm(cross(ẑ, H))
        Ω = atan(dot(ŷ, nΩ), dot(x̂, nΩ))
        nΩperp = cross(H, nΩ/norm(H))

        if norm(e) > eps(T)
            ω = atan(dot(e, nΩperp), dot(e, nΩ))

            f = L - (Ω + ω)
        else
            f = L - Ω
        end
    else
        if norm(e) > eps(T)

            ω_hat = atan(dot(e, ŷ), dot(e, x̂))

            f = L - ω_hat
        else
            f = L
        end
    end

    p = norm(H)^2/μ
    rmag = p/(1.0 + norm(e)*cos(f))

    if abs(H[1]) > eps(T) || abs(H[2]) > eps(T)
        if norm(e) > eps(T)
            ê = e/norm(e)
            eperp = cross(H, e)/(norm(H)*norm(e))
        else
            ê = nΩ
            eperp = nΩperp
        end
    else
        if norm(e) > eps(T)
            ê = e/norm(e)
            eperp = cross(H, e)/(norm(H)*norm(e))
        else
            ê = x̂
            eperp = ŷ
        end
    end

    r = rmag * (cos(f)*ê + sin(f)*eperp)
    v = √(μ/p) * (-sin(f)*ê + (norm(e) + cos(f))*eperp)

    return [r; v]

end

export koeM2cart
"""
Converts Alternative Keplerian State Vector into the Cartesian State

Arugments:
-'u::AbstractVector{<:Number}': Alternative Keplerian State Vector [a; e; i; Ω(RAAN); ω(AOP); M(Mean Anomaly)]
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_cart::Vector{<:Number}': Cartesian State Vector [x; y; z; ẋ; ẏ; ż] 
"""
function koeM2cart(u::AbstractVector{<:Number}, μ::Number)
    a, e, i, Ω, ω, M = u

    f = KeplerSolver(M, e)

    return koe2cart([a; e; i; Ω; ω; f], μ)
end

export cart2koeM
"""
Converts Cartesian State Vector into the Alternative Keplerian State Vector

Arugments:
-'u::AbstractVector{<:Number}':  Cartesian State Vector [x; y; z; ẋ; ẏ; ż]  
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_cart::Vector{<:Number}': Alternative Keplerian State Vector [a; e; i; Ω(RAAN); ω(AOP); M(Mean Anomaly)]
"""
function cart2koeM(u::AbstractVector{<:Number}, μ::Number)
    
    (a, e, i, Ω, ω, f) = cart2koe(u, μ)

    M = trueAnomaly2MeanAnomaly(f, e)

    return [a, e, i, Ω, ω, M]

end

export cart2cylind
"""
Computes the Cylindrical Orbital Elements from a Cartesian Set

Arguments:
-'u::AbstractVector{<:Number}': Cartesian State Vector [x; y; z; ẋ; ẏ; ż]  
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_cylind::Vector{<:Number}': Cylindrical Orbital Element Vector [r; θ; z; ṙ; θdot; ż]
*All Angles are in Radians
"""
function cart2cylind(u::AbstractVector{<:Number}, μ::Number)
    x, y, z, ẋ, ẏ, ż = u

    r = norm([x; y; z])
    θ = atan(y, x)

    ṙ = (x*ẋ + y*ẏ) / r
    θdot = (x*ẏ - y*ẋ) / r

    return [r; θ; z; ṙ; θdot; ż]
end

export cylind2cart
"""
Computes the Cartesian Orbital Elements from a Cylindrical Set

Arguments:
-'u::AbstractVector{<:Number}': Cylindrical State Vector [r; θ; z; ṙ; θdot; ż]
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_cart::Vector{<:Real}': Cartesian Orbital Element Vector [x; y; z; ẋ; ẏ; ż]  
*All Angles are in Radians
"""
function cylind2cart(u::AbstractVector{<:Number}, μ::Number)
    r, θ, z, ṙ, θdot, ż = u

    x = r * cos(θ)
    y = r * sin(θ)

    ẋ = ṙ*cos(θ) - θdot*sin(θ)
    ẏ = ṙ*sin(θ) + θdot*cos(θ)

    return [x; y; z; ẋ; ẏ; ż]
end

export cart2sphere
"""
Computes the Spherical Orbital Elements from a Spherical Set

Arguments:
-'u::AbstractVector{<:Number}': Cartesian State Vector [x; y; z; ẋ; ẏ; ż]  
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_sphere::Vector{<:Number}': Spherical Orbital Element Vector [r; θ; ϕ; ṙ; θdot; ϕdot]
*All Angles are in Radians
"""
function cart2sphere(u::AbstractVector{<:Real}, μ::Number)
    x, y, z, ẋ, ẏ, ż = u

    r = norm([x; y; z])
    θ = atan(y, x)
    ϕ = asin(z)

    ṙ = (x*ẋ + y*ẏ + z*ż) / r
    θdot = (x*ẏ - ẋ*y) / norm([x; y])
    ϕdot = (ż*x^2 - ẋ*z*x + ż*y^2 - ẏ*z*y) / (√(1.0 - (z^2)/r) * r^1.5)

    return [r; θ; ϕ; ṙ; θdot; ϕdot]
end

export sphere2cart
"""
Computes the Cartesian Orbital Elements from a Spherical Set

Arguments:
-'u::AbstractVector{<:Number}': Spherical Orbital Element Vector [r; θ; ϕ; ṙ; θdot; ϕdot]
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_cart::Vector{<:Number}': Cartesian Orbital Element Vector [x; y; z; ẋ; ẏ; ż]
*All Angles are in Radians
"""
function sphere2cart(u::AbstractVector{<:Number}, μ::Number)
    r, θ, ϕ, ṙ, θdot, ϕdot = u

    x = r * cos(θ) * cos(ϕ)
    y = r * sin(θ) * cos(ϕ)
    z = r * sin(ϕ)

    ẋ = ṙ*cos(θ)*cos(ϕ) - r*θdot*sin(θ)*cos(ϕ) - ϕdot*r*cos(θ)*sin(ϕ)
    ẏ = ṙ*sin(θ)*cos(ϕ) - r*ϕdot*sin(θ)*sin(ϕ) + r*θdot*cos(θ)*cos(ϕ)
    ż = ṙ*sin(ϕ) + r*ϕdot*cos(ϕ)

    return [x; y; z; ẋ; ẏ; ż]
end

export cart2delauney
"""
Computes the Delauney Orbital Elements from a Cartesian Set

Arguments:
-'u::AbstractVector{<:Number}': Cartesian Orbital Element Vector [x; y; z; ẋ; ẏ; ż]
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_cart::Vector{<:Number}': Delauney Orbital Element Vector [L; G; H; M; ω; Ω]
*All Angles are in Radians
"""
function cart2delauney(u::AbstractVector{<:Number}, μ::Number)

    (a, e, i, Ω, ω, M) = cart2koeM(u, μ)

    L = √(μ * a)
    G = L*√(1.0 - e^2)
    H = G * cos(i)

    return [L; G; H; M; ω; Ω]

end

export delauney2cart
"""
Computes the Cartesian Orbital Elements from a Delauney Set
Laskar, Jacques. "Andoyer construction for Hill and Delaunay variables." Celestial Mechanics and Dynamical Astronomy 128.4 (2017): 475-482.

Arguments:
-'u::AbstractVector{<:Number}': Delauney Orbital Element Vector [L; G; H; M; ω; Ω]
-'μ::Number': Standard Graviational Parameter of Central Body

Returns:
-'u_cart::Vector{<:Number}': Cartesian Orbital Element Vector [x; y; z; ẋ; ẏ; ż]
*All Angles are in Radians
"""
function delauney2cart(u::AbstractVector{<:Number}, μ::Number)

    L, G, H, M, ω, Ω = u

    a = L^2/μ
    e = √(1.0 - (G/L)^2)
    i = acos(H/G)

    return koeM2cart([a; e; i; Ω; ω; M], μ)

end
