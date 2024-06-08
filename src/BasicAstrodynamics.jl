export meanMotion
"""
Computes the Keplerian Mean Motion About a Central Body

Arguments:
- 'a::Number': Semi-Major Axis
- 'μ::Number': Standard Graviational Parameter of Central Body

Returns:
- 'n::Number': Orbital Mean Motion
"""
function meanMotion(a::Number, μ::Number)

    return √(μ/(a^3))
end


export orbitalPeriod
"""
Computes the Keplerian Orbital Period About a Central Body

Arguments:
- 'a::Number': Semi-Major Axis
- 'μ::Number': Standard Graviational Parameter of Central Body

Returns:
- 'T::Number': Orbital Period
"""
function orbitalPeriod(a::Number, μ::Number)

    return 2*π / √(μ/(a^3))
end

export orbitalNRG
"""
Computes the Keplerian Orbital Energy

Arguments:
-'a::Number': AstroCoord State Vector
-'μ::Number': Standard Graviational Parameter of Central Body

Returns
-'NRG::Number' - Orbital energy 

"""
function orbitalNRG(a::Number, μ::Number)

    return -μ/(2.0*a)

end

export angularMomentumVector
"""
Computes the Instantaneous Angular Velocity Vector

Arguments:
-'u::AbstractVector{<:Number}': Cartesian State Vector [x; y; z; ẋ; ẏ; ż]

Returns
-'angular_momentum::Vector{<:Number}' - 3-Dimensional Angular Momemtum Vector
"""
function angularMomentumVector(u::AbstractVector{<:Number})

    return cross(@view(u[1:3]), @view(u[4:6]))

end

export angularMomentumQuantity
"""
Computes the Instantaneous Angular Velocity Quantity

Arguments:
-'u::AbstractVector{<:Number}': Cartesian State Vector [x; y; z; ẋ; ẏ; ż]

Returns
-'angular_momentum::Number' - Norm of Angular Momentum Vector
"""
function angularMomentumQuantity(u::AbstractVector{<:Number})

    return norm(angularMomentumVector(u))

end  

export meanMotion
"""
Computes the Keplerian Mean Motion About a Central Body

Arguments:
- 'X::AstroCoord': Astro Coordinate 
- 'μ::Number': Standard Graviational Parameter of Central Body

Returns:
- 'n::Number': Orbital Mean Motion
"""
function meanMotion(X::AstroCoord, μ::Number)

    kep = Keplerian(X, μ)

    return meanMotion(kep.a, μ)
    
end

export orbitalPeriod
"""
Computes the Keplerian Orbital Period About a Central Body

Arguments:
- 'X::AstroCoord': Astro Coordinate 
- 'μ::Number': Standard Graviational Parameter of Central Body

Returns:
- 'T::Number': Orbital Period
"""
function orbitalPeriod(X::AstroCoord, μ::Number)

    kep = Keplerian(X, μ)

    return orbitalPeriod(kep.a, μ)
    
end

export orbitalNRG
"""
Computes the Keplerian Orbital Energy

Arguments:
-'X::AstroCoord': Astro Coordinate 
-'μ::Number': Standard Graviational Parameter of Central Body

Returns
-'NRG::Number' - Orbital energy 

"""
function orbitalNRG(X::AstroCoord, μ::Number)

    kep = Keplerian(X, μ)

    return orbitalNRG(kep.a, μ)

end

export angularMomentumVector
"""
Computes the Instantaneous Angular Velocity Vector

Arguments:
-'X::AstroCoord': Astro Coordinate 
-'μ::Number': Standard Graviational Parameter of Central Body

Returns
-'angular_momentum::Vector{<:Number}' - 3-Dimensional Angular Momemtum Vector
"""
function angularMomentumVector(X::AstroCoord, μ::Number)

    cart = Cartesian(X, μ)

    return angularMomentumVector(params(cart))
    
end

export angularMomentumQuantity
"""
Computes the Instantaneous Angular Velocity Quantity

Arguments:
-'X::AstroCoord': Astro Coordinate 
-'μ::Number': Standard Graviational Parameter of Central Body

Returns
-'angular_momentum::Number' - Norm of Angular Momentum Vector
"""
function angularMomentumQuantity(X::AstroCoord, μ::Number)

    cart = Cartesian(X, μ)

    return angularMomentumQuantity(params(cart))

end