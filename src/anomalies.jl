#TODO: DEFINE ADJOINT (IMPLICIT?)
export KeplerSolver
"""
Solves for True Anomaly given the Mean Anomaly and eccentricity of an Orbit

Arguments:
-'M::Number': Mean Anomaly of the orbit
-'e::Number': Eccentricity of the orbit

Returns
-'f::Number': True Anomaly of the orbit
*All angles are in radians
"""
function KeplerSolver(M::T, e::Number; tol=10*eps(T)) where {T<:Number}

    if e < 1.0
        E_guess = M
        fE = M - (E_guess - e*sin(E_guess))
        dfE = (e*cos(E_guess) - 1)

        while abs(fE) > tol
            E_guess -= fE/dfE
            fE = M - (E_guess - e*sin(E_guess))
            dfE = (e*cos(E_guess) - 1)
        end
        E_guess = mod(E_guess, 2*π)

        f = 2*atan(√(1+e)*sin(E_guess/2), √(1-e)*cos(E_guess/2))

        if f < 0.0
            f = 2*π - abs(f)
        end
    else
        F_guess = M
        fF = M + (F_guess - e*sinh(F_guess))
        dfF = (1 - e*cosh(F_guess))

        while abs(fF) > tol
            F_guess -= fF/dfF
            fF = M + (F_guess - e*sinh(F_guess))
            dfF = (1 - e*cosh(F_guess))
        end

        F_guess = mod(F_guess, 2*π)

        f = 2.0*atan(√(1+e)*sinh(F_guess/2), √(e-1)*cosh(F_guess/2))

        if f < 0.0
            f = 2*π - abs(f)
        end
    end

    return f
end

export trueAnomaly2MeanAnomaly
"""
Converts the True Anomaly into the Mean Anomaly

Arguments:
-'f::Number': True Anomaly of the orbit
-'e::Number': Eccentricity of the orbit

Returns
-'M::Number': Mean Anomaly of the orbit
*All angles are in radians
"""
function trueAnomaly2MeanAnomaly(f::Number, e::Number)
    
    E = atan((sin(f)*√(1-e^2))/(1.0+e*cos(f)), (e+cos(f))/(1.0+e*cos(f)))

    M = E - e*sin(E)

    M = M % (2.0*π)

    return M

end 

export trueAnomaly2EccentricAnomaly
"""
Converts the True Anomaly into the Mean Anomaly

Arguments:
-'f::Number': True Anomaly of the orbit
-'e::Number': Eccentricity of the orbit

Returns
-'E::Number': Eccentric Anomaly of the orbit
*All angles are in radians
"""
function trueAnomaly2EccentricAnomaly(f::Number, e::Number)
    
    E = atan((sin(f)*√(1-e^2))/(1.0+e*cos(f)), (e+cos(f))/(1.0+e*cos(f)))

end 

export eccentricAnomaly2MeanAnomaly
"""
Converts the True Anomaly into the Mean Anomaly

Arguments:
-'E::Number': Eccentric Anomaly of the orbit
-'e::Number': Eccentricity of the orbit

Returns
-'M::Number': Mean Anomaly of the orbit
*All angles are in radians
"""
function eccentricAnomaly2MeanAnomaly(E::Number, e::Number)

    M = E - e*sin(E)

    M = M % (2.0*π)

    return M

end 

export eccentricAnomaly2TrueAnomaly
"""
Converts the Eccentric Anomaly into the True Anomaly

Arguments:
-'E::Number': Eccentric Anomaly of the orbit
-'e::Number': Eccentricity of the orbit

Returns
-'f::Number': True Anomaly of the orbit
*All angles are in radians
"""
function eccentricAnomaly2TrueAnomaly(E::Number, e::Number)
    
    return 2*atan(√(1.0+e)*sin(E/2.0), √(1.0-e)*cos(E/2.0))

end 

export meanAnomaly2TrueAnomaly
"""
Converts the Mean Anomaly into the True Anomaly

Arguments:
-'M::Number': Mean Anomaly of the orbit
-'e::Number': Eccentricity of the orbit

Returns
-'f::Number': Mean Anomaly of the orbit
*All angles are in radians
"""
function meanAnomaly2TrueAnomaly(M::T, e::Number; tol=10*eps(T)) where {T<:Number}
    
    return KeplerSolver(M, e; tol=tol)

end 

export meanAnomaly2EccentricAnomaly
"""
Converts the Mean Anomaly into the Eccentric Anomaly

Arguments:
-'M::Number': Mean Anomaly of the orbit
-'e::Number': Eccentricity of the orbit

Returns
-'E::Number': Eccentric Anomaly of the orbit
*All angles are in radians
"""
function meanAnomaly2EccentricAnomaly(M::T, e::Number; tol=10*eps(T)) where {T<:Number}
    
    return trueAnomaly2EccentricAnomaly(KeplerSolver(M, e; tol=tol), e)

end 