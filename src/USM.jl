export USM7
"""
    USM7{T} <: AstroCoord

Unified State Model Orbital Elements. 7D parameterziation of the orbit using Velocity Hodograph and Quaternions
C - 
Rf1 - 
Rf2 -
ϵO1 -
ϵO2 -
ϵO3 -
η0 -

Constructors
USM7(C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0)
USM7(X::AbstractArray)
USM7(X::AstroCoord, μ::Number)

"""
struct USM7{T} <: AstroCoord{7, T}
    C::T
    Rf1::T
    Rf2::T
    ϵO1::T
    ϵO2::T
    ϵO3::T
    η0::T
    @inline USM7{T}(C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0) where T = new{T}(C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0)
    @inline USM7{T}(p::USM7) where T = new{T}(p.C, p.Rf1, p.Rf2, p.ϵO1, p.ϵO2, p.ϵO3, p.η0)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
USM7(X::AbstractArray{T, 1}) where T = USM7{T}(X...)
USM7(C::CT, Rf1::R1, Rf2::R2, ϵO1::E1, ϵO2::E2, ϵO3::E3, η0::N) where{CT,R1,R2,E1,E2,E3,N} = USM7{promote_type(CT,R1,R2,E1,E2,E3,N)}(C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0)
(::Type{U7})(g::StaticVector) where U7<:USM7 = U7(g[1], g[2], g[3], g[4], g[5], g[6], g[7])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::USM7) = SVector{7}(g.C, g.Rf1, g.Rf2, g.ϵO1, g.ϵO2, g.ϵO3, g.η0)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
Base.one(::Type{U7}) where U7 <: USM7 = U7(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::USM7, i::Int)
    if i == 1 return p.C
    elseif i == 2 return p.Rf1
    elseif i == 3 return p.Rf2
    elseif i == 4 return p.ϵO1
    elseif i == 5 return p.ϵO2
    elseif i == 6 return p.ϵO3
    elseif i == 7 return p.η0
    else throw(BoundsError(r,i))
    end
end

#TODO: CITE PAPER
#TODO: I DONT THINK THIS BELONGS HERE
export USM7_ODE
"""
Unified State Model ODE System

Arguments:
-'u::AbstractArray{AbstractFloat, 1}': USM State Vector [C; Rf1; Rf2; ϵO1; ϵO2; ϵO3; η0]
-'p::AbstractArray{AbstractFloat, 1}': Parameter Vector
-'t::AbstractFloat': Time
-'accel::Function': Function computing perturbing acceleration in the orbital frame 

Returns:
-'du::AbstractArray{AbstractFloat, 1}': USM ODE Change in State Vector [dC; dRf1; dRf2; dϵO1; dϵO2; dϵO3; dη0]
"""
function USM7_ODE(u::AbstractVector{<:Number}, p::AbstractVector{<:Number}, t::Number, accel::Function)

    C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0 = u
    μ = p[1]

    sinλ = (2*ϵO3*η0) / (ϵO3^2 + η0^2)
    cosλ = (η0^2 - ϵO3^2) / (ϵO3^2 + η0^2)
        
    l = (ϵO1*ϵO3 - ϵO2*η0) / (ϵO3^2 + η0^2)

    ve2 = C - Rf1*sinλ + Rf2*cosλ

    ω3 = (C*ve2^2)/μ

    ρ = C / ve2

    input = [C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0, sinλ, cosλ, l, ve2, ω3, ρ, 0.0, 0.0, 0.0]
        
    fe = accel(input, p, t)

    ω1 = fe[3] / ve2    

    dC = -ρ*fe[2]
    dRf1 = fe[1]*cosλ - fe[2]*(1. + ρ)*sinλ - fe[3]*l*(Rf2/ve2)
    dRf2 = fe[1]*sinλ + fe[2]*(1. + ρ)*cosλ + fe[3]*l*(Rf1/ve2)
    dϵO1 = .5*(ω3*ϵO2 + ω1*η0)
    dϵO2 = .5*(-ω3*ϵO1 + ω1*ϵO3)
    dϵO3 = .5*(-ω1*ϵO2 + ω3*η0)
    dη0 = .5*(-ω1*ϵO1 - ω3*ϵO3)

    return [dC; dRf1; dRf2; dϵO1; dϵO2; dϵO3; dη0]
    
end

export USM6
"""
    USM6{T} <: AstroCoord

Unified State Model Orbital Elements. 6D parameterziation of the orbit using Velocity Hodograph and MRP's
C - 
Rf1 - 
Rf2 -
σ1 -
σ2 -
σ3 -

Constructors
USM6(C, Rf1, Rf2, σ1, σ2, σ3)
USM6(X::AbstractArray)
USM6(X::AstroCoord, μ::Number)

"""
struct USM6{T} <: AstroCoord{6, T}
    C::T
    Rf1::T
    Rf2::T
    σ1::T
    σ2::T
    σ3::T
    @inline USM6{T}(C, Rf1, Rf2, σ1, σ2, σ3) where T = new{T}(C, Rf1, Rf2, σ1, σ2, σ3)
    @inline USM6{T}(p::USM6) where T = new{T}(p.C, p.Rf1, p.Rf2, p.σ1, p.σ2, p.σ3)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
USM6(X::AbstractArray{T, 1}) where T = USM6{T}(X...)
USM6(C::CT, Rf1::R1, Rf2::R2, σ1::S1, σ2::S2, σ3::S3) where{CT,R1,R2,S1,S2,S3} = USM6{promote_type(CT,R1,R2,S1,S2,S3)}(C, Rf1, Rf2, σ1, σ2, σ3)
(::Type{U6})(g::StaticVector) where U6<:USM6 = U6(g[1], g[2], g[3], g[4], g[5], g[6])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::USM6) = SVector{6}(g.C, g.Rf1, g.Rf2, g.σ1, g.σ2, g.σ3)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
Base.one(::Type{U6}) where U6 <: USM6 = U6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::USM6, i::Int)
    if i == 1 return p.C
    elseif i == 2 return p.Rf1
    elseif i == 3 return p.Rf2
    elseif i == 4 return p.σ1
    elseif i == 5 return p.σ2
    elseif i == 6 return p.σ3
    else throw(BoundsError(r,i))
    end
end

export USM6_ODE
"""
USM6 ODE System

Arguments:
-'u::AbstractArray{AbstractFloat, 1}': USM6 State Vector [C; Rf1; Rf2; σ1; σ2; σ3]
-'p::AbstractArray{AbstractFloat, 1}': Parameter Vector
-'t::AbstractFloat': Time
-'accel::Function': Function computing perturbing acceleration in the orbital frame 

Returns:
-'du::AbstractArray{AbstractFloat, 1}': USM6 ODE Change in State Vector [C; Rf1; Rf2; σ1; σ2; σ3]
"""
function USM6_ODE(u::AbstractVector{<:Number}, p::AbstractVector{<:Number}, t::Number, accel::Function)

    C, Rf1, Rf2, σ1, σ2, σ3 = u
    σ = norm([σ1; σ2; σ3])

    if σ > 1.0
        u[4:6] = -[σ1; σ2; σ3]/σ
    end

    C, Rf1, Rf2, σ1, σ2, σ3 = u
    σ = norm([σ1; σ2; σ3])

    μ = p[1]

    sinλ = (4.0*σ3*(1.0 - σ^2)) / (4.0*σ3^2 + (1 - σ^2)^2)
    cosλ = ((1.0 - σ^2)^2 - 4.0*σ3^2) / (4.0*σ3^2 + (1 - σ^2)^2)
        
    l = (ϵO1*ϵO3 - ϵO2*η0) / (ϵO3^2 + η0^2)

    ve2 = C - Rf1*sinλ + Rf2*cosλ

    ω3 = (C*ve2^2)/μ

    ρ = C / ve2

    input = [C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0, sinλ, cosλ, l, ve2, ω3, ρ, 0.0, 0.0, 0.0]
        
    fe = accel(input, p, t)

    ω1 = fe[3] / ve2    

    dC = -ρ*fe[2]
    dRf1 = fe[1]*cosλ - fe[2]*(1. + ρ)*sinλ - fe[3]*l*(Rf2/ve2)
    dRf2 = fe[1]*sinλ + fe[2]*(1. + ρ)*cosλ + fe[3]*l*(Rf1/ve2)
    dσ1 = .25*((1.0 - σ^2 + 2.0*σ1^2)*ω1 + 2.0*(σ1*σ3 + σ2)*ω3)
    dσ2 = .25*(2.0*(σ2*σ1 + σ3)*ω1 + 2.0*(σ2*σ3 - σ1)*ω3)  
    dσ3 = .25*(2.0*(σ3*σ1 - σ2)*ω1 + (1.0 - σ^2 + 2.0*σ3^2)*ω3)

    return [dC; dRf1; dRf2; dσ1; dσ2; dσ3]
    
end

export USMEM
"""
    USMEM{T} <: AstroCoord

Unified State Model Orbital Elements. 6D parameterziation of the orbit using Velocity Hodograph and Exponential Mapping
C - 
Rf1 - 
Rf2 -
a1 -
a2 -
a3 -

Constructors
USMEM(C, Rf1, Rf2, a1, a2, a3)
USMEM(X::AbstractArray)
USMEM(X::AstroCoord, μ::Number)

"""
struct USMEM{T} <: AstroCoord{6, T}
    C::T
    Rf1::T
    Rf2::T
    a1::T
    a2::T
    a3::T
    @inline USMEM{T}(C, Rf1, Rf2, a1, a2, a3) where T = new{T}(C, Rf1, Rf2, a1, a2, a3)
    @inline USMEM{T}(p::USMEM) where T = new{T}(p.C, p.Rf1, p.Rf2, p.a1, p.a2, p.a3)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
USMEM(X::AbstractArray{T, 1}) where T = USMEM{T}(X...)
USMEM(C::CT, Rf1::R1, Rf2::R2, a1::A1, a2::A2, a3::A3) where{CT,R1,R2,A1,A2,A3} = USMEM{promote_type(CT,R1,R2,A1,A2,A3)}(C, Rf1, Rf2, a1, a2, a3)
(::Type{UE})(g::StaticVector) where UE<:USMEM = UE(g[1], g[2], g[3], g[4], g[5], g[6])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::USMEM) = SVector{6}(g.C, g.Rf1, g.Rf2, g.a1, g.a2, g.a3)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
Base.one(::Type{UE}) where UE <: USMEM = UE(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::USMEM, i::Int)
    if i == 1 return p.C
    elseif i == 2 return p.Rf1
    elseif i == 3 return p.Rf2
    elseif i == 4 return p.a1
    elseif i == 5 return p.a2
    elseif i == 6 return p.a3
    else throw(BoundsError(r,i))
    end
end

#TODO: FINISH EM
export USMEM_ODE
"""
USM6 ODE System using Exponential Mapping

Arguments:
-'u::AbstractArray{AbstractFloat, 1}': USM6 State Vector [C; Rf1; Rf2; σ1; σ2; σ3]
-'p::AbstractArray{AbstractFloat, 1}': Parameter Vector
-'t::AbstractFloat': Time
-'accel::Function': Function computing perturbing acceleration in the orbital frame 

Returns:
-'du::AbstractArray{AbstractFloat, 1}': USM6 ODE Change in State Vector [C; Rf1; Rf2; σ1; σ2; σ3]
"""
function USMEM_ODE(u::AbstractVector{<:Number}, p::AbstractVector{<:Number}, t::Number, accel::Function)

    C, Rf1, Rf2, a1, a2, a3 = u
    μ = p[1]

    Φ = norm([a1; a2; a3])

    u_USM = USMEM2USM(u, μ)

    C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0 = u_USM

    sinλ = (2*ϵO3*η0) / (ϵO3^2 + η0^2)
    cosλ = (η0^2 - ϵO3^2) / (ϵO3^2 + η0^2)
        
    l = (ϵO1*ϵO3 - ϵO2*η0) / (ϵO3^2 + η0^2)

    ve2 = C - Rf1*sinλ + Rf2*cosλ

    ω3 = (C*ve2^2)/μ

    ρ = C / ve2

    input = [C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0, sinλ, cosλ, l, ve2, ω3, ρ, 0.0, 0.0, 0.0]
        
    fe = accel(input, p, t)

    ω1 = fe[3] / ve2    

    ax = skew_sym([a1; a2; a3])
    ω = [ω1; 0.0; ω3]

    dC = -ρ*fe[2]
    dRf1 = fe[1]*cosλ - fe[2]*(1. + ρ)*sinλ - fe[3]*l*(Rf2/ve2)
    dRf2 = fe[1]*sinλ + fe[2]*(1. + ρ)*cosλ + fe[3]*l*(Rf1/ve2)
    if abs(Φ) < eps(T)^.25 
        da = (I + ax/2.0 + (1.0/(Φ^2))*(1.0 - (Φ/2.0)*cot(Φ/2.0))*ax*ax)*ω
    else
        da = .5*(((12.0 - Φ^2)/6)*ω - cross(ω, a) - dot(ω, a)*((60.0 + Φ^2)/360.0)*a)
    end

    return [dC; dRf1; dRf2; da]
    
end