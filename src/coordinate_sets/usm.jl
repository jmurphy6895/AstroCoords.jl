export USM7
"""
    USM7{T} <: AstroCoord

Unified State Model Orbital Elements. 7D parameterziation of the orbit using Velocity Hodograph and Quaternions
C - Velocity Hodograph Component Normal to the Radial Vector Laying in the Orbital Plane
Rf1 - Velocity Hodograph Component 90 degrees ahead of the Eccentricity Vector - Along the Intermediate Rotating Frame X-Axis
Rf2 - Velocity Hodograph Component 90 degrees ahead of the Eccentricity Vector - Along the Intermediate Rotating Frame Y-Axis
ϵO1 - First Imaginary Quaternion Component
ϵO2 - Second Imaginary Quaternion Component
ϵO3 - Third Imaginary Quaternion Component
η0 - Real Quaternion Component

Constructors
USM7(C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0)
USM7(X::AbstractArray)
USM7(X::AstroCoord, μ::Number)

"""
struct USM7{T} <: AstroCoord{7,T}
    C::T
    Rf1::T
    Rf2::T
    ϵO1::T
    ϵO2::T
    ϵO3::T
    η0::T
    @inline USM7{T}(C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0) where {T} =
        new{T}(C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0)
    @inline USM7{T}(p::USM7) where {T} =
        new{T}(p.C, p.Rf1, p.Rf2, p.ϵO1, p.ϵO2, p.ϵO3, p.η0)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
USM7(X::AbstractArray{T,1}) where {T} = USM7{T}(X...)
function USM7(
    C::CT, Rf1::R1, Rf2::R2, ϵO1::E1, ϵO2::E2, ϵO3::E3, η0::N
) where {CT,R1,R2,E1,E2,E3,N}
    return USM7{promote_type(CT, R1, R2, E1, E2, E3, N)}(C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0)
end
function (::Type{U7})(g::StaticVector) where {U7<:USM7}
    return U7(g[1], g[2], g[3], g[4], g[5], g[6], g[7])
end

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::USM7) = SVector{7}(g.C, g.Rf1, g.Rf2, g.ϵO1, g.ϵO2, g.ϵO3, g.η0)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
Base.one(::Type{U7}) where {U7<:USM7} = U7(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::USM7, i::Int)
    if i == 1
        return p.C
    elseif i == 2
        return p.Rf1
    elseif i == 3
        return p.Rf2
    elseif i == 4
        return p.ϵO1
    elseif i == 5
        return p.ϵO2
    elseif i == 6
        return p.ϵO3
    elseif i == 7
        return p.η0
    else
        throw(BoundsError(r, i))
    end
end

export USM6
"""
    USM6{T} <: AstroCoord

Unified State Model Orbital Elements. 6D parameterziation of the orbit using Velocity Hodograph and MRP's
C - Velocity Hodograph Component Normal to the Radial Vector Laying in the Orbital Plane
Rf1 - Velocity Hodograph Component 90 degrees ahead of the Eccentricity Vector - Along the Intermediate Rotating Frame X-Axis
Rf2 - Velocity Hodograph Component 90 degrees ahead of the Eccentricity Vector - Along the Intermediate Rotating Frame Y-Axis
σ1 - First Modified Rodriguez Parameter 
σ2 - Second Modified Rodriguez Parameter 
σ3 - Third Modified Rodriguez Parameter 

Constructors
USM6(C, Rf1, Rf2, σ1, σ2, σ3)
USM6(X::AbstractArray)
USM6(X::AstroCoord, μ::Number)

"""
struct USM6{T} <: AstroCoord{6,T}
    C::T
    Rf1::T
    Rf2::T
    σ1::T
    σ2::T
    σ3::T
    @inline USM6{T}(C, Rf1, Rf2, σ1, σ2, σ3) where {T} = new{T}(C, Rf1, Rf2, σ1, σ2, σ3)
    @inline USM6{T}(p::USM6) where {T} = new{T}(p.C, p.Rf1, p.Rf2, p.σ1, p.σ2, p.σ3)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
USM6(X::AbstractArray{T,1}) where {T} = USM6{T}(X...)
function USM6(C::CT, Rf1::R1, Rf2::R2, σ1::S1, σ2::S2, σ3::S3) where {CT,R1,R2,S1,S2,S3}
    return USM6{promote_type(CT, R1, R2, S1, S2, S3)}(C, Rf1, Rf2, σ1, σ2, σ3)
end
(::Type{U6})(g::StaticVector) where {U6<:USM6} = U6(g[1], g[2], g[3], g[4], g[5], g[6])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::USM6) = SVector{6}(g.C, g.Rf1, g.Rf2, g.σ1, g.σ2, g.σ3)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
Base.one(::Type{U6}) where {U6<:USM6} = U6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::USM6, i::Int)
    if i == 1
        return p.C
    elseif i == 2
        return p.Rf1
    elseif i == 3
        return p.Rf2
    elseif i == 4
        return p.σ1
    elseif i == 5
        return p.σ2
    elseif i == 6
        return p.σ3
    else
        throw(BoundsError(r, i))
    end
end

export USMEM
"""
    USMEM{T} <: AstroCoord

Unified State Model Orbital Elements. 6D parameterziation of the orbit using Velocity Hodograph and Exponential Mapping
C - Velocity Hodograph Component Normal to the Radial Vector Laying in the Orbital Plane
Rf1 - Velocity Hodograph Component 90 degrees ahead of the Eccentricity Vector - Along the Intermediate Rotating Frame X-Axis
Rf2 - Velocity Hodograph Component 90 degrees ahead of the Eccentricity Vector - Along the Intermediate Rotating Frame Y-Axis
a1 - First Exponential Mapping Component
a2 - First Exponential Mapping Component
a3 - First Exponential Mapping Component

Constructors
USMEM(C, Rf1, Rf2, a1, a2, a3)
USMEM(X::AbstractArray)
USMEM(X::AstroCoord, μ::Number)

"""
struct USMEM{T} <: AstroCoord{6,T}
    C::T
    Rf1::T
    Rf2::T
    a1::T
    a2::T
    a3::T
    @inline USMEM{T}(C, Rf1, Rf2, a1, a2, a3) where {T} = new{T}(C, Rf1, Rf2, a1, a2, a3)
    @inline USMEM{T}(p::USMEM) where {T} = new{T}(p.C, p.Rf1, p.Rf2, p.a1, p.a2, p.a3)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
USMEM(X::AbstractArray{T,1}) where {T} = USMEM{T}(X...)
function USMEM(C::CT, Rf1::R1, Rf2::R2, a1::A1, a2::A2, a3::A3) where {CT,R1,R2,A1,A2,A3}
    return USMEM{promote_type(CT, R1, R2, A1, A2, A3)}(C, Rf1, Rf2, a1, a2, a3)
end
(::Type{UET})(g::StaticVector) where {UET<:USMEM} = UET(g[1], g[2], g[3], g[4], g[5], g[6])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::USMEM) = SVector{6}(g.C, g.Rf1, g.Rf2, g.a1, g.a2, g.a3)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
Base.one(::Type{UET}) where {UET<:USMEM} = UET(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::USMEM, i::Int)
    if i == 1
        return p.C
    elseif i == 2
        return p.Rf1
    elseif i == 3
        return p.Rf2
    elseif i == 4
        return p.a1
    elseif i == 5
        return p.a2
    elseif i == 6
        return p.a3
    else
        throw(BoundsError(r, i))
    end
end
