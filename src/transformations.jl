"""
    abstract type AstrodynamicsTransformation <: Transformation

An abstract type representing a Transformation Between Astrodynamics Coordinates
"""
abstract type AstrodynamicsTransformation <: Transformation end

"""
    abstract type AstroCoordTransformation <: AstrodynamicsTransformation

An abstract type representing a Transformation Between Astrodynamics Coordinates with no Time Regularization
"""
abstract type AstroCoordTransformation <: AstrodynamicsTransformation end

# ~~~~~~~~~~~~~~~ Cartesian <=> Keplerian ~~~~~~~~~~~~~~~ #
struct CartesiantoKeplerianTransform <: AstroCoordTransformation end

@inline function (::CartesiantoKeplerianTransform)(
    x::Cartesian{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Keplerian{RT}(cart2koe(params(x), μ))
end

const CartesiantoKeplerian = CartesiantoKeplerianTransform()

struct KepleriantoCartesianTransform <: AstroCoordTransformation end

@inline function (::KepleriantoCartesianTransform)(
    x::Keplerian{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Cartesian{RT}(koe2cart(params(x), μ))
end

const KepleriantoCartesian = KepleriantoCartesianTransform()

Base.inv(::CartesiantoKeplerianTransform) = KepleriantoCartesianTransform()
Base.inv(::KepleriantoCartesianTransform) = CartesiantoKeplerianTransform()

# ~~~~~~~~~~~~~~~ Keplerian <=> USM7 ~~~~~~~~~~~~~~~ #
struct KepleriantoUSM7Transform <: AstroCoordTransformation end

@inline function (::KepleriantoUSM7Transform)(
    x::Keplerian{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return USM7{RT}(koe2USM7(params(x), μ))
end

const KepleriantoUSM7 = KepleriantoUSM7Transform()

struct USM7toKeplerianTransform <: AstroCoordTransformation end

@inline function (::USM7toKeplerianTransform)(x::USM7{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Keplerian{RT}(USM72koe(params(x), μ))
end

const USM7toKeplerian = USM7toKeplerianTransform()

Base.inv(::KepleriantoUSM7Transform) = USM7toKeplerianTransform()
Base.inv(::USM7toKeplerianTransform) = KepleriantoUSM7Transform()

# ~~~~~~~~~~~~~~~ Cartesian <=> USM7 ~~~~~~~~~~~~~~~ #
const CartesiantoUSM7 = KepleriantoUSM7 ∘ CartesiantoKeplerian
const USM7toCartesian = KepleriantoCartesian ∘ USM7toKeplerian

# ~~~~~~~~~~~~~~~ USM6 <=> USM7 ~~~~~~~~~~~~~~~ #
struct USM6toUSM7Transform <: AstroCoordTransformation end

@inline function (::USM6toUSM7Transform)(x::USM6{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return USM7{RT}(USM62USM7(params(x), μ))
end

const USM6toUSM7 = USM6toUSM7Transform()

struct USM7toUSM6Transform <: AstroCoordTransformation end

@inline function (::USM7toUSM6Transform)(x::USM7{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return USM6{RT}(USM72USM6(params(x), μ))
end

const USM7toUSM6 = USM7toUSM6Transform()

Base.inv(::USM6toUSM7Transform) = USM7toUSM6Transform()
Base.inv(::USM7toUSM6Transform) = USM6toUSM7Transform()

# ~~~~~~~~~~~~~~~ Cartesian <=> USM6 ~~~~~~~~~~~~~~~ #
const CartesiantoUSM6 = USM7toUSM6 ∘ CartesiantoUSM7
const USM6toCartesian = USM7toCartesian ∘ USM6toUSM7

# ~~~~~~~~~~~~~~~ Keplerian <=> USM6 ~~~~~~~~~~~~~~~ #
const KepleriantoUSM6 = USM7toUSM6 ∘ KepleriantoUSM7
const USM6toKeplerian = USM7toKeplerian ∘ USM6toUSM7

# ~~~~~~~~~~~~~~~ USMEM <=> USM7 ~~~~~~~~~~~~~~~ #
struct USMEMtoUSM7Transform <: AstroCoordTransformation end

@inline function (::USMEMtoUSM7Transform)(x::USMEM{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return USM7{RT}(USMEM2USM7(params(x), μ))
end

const USMEMtoUSM7 = USMEMtoUSM7Transform()

struct USM7toUSMEMTransform <: AstroCoordTransformation end

@inline function (::USM7toUSMEMTransform)(x::USM7{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return USMEM{RT}(USM72USMEM(params(x), μ))
end

const USM7toUSMEM = USM7toUSMEMTransform()

Base.inv(::USMEMtoUSM7Transform) = USM7toUSMEMTransform()
Base.inv(::USM7toUSMEMTransform) = USMEMtoUSM7Transform()

# ~~~~~~~~~~~~~~~ Cartesian <=> USMEM ~~~~~~~~~~~~~~~ #
const CartesiantoUSMEM = USM7toUSMEM ∘ CartesiantoUSM7
const USMEMtoCartesian = USM7toCartesian ∘ USMEMtoUSM7

# ~~~~~~~~~~~~~~~ Keplerian <=> USMEM ~~~~~~~~~~~~~~~ #
const KepleriantoUSMEM = USM7toUSMEM ∘ KepleriantoUSM7
const USMEMtoKeplerian = USM7toKeplerian ∘ USMEMtoUSM7

# ~~~~~~~~~~~~~~~ USM6 <=> USMEM ~~~~~~~~~~~~~~~ #
const USM6toUSMEM = USM7toUSMEM ∘ USM6toUSM7
const USMEMtoUSM6 = USM7toUSM6 ∘ USMEMtoUSM7

# ~~~~~~~~~~~~~~~ Cartesian <=> Milankovich ~~~~~~~~~~~~~~~ #
struct CartesiantoMilankovichTransform <: AstroCoordTransformation end

@inline function (::CartesiantoMilankovichTransform)(
    x::Cartesian{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Milankovich{RT}(cart2Mil(params(x), μ))
end

const CartesiantoMilankovich = CartesiantoMilankovichTransform()

struct MilankovichtoCartesianTransform <: AstroCoordTransformation end

@inline function (::MilankovichtoCartesianTransform)(
    x::Milankovich{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Cartesian{RT}(Mil2cart(params(x), μ))
end

const MilankovichtoCartesian = MilankovichtoCartesianTransform()

Base.inv(::MilankovichtoCartesianTransform) = CartesiantoMilankovichTransform()
Base.inv(::CartesiantoMilankovichTransform) = MilankovichtoCartesianTransform()

# ~~~~~~~~~~~~~~~ Keplerian <=> Milankovich ~~~~~~~~~~~~~~~ #
const KepleriantoMilankovich = CartesiantoMilankovich ∘ KepleriantoCartesian
const MilankovichtoKeplerian = CartesiantoKeplerian ∘ MilankovichtoCartesian

# ~~~~~~~~~~~~~~~ USM7 <=> Milankovich ~~~~~~~~~~~~~~~ #
const USM7toMilankovich = CartesiantoMilankovich ∘ USM7toCartesian
const MilankovichtoUSM7 = CartesiantoUSM7 ∘ MilankovichtoCartesian

# ~~~~~~~~~~~~~~~ USM6 <=> Milankovich ~~~~~~~~~~~~~~~ #
const USM6toMilankovich = CartesiantoMilankovich ∘ USM6toCartesian
const MilankovichtoUSM6 = CartesiantoUSM6 ∘ MilankovichtoCartesian

# ~~~~~~~~~~~~~~~ USMEM <=> Milankovich ~~~~~~~~~~~~~~~ #
const USMEMtoMilankovich = CartesiantoMilankovich ∘ USMEMtoCartesian
const MilankovichtoUSMEM = CartesiantoUSMEM ∘ MilankovichtoCartesian

# ~~~~~~~~~~~~~~~ Keplerian to ModifiedEquinoctial ~~~~~~~~~~~~~~~ #
struct KepleriantoModifiedEquinoctialTransform <: AstroCoordTransformation end

@inline function (::KepleriantoModifiedEquinoctialTransform)(
    x::Keplerian{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return ModEq{RT}(koe2ModEq(params(x), μ))
end

const KepleriantoModifiedEquinoctial = KepleriantoModifiedEquinoctialTransform()

struct ModifiedEquinoctialtoKeplerianTransform <: AstroCoordTransformation end

@inline function (::ModifiedEquinoctialtoKeplerianTransform)(
    x::ModEq{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Keplerian{RT}(ModEq2koe(params(x), μ))
end

const ModifiedEquinoctialtoKeplerian = ModifiedEquinoctialtoKeplerianTransform()

function Base.inv(::KepleriantoModifiedEquinoctialTransform)
    return ModifiedEquinoctialtoKeplerianTransform()
end
function Base.inv(::ModifiedEquinoctialtoKeplerianTransform)
    return KepleriantoModifiedEquinoctialTransform()
end

# ~~~~~~~~~~~~~~~ Cartesian <=> Modified Equinoctial ~~~~~~~~~~~~~~~ #
const CartesiantoModifiedEquinoctial = KepleriantoModifiedEquinoctial ∘ CartesiantoKeplerian
const ModifiedEquinoctialtoCartesian = KepleriantoCartesian ∘ ModifiedEquinoctialtoKeplerian

# ~~~~~~~~~~~~~~~ USM7 <=> Modified Equinoctial ~~~~~~~~~~~~~~~ #
const USM7toModifiedEquinoctial = KepleriantoModifiedEquinoctial ∘ USM7toKeplerian
const ModifiedEquinoctialtoUSM7 = KepleriantoUSM7 ∘ ModifiedEquinoctialtoKeplerian

# ~~~~~~~~~~~~~~~ USM6 <=> Modified Equinoctial ~~~~~~~~~~~~~~~ #
const USM6toModifiedEquinoctial = KepleriantoModifiedEquinoctial ∘ USM6toKeplerian
const ModifiedEquinoctialtoUSM6 = KepleriantoUSM6 ∘ ModifiedEquinoctialtoKeplerian

# ~~~~~~~~~~~~~~~ USMEM <=> Modified Equinoctial ~~~~~~~~~~~~~~~ #
const USMEMtoModifiedEquinoctial = KepleriantoModifiedEquinoctial ∘ USMEMtoKeplerian
const ModifiedEquinoctialtoUSMEM = KepleriantoUSMEM ∘ ModifiedEquinoctialtoKeplerian

# ~~~~~~~~~~~~~~~ Milankovich <=> Modified Equinoctial ~~~~~~~~~~~~~~~ #
const MilankovichtoModifiedEquinoctial =
    KepleriantoModifiedEquinoctial ∘ MilankovichtoKeplerian
const ModifiedEquinoctialtoMilankovich =
    KepleriantoMilankovich ∘ ModifiedEquinoctialtoKeplerian

# ~~~~~~~~~~~~~~~ Cartesian to Cylindrical ~~~~~~~~~~~~~~~ #
struct CartesiantoCylindricalTransform <: AstroCoordTransformation end

@inline function (::CartesiantoCylindricalTransform)(
    x::Cartesian{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Cylindrical{RT}(cart2cylind(params(x), μ))
end

const CartesiantoCylindrical = CartesiantoCylindricalTransform()

struct CylindricaltoCartesianTransform <: AstroCoordTransformation end

@inline function (::CylindricaltoCartesianTransform)(
    x::Cylindrical{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Cartesian{RT}(cylind2cart(params(x), μ))
end

const CylindricaltoCartesian = CylindricaltoCartesianTransform()

Base.inv(::CylindricaltoCartesianTransform) = CartesiantoCylindricalTransform()
Base.inv(::CartesiantoCylindricalTransform) = CylindricaltoCartesianTransform()

# ~~~~~~~~~~~~~~~ Keplerian <=> Cylindrical ~~~~~~~~~~~~~~~ #
const KepleriantoCylindrical = CartesiantoCylindrical ∘ KepleriantoCartesian
const CylindricaltoKeplerian = CartesiantoKeplerian ∘ CylindricaltoCartesian

# ~~~~~~~~~~~~~~~ USM7 <=> Cylindrical ~~~~~~~~~~~~~~~ #
const USM7toCylindrical = CartesiantoCylindrical ∘ USM7toCartesian
const CylindricaltoUSM7 = CartesiantoUSM7 ∘ CylindricaltoCartesian

# ~~~~~~~~~~~~~~~ USM6 <=> Cylindrical ~~~~~~~~~~~~~~~ #
const USM6toCylindrical = CartesiantoCylindrical ∘ USM6toCartesian
const CylindricaltoUSM6 = CartesiantoUSM6 ∘ CylindricaltoCartesian

# ~~~~~~~~~~~~~~~ USMEM <=> Cylindrical ~~~~~~~~~~~~~~~ #
const USMEMtoCylindrical = CartesiantoCylindrical ∘ USMEMtoCartesian
const CylindricaltoUSMEM = CartesiantoUSMEM ∘ CylindricaltoCartesian

# ~~~~~~~~~~~~~~~ Milankovich <=> Cylindrical ~~~~~~~~~~~~~~~ #
const MilankovichtoCylindrical = CartesiantoCylindrical ∘ MilankovichtoCartesian
const CylindricaltoMilankovich = CartesiantoMilankovich ∘ CylindricaltoCartesian

# ~~~~~~~~~~~~~~~ Modified Equinoctial <=> Cylindrical ~~~~~~~~~~~~~~~ #
const ModifiedEquinoctialtoCylindrical =
    CartesiantoCylindrical ∘ ModifiedEquinoctialtoCartesian
const CylindricaltoModifiedEquinoctial =
    CartesiantoModifiedEquinoctial ∘ CylindricaltoCartesian

# ~~~~~~~~~~~~~~~ Cartesian to Spherical ~~~~~~~~~~~~~~~ #
struct CartesiantoSphericalTransform <: AstroCoordTransformation end

@inline function (::CartesiantoSphericalTransform)(
    x::Cartesian{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Spherical{RT}(cart2sphere(params(x), μ))
end

const CartesiantoSpherical = CartesiantoSphericalTransform()

struct SphericaltoCartesianTransform <: AstroCoordTransformation end

@inline function (::SphericaltoCartesianTransform)(
    x::Spherical{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Cartesian{RT}(sphere2cart(params(x), μ))
end

const SphericaltoCartesian = SphericaltoCartesianTransform()

Base.inv(::SphericaltoCartesianTransform) = CartesiantoSphericalTransform()
Base.inv(::CartesiantoSphericalTransform) = SphericaltoCartesianTransform()

# ~~~~~~~~~~~~~~~ Keplerian <=> Spherical ~~~~~~~~~~~~~~~ #
const KepleriantoSpherical = CartesiantoSpherical ∘ KepleriantoCartesian
const SphericaltoKeplerian = CartesiantoKeplerian ∘ SphericaltoCartesian

# ~~~~~~~~~~~~~~~ USM7 <=> Spherical ~~~~~~~~~~~~~~~ #
const USM7toSpherical = CartesiantoSpherical ∘ USM7toCartesian
const SphericaltoUSM7 = CartesiantoUSM7 ∘ SphericaltoCartesian

# ~~~~~~~~~~~~~~~ USM6 <=> Spherical ~~~~~~~~~~~~~~~ #
const USM6toSpherical = CartesiantoSpherical ∘ USM6toCartesian
const SphericaltoUSM6 = CartesiantoUSM6 ∘ SphericaltoCartesian

# ~~~~~~~~~~~~~~~ USMEM <=> Spherical ~~~~~~~~~~~~~~~ #
const USMEMtoSpherical = CartesiantoSpherical ∘ USMEMtoCartesian
const SphericaltoUSMEM = CartesiantoUSMEM ∘ SphericaltoCartesian

# ~~~~~~~~~~~~~~~ Milankovich <=> Spherical ~~~~~~~~~~~~~~~ #
const MilankovichtoSpherical = CartesiantoSpherical ∘ MilankovichtoCartesian
const SphericaltoMilankovich = CartesiantoMilankovich ∘ SphericaltoCartesian

# ~~~~~~~~~~~~~~~ Modified Equinoctial <=> Spherical ~~~~~~~~~~~~~~~ #
const ModifiedEquinoctialtoSpherical = CartesiantoSpherical ∘ ModifiedEquinoctialtoCartesian
const SphericaltoModifiedEquinoctial = CartesiantoModifiedEquinoctial ∘ SphericaltoCartesian

# ~~~~~~~~~~~~~~~ Cylindrical <=> Spherical ~~~~~~~~~~~~~~~ #
const CylindricaltoSpherical = CartesiantoSpherical ∘ CylindricaltoCartesian
const SphericaltoCylindrical = CartesiantoCylindrical ∘ SphericaltoCartesian

# ~~~~~~~~~~~~~~~ Cartesian to Delaunay ~~~~~~~~~~~~~~~ #
struct CartesiantoDelaunayTransform <: AstroCoordTransformation end

@inline function (::CartesiantoDelaunayTransform)(
    x::Cartesian{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Delaunay{RT}(cart2delaunay(params(x), μ))
end

const CartesiantoDelaunay = CartesiantoDelaunayTransform()

struct DelaunaytoCartesianTransform <: AstroCoordTransformation end

@inline function (::DelaunaytoCartesianTransform)(
    x::Delaunay{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Cartesian{RT}(delaunay2cart(params(x), μ))
end

const DelaunaytoCartesian = DelaunaytoCartesianTransform()

Base.inv(::DelaunaytoCartesianTransform) = CartesiantoDelaunayTransform()
Base.inv(::CartesiantoDelaunayTransform) = DelaunaytoCartesianTransform()

# ~~~~~~~~~~~~~~~ Keplerian <=> Delaunay ~~~~~~~~~~~~~~~ #
const KepleriantoDelaunay = CartesiantoDelaunay ∘ KepleriantoCartesian
const DelaunaytoKeplerian = CartesiantoKeplerian ∘ DelaunaytoCartesian

# ~~~~~~~~~~~~~~~ USM7 <=> Delaunay ~~~~~~~~~~~~~~~ #
const USM7toDelaunay = CartesiantoDelaunay ∘ USM7toCartesian
const DelaunaytoUSM7 = CartesiantoUSM7 ∘ DelaunaytoCartesian

# ~~~~~~~~~~~~~~~ USM6 <=> Delaunay ~~~~~~~~~~~~~~~ #
const USM6toDelaunay = CartesiantoDelaunay ∘ USM6toCartesian
const DelaunaytoUSM6 = CartesiantoUSM6 ∘ DelaunaytoCartesian

# ~~~~~~~~~~~~~~~ USMEM <=> Delaunay ~~~~~~~~~~~~~~~ #
const USMEMtoDelaunay = CartesiantoDelaunay ∘ USMEMtoCartesian
const DelaunaytoUSMEM = CartesiantoUSMEM ∘ DelaunaytoCartesian

# ~~~~~~~~~~~~~~~ Milankovich <=> Delaunay ~~~~~~~~~~~~~~~ #
const MilankovichtoDelaunay = CartesiantoDelaunay ∘ MilankovichtoCartesian
const DelaunaytoMilankovich = CartesiantoMilankovich ∘ DelaunaytoCartesian

# ~~~~~~~~~~~~~~~ Modified Equinoctial <=> Delaunay ~~~~~~~~~~~~~~~ #
const ModifiedEquinoctialtoDelaunay = CartesiantoDelaunay ∘ ModifiedEquinoctialtoCartesian
const DelaunaytoModifiedEquinoctial = CartesiantoModifiedEquinoctial ∘ DelaunaytoCartesian

# ~~~~~~~~~~~~~~~ Cylindrical <=> Delaunay ~~~~~~~~~~~~~~~ #
const CylindricaltoDelaunay = CartesiantoDelaunay ∘ CylindricaltoCartesian
const DelaunaytoCylindrical = CartesiantoCylindrical ∘ DelaunaytoCartesian

# ~~~~~~~~~~~~~~~ Spherical <=> Delaunay ~~~~~~~~~~~~~~~ #
const SphericaltoDelaunay = CartesiantoDelaunay ∘ SphericaltoCartesian
const DelaunaytoSpherical = CartesiantoSpherical ∘ DelaunaytoCartesian

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~ Additional Constructors ~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#TODO: NOT MAINTAINABLE
#TODO: CAN I DO THIS WITH METAPROGRAMMING?
# ~~~~~~~~~~~~~~~ Cartesian ~~~~~~~~~~~~~~~ #
Cartesian(X::Cartesian{T}, μ::Number) where {T<:Number} = X
Cartesian(X::Keplerian{T}, μ::Number) where {T<:Number} = KepleriantoCartesian(X, μ)
Cartesian(X::USM7{T}, μ::Number) where {T<:Number} = USM7toCartesian(X, μ)
Cartesian(X::USM6{T}, μ::Number) where {T<:Number} = USM6toCartesian(X, μ)
Cartesian(X::USMEM{T}, μ::Number) where {T<:Number} = USMEMtoCartesian(X, μ)
Cartesian(X::Milankovich{T}, μ::Number) where {T<:Number} = MilankovichtoCartesian(X, μ)
Cartesian(X::ModEq{T}, μ::Number) where {T<:Number} = ModifiedEquinoctialtoCartesian(X, μ)
Cartesian(X::Cylindrical{T}, μ::Number) where {T<:Number} = CylindricaltoCartesian(X, μ)
Cartesian(X::Spherical{T}, μ::Number) where {T<:Number} = SphericaltoCartesian(X, μ)
Cartesian(X::Delaunay{T}, μ::Number) where {T<:Number} = DelaunaytoCartesian(X, μ)

# ~~~~~~~~~~~~~~~ Keplerian ~~~~~~~~~~~~~~~ #
Keplerian(X::Keplerian{T}, μ::Number) where {T<:Number} = X
Keplerian(X::Cartesian{T}, μ::Number) where {T<:Number} = CartesiantoKeplerian(X, μ)
Keplerian(X::USM7{T}, μ::Number) where {T<:Number} = USM7toKeplerian(X, μ)
Keplerian(X::USM6{T}, μ::Number) where {T<:Number} = USM6toKeplerian(X, μ)
Keplerian(X::USMEM{T}, μ::Number) where {T<:Number} = USMEMtoKeplerian(X, μ)
Keplerian(X::Milankovich{T}, μ::Number) where {T<:Number} = MilankovichtoKeplerian(X, μ)
Keplerian(X::ModEq{T}, μ::Number) where {T<:Number} = ModifiedEquinoctialtoKeplerian(X, μ)
Keplerian(X::Cylindrical{T}, μ::Number) where {T<:Number} = CylindricaltoKeplerian(X, μ)
Keplerian(X::Spherical{T}, μ::Number) where {T<:Number} = SphericaltoKeplerian(X, μ)
Keplerian(X::Delaunay{T}, μ::Number) where {T<:Number} = DelaunaytoKeplerian(X, μ)

# ~~~~~~~~~~~~~~~ USM7 ~~~~~~~~~~~~~~~ #
USM7(X::USM7{T}, μ::Number) where {T<:Number} = X
USM7(X::Cartesian{T}, μ::Number) where {T<:Number} = CartesiantoUSM7(X, μ)
USM7(X::Keplerian{T}, μ::Number) where {T<:Number} = KepleriantoUSM7(X, μ)
USM7(X::USM6{T}, μ::Number) where {T<:Number} = USM6toUSM7(X, μ)
USM7(X::USMEM{T}, μ::Number) where {T<:Number} = USMEMtoUSM7(X, μ)
USM7(X::Milankovich{T}, μ::Number) where {T<:Number} = MilankovichtoUSM7(X, μ)
USM7(X::ModEq{T}, μ::Number) where {T<:Number} = ModifiedEquinoctialtoUSM7(X, μ)
USM7(X::Cylindrical{T}, μ::Number) where {T<:Number} = CylindricaltoUSM7(X, μ)
USM7(X::Spherical{T}, μ::Number) where {T<:Number} = SphericaltoUSM7(X, μ)
USM7(X::Delaunay{T}, μ::Number) where {T<:Number} = DelaunaytoUSM7(X, μ)

# ~~~~~~~~~~~~~~~ USM6 ~~~~~~~~~~~~~~~ #
USM6(X::USM6{T}, μ::Number) where {T<:Number} = X
USM6(X::Cartesian{T}, μ::Number) where {T<:Number} = CartesiantoUSM6(X, μ)
USM6(X::Keplerian{T}, μ::Number) where {T<:Number} = KepleriantoUSM6(X, μ)
USM6(X::USM7{T}, μ::Number) where {T<:Number} = USM7toUSM6(X, μ)
USM6(X::USMEM{T}, μ::Number) where {T<:Number} = USMEMtoUSM6(X, μ)
USM6(X::Milankovich{T}, μ::Number) where {T<:Number} = MilankovichtoUSM6(X, μ)
USM6(X::ModEq{T}, μ::Number) where {T<:Number} = ModifiedEquinoctialtoUSM6(X, μ)
USM6(X::Cylindrical{T}, μ::Number) where {T<:Number} = CylindricaltoUSM6(X, μ)
USM6(X::Spherical{T}, μ::Number) where {T<:Number} = SphericaltoUSM6(X, μ)
USM6(X::Delaunay{T}, μ::Number) where {T<:Number} = DelaunaytoUSM6(X, μ)

# ~~~~~~~~~~~~~~~ USMEM ~~~~~~~~~~~~~~~ #
USMEM(X::USMEM{T}, μ::Number) where {T<:Number} = X
USMEM(X::Cartesian{T}, μ::Number) where {T<:Number} = CartesiantoUSMEM(X, μ)
USMEM(X::Keplerian{T}, μ::Number) where {T<:Number} = KepleriantoUSMEM(X, μ)
USMEM(X::USM7{T}, μ::Number) where {T<:Number} = USM7toUSMEM(X, μ)
USMEM(X::USM6{T}, μ::Number) where {T<:Number} = USM6toUSMEM(X, μ)
USMEM(X::Milankovich{T}, μ::Number) where {T<:Number} = MilankovichtoUSMEM(X, μ)
USMEM(X::ModEq{T}, μ::Number) where {T<:Number} = ModifiedEquinoctialtoUSMEM(X, μ)
USMEM(X::Cylindrical{T}, μ::Number) where {T<:Number} = CylindricaltoUSMEM(X, μ)
USMEM(X::Spherical{T}, μ::Number) where {T<:Number} = SphericaltoUSMEM(X, μ)
USMEM(X::Delaunay{T}, μ::Number) where {T<:Number} = DelaunaytoUSMEM(X, μ)

# ~~~~~~~~~~~~~~~ Milankovich ~~~~~~~~~~~~~~~ #
Milankovich(X::Milankovich{T}, μ::Number) where {T<:Number} = X
Milankovich(X::Cartesian{T}, μ::Number) where {T<:Number} = CartesiantoMilankovich(X, μ)
Milankovich(X::Keplerian{T}, μ::Number) where {T<:Number} = KepleriantoMilankovich(X, μ)
Milankovich(X::USM7{T}, μ::Number) where {T<:Number} = USM7toMilankovich(X, μ)
Milankovich(X::USM6{T}, μ::Number) where {T<:Number} = USM6toMilankovich(X, μ)
Milankovich(X::USMEM{T}, μ::Number) where {T<:Number} = USMEMtoMilankovich(X, μ)
function Milankovich(X::ModEq{T}, μ::Number) where {T<:Number}
    return ModifiedEquinoctialtoMilankovich(X, μ)
end
Milankovich(X::Cylindrical{T}, μ::Number) where {T<:Number} = CylindricaltoMilankovich(X, μ)
Milankovich(X::Spherical{T}, μ::Number) where {T<:Number} = SphericaltoMilankovich(X, μ)
Milankovich(X::Delaunay{T}, μ::Number) where {T<:Number} = DelaunaytoMilankovich(X, μ)

# ~~~~~~~~~~~~~~~ ModEq ~~~~~~~~~~~~~~~ #
ModEq(X::ModEq{T}, μ::Number) where {T<:Number} = X
ModEq(X::Cartesian{T}, μ::Number) where {T<:Number} = CartesiantoModifiedEquinoctial(X, μ)
ModEq(X::Keplerian{T}, μ::Number) where {T<:Number} = KepleriantoModifiedEquinoctial(X, μ)
ModEq(X::USM7{T}, μ::Number) where {T<:Number} = USM7toModifiedEquinoctial(X, μ)
ModEq(X::USM6{T}, μ::Number) where {T<:Number} = USM6toModifiedEquinoctial(X, μ)
ModEq(X::USMEM{T}, μ::Number) where {T<:Number} = USMEMtoModifiedEquinoctial(X, μ)
function ModEq(X::Milankovich{T}, μ::Number) where {T<:Number}
    return MilankovichtoModifiedEquinoctial(X, μ)
end
function ModEq(X::Cylindrical{T}, μ::Number) where {T<:Number}
    return CylindricaltoModifiedEquinoctial(X, μ)
end
ModEq(X::Spherical{T}, μ::Number) where {T<:Number} = SphericaltoModifiedEquinoctial(X, μ)
ModEq(X::Delaunay{T}, μ::Number) where {T<:Number} = DelaunaytoModifiedEquinoctial(X, μ)

# ~~~~~~~~~~~~~~~ Cylindrical ~~~~~~~~~~~~~~~ #
Cylindrical(X::Cylindrical{T}, μ::Number) where {T<:Number} = X
Cylindrical(X::Cartesian{T}, μ::Number) where {T<:Number} = CartesiantoCylindrical(X, μ)
Cylindrical(X::Keplerian{T}, μ::Number) where {T<:Number} = KepleriantoCylindrical(X, μ)
Cylindrical(X::USM7{T}, μ::Number) where {T<:Number} = USM7toCylindrical(X, μ)
Cylindrical(X::USM6{T}, μ::Number) where {T<:Number} = USM6toCylindrical(X, μ)
Cylindrical(X::USMEM{T}, μ::Number) where {T<:Number} = USMEMtoCylindrical(X, μ)
Cylindrical(X::Milankovich{T}, μ::Number) where {T<:Number} = MilankovichtoCylindrical(X, μ)
function Cylindrical(X::ModEq{T}, μ::Number) where {T<:Number}
    return ModifiedEquinoctialtoCylindrical(X, μ)
end
Cylindrical(X::Spherical{T}, μ::Number) where {T<:Number} = SphericaltoCylindrical(X, μ)
Cylindrical(X::Delaunay{T}, μ::Number) where {T<:Number} = DelaunaytoCylindrical(X, μ)

# ~~~~~~~~~~~~~~~ Spherical ~~~~~~~~~~~~~~~ #
Spherical(X::Spherical{T}, μ::Number) where {T<:Number} = X
Spherical(X::Cartesian{T}, μ::Number) where {T<:Number} = CartesiantoSpherical(X, μ)
Spherical(X::Keplerian{T}, μ::Number) where {T<:Number} = KepleriantoSpherical(X, μ)
Spherical(X::USM7{T}, μ::Number) where {T<:Number} = USM7toSpherical(X, μ)
Spherical(X::USM6{T}, μ::Number) where {T<:Number} = USM6toSpherical(X, μ)
Spherical(X::USMEM{T}, μ::Number) where {T<:Number} = USMEMtoSpherical(X, μ)
Spherical(X::Milankovich{T}, μ::Number) where {T<:Number} = MilankovichtoSpherical(X, μ)
Spherical(X::ModEq{T}, μ::Number) where {T<:Number} = ModifiedEquinoctialtoSpherical(X, μ)
Spherical(X::Cylindrical{T}, μ::Number) where {T<:Number} = CylindricaltoSpherical(X, μ)
Spherical(X::Delaunay{T}, μ::Number) where {T<:Number} = DelaunaytoSpherical(X, μ)

# ~~~~~~~~~~~~~~~ Delaunay ~~~~~~~~~~~~~~~ #
Delaunay(X::Delaunay{T}, μ::Number) where {T<:Number} = X
Delaunay(X::Cartesian{T}, μ::Number) where {T<:Number} = CartesiantoDelaunay(X, μ)
Delaunay(X::Keplerian{T}, μ::Number) where {T<:Number} = KepleriantoDelaunay(X, μ)
Delaunay(X::USM7{T}, μ::Number) where {T<:Number} = USM7toDelaunay(X, μ)
Delaunay(X::USM6{T}, μ::Number) where {T<:Number} = USM6toDelaunay(X, μ)
Delaunay(X::USMEM{T}, μ::Number) where {T<:Number} = USMEMtoDelaunay(X, μ)
Delaunay(X::Milankovich{T}, μ::Number) where {T<:Number} = MilankovichtoDelaunay(X, μ)
Delaunay(X::ModEq{T}, μ::Number) where {T<:Number} = ModifiedEquinoctialtoDelaunay(X, μ)
Delaunay(X::Cylindrical{T}, μ::Number) where {T<:Number} = CylindricaltoDelaunay(X, μ)
Delaunay(X::Spherical{T}, μ::Number) where {T<:Number} = SphericaltoDelaunay(X, μ)
