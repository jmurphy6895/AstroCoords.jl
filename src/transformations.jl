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

@inline function (::CartesiantoKeplerianTransform)(x::Cartesian, μ::Number)
    return Keplerian(cart2koe(params(x), μ))
end

CartesiantoKeplerian = CartesiantoKeplerianTransform()

struct KepleriantoCartesianTransform <: AstroCoordTransformation end

@inline function (::KepleriantoCartesianTransform)(x::Keplerian, μ::Number)
    return Cartesian(koe2cart(params(x), μ))
end

KepleriantoCartesian = KepleriantoCartesianTransform()

Base.inv(::CartesiantoKeplerianTransform) = KepleriantoCartesianTransform()
Base.inv(::KepleriantoCartesianTransform) = CartesiantoKeplerianTransform()

# ~~~~~~~~~~~~~~~ Keplerian <=> USM7 ~~~~~~~~~~~~~~~ #
struct KepleriantoUSM7Transform <: AstroCoordTransformation end

@inline function (::KepleriantoUSM7Transform)(x::Keplerian, μ::Number)
    return USM7(koe2USM7(params(x), μ))
end

KepleriantoUSM7 = KepleriantoUSM7Transform()

struct USM7toKeplerianTransform <: AstroCoordTransformation end

@inline function (::USM7toKeplerianTransform)(x::USM7, μ::Number)
    return Keplerian(USM72koe(params(x), μ))
end

USM7toKeplerian = USM7toKeplerianTransform()

Base.inv(::KepleriantoUSM7Transform) = USM7toKeplerianTransform()
Base.inv(::USM7toKeplerianTransform) = KepleriantoUSM7Transform()

# ~~~~~~~~~~~~~~~ Cartesian <=> USM7 ~~~~~~~~~~~~~~~ #
CartesiantoUSM7 = KepleriantoUSM7 ∘ CartesiantoKeplerian
USM7toCartesian = KepleriantoCartesian ∘ USM7toKeplerian

# ~~~~~~~~~~~~~~~ USM6 <=> USM7 ~~~~~~~~~~~~~~~ #
struct USM6toUSM7Transform <: AstroCoordTransformation end

@inline function (::USM6toUSM7Transform)(x::USM6, μ::Number)
    return USM7(USM62USM7(params(x), μ))
end

USM6toUSM7 = USM6toUSM7Transform()

struct USM7toUSM6Transform <: AstroCoordTransformation end

@inline function (::USM7toUSM6Transform)(x::USM7, μ::Number)
    return USM6(USM72USM6(params(x), μ))
end

USM7toUSM6 = USM7toUSM6Transform()

Base.inv(::USM6toUSM7Transform) = USM7toUSM6Transform()
Base.inv(::USM7toUSM6Transform) = USM6toUSM7Transform()

# ~~~~~~~~~~~~~~~ Cartesian <=> USM6 ~~~~~~~~~~~~~~~ #
CartesiantoUSM6 = USM7toUSM6 ∘ CartesiantoUSM7
USM6toCartesian = USM7toCartesian ∘ USM6toUSM7

# ~~~~~~~~~~~~~~~ Keplerian <=> USM6 ~~~~~~~~~~~~~~~ #
KepleriantoUSM6 = USM7toUSM6 ∘ KepleriantoUSM7
USM6toKeplerian = USM7toKeplerian ∘ USM6toUSM7

# ~~~~~~~~~~~~~~~ USMEM <=> USM7 ~~~~~~~~~~~~~~~ #
struct USMEMtoUSM7Transform <: AstroCoordTransformation end

@inline function (::USMEMtoUSM7Transform)(x::USMEM, μ::Number)
    return USM7(USMEM2USM7(params(x), μ))
end

USMEMtoUSM7 = USMEMtoUSM7Transform()

struct USM7toUSMEMTransform <: AstroCoordTransformation end

@inline function (::USM7toUSMEMTransform)(x::USM7, μ::Number)
    return USMEM(USM72USMEM(params(x), μ))
end

USM7toUSMEM = USM7toUSMEMTransform()

Base.inv(::USMEMtoUSM7Transform) = USM7toUSMEMTransform()
Base.inv(::USM7toUSMEMTransform) = USMEMtoUSM7Transform()

# ~~~~~~~~~~~~~~~ Cartesian <=> USMEM ~~~~~~~~~~~~~~~ #
CartesiantoUSMEM = USM7toUSMEM ∘ CartesiantoUSM7
USMEMtoCartesian = USM7toCartesian ∘ USMEMtoUSM7

# ~~~~~~~~~~~~~~~ Keplerian <=> USMEM ~~~~~~~~~~~~~~~ #
KepleriantoUSMEM = USM7toUSMEM ∘ KepleriantoUSM7
USMEMtoKeplerian = USM7toKeplerian ∘ USMEMtoUSM7

# ~~~~~~~~~~~~~~~ USM6 <=> USMEM ~~~~~~~~~~~~~~~ #
USM6toUSMEM = USM7toUSMEM ∘ USM6toUSM7
USMEMtoUSM6 = USM7toUSM6 ∘ USMEMtoUSM7

# ~~~~~~~~~~~~~~~ Cartesian <=> Milankovich ~~~~~~~~~~~~~~~ #
struct CartesiantoMilankovichTransform <: AstroCoordTransformation end

@inline function (::CartesiantoMilankovichTransform)(x::Cartesian, μ::Number)
    return Milankovich(cart2Mil(params(x), μ))
end

CartesiantoMilankovich = CartesiantoMilankovichTransform()

struct MilankovichtoCartesianTransform <: AstroCoordTransformation end

@inline function (::MilankovichtoCartesianTransform)(x::Milankovich, μ::Number)
    return Cartesian(Mil2cart(params(x), μ))
end

MilankovichtoCartesian = MilankovichtoCartesianTransform()

Base.inv(::MilankovichtoCartesianTransform) = CartesiantoMilankovichTransform()
Base.inv(::CartesiantoMilankovichTransform) = MilankovichtoCartesianTransform()

# ~~~~~~~~~~~~~~~ Keplerian <=> Milankovich ~~~~~~~~~~~~~~~ #
KepleriantoMilankovich = CartesiantoMilankovich ∘ KepleriantoCartesian
MilankovichtoKeplerian = CartesiantoKeplerian ∘ MilankovichtoCartesian

# ~~~~~~~~~~~~~~~ USM7 <=> Milankovich ~~~~~~~~~~~~~~~ #
USM7toMilankovich = CartesiantoMilankovich ∘ USM7toCartesian
MilankovichtoUSM7 = CartesiantoUSM7 ∘ MilankovichtoCartesian

# ~~~~~~~~~~~~~~~ USM6 <=> Milankovich ~~~~~~~~~~~~~~~ #
USM6toMilankovich = CartesiantoMilankovich ∘ USM6toCartesian
MilankovichtoUSM6 = CartesiantoUSM6 ∘ MilankovichtoCartesian

# ~~~~~~~~~~~~~~~ USMEM <=> Milankovich ~~~~~~~~~~~~~~~ #
USMEMtoMilankovich = CartesiantoMilankovich ∘ USMEMtoCartesian
MilankovichtoUSMEM = CartesiantoUSMEM ∘ MilankovichtoCartesian

# ~~~~~~~~~~~~~~~ Keplerian to ModifiedEquinoctial ~~~~~~~~~~~~~~~ #
struct KepleriantoModifiedEquinoctialTransform <: AstroCoordTransformation end

@inline function (::KepleriantoModifiedEquinoctialTransform)(x::Keplerian, μ::Number)
    return ModEq(koe2ModEq(params(x), μ))
end

KepleriantoModifiedEquinoctial = KepleriantoModifiedEquinoctialTransform()

struct ModifiedEquinoctialtoKeplerianTransform <: AstroCoordTransformation end

@inline function (::ModifiedEquinoctialtoKeplerianTransform)(x::ModEq, μ::Number)
    return Keplerian(ModEq2koe(params(x), μ))
end

ModifiedEquinoctialtoKeplerian = ModifiedEquinoctialtoKeplerianTransform()

function Base.inv(::KepleriantoModifiedEquinoctialTransform)
    return ModifiedEquinoctialtoKeplerianTransform()
end
function Base.inv(::ModifiedEquinoctialtoKeplerianTransform)
    return KepleriantoModifiedEquinoctialTransform()
end

# ~~~~~~~~~~~~~~~ Cartesian <=> Modified Equinoctial ~~~~~~~~~~~~~~~ #
CartesiantoModifiedEquinoctial = KepleriantoModifiedEquinoctial ∘ CartesiantoKeplerian
ModifiedEquinoctialtoCartesian = KepleriantoCartesian ∘ ModifiedEquinoctialtoKeplerian

# ~~~~~~~~~~~~~~~ USM7 <=> Modified Equinoctial ~~~~~~~~~~~~~~~ #
USM7toModifiedEquinoctial = KepleriantoModifiedEquinoctial ∘ USM7toKeplerian
ModifiedEquinoctialtoUSM7 = KepleriantoUSM7 ∘ ModifiedEquinoctialtoKeplerian

# ~~~~~~~~~~~~~~~ USM6 <=> Modified Equinoctial ~~~~~~~~~~~~~~~ #
USM6toModifiedEquinoctial = KepleriantoModifiedEquinoctial ∘ USM6toKeplerian
ModifiedEquinoctialtoUSM6 = KepleriantoUSM6 ∘ ModifiedEquinoctialtoKeplerian

# ~~~~~~~~~~~~~~~ USMEM <=> Modified Equinoctial ~~~~~~~~~~~~~~~ #
USMEMtoModifiedEquinoctial = KepleriantoModifiedEquinoctial ∘ USMEMtoKeplerian
ModifiedEquinoctialtoUSMEM = KepleriantoUSMEM ∘ ModifiedEquinoctialtoKeplerian

# ~~~~~~~~~~~~~~~ Milankovich <=> Modified Equinoctial ~~~~~~~~~~~~~~~ #
MilankovichtoModifiedEquinoctial = KepleriantoModifiedEquinoctial ∘ MilankovichtoKeplerian
ModifiedEquinoctialtoMilankovich = KepleriantoMilankovich ∘ ModifiedEquinoctialtoKeplerian

# ~~~~~~~~~~~~~~~ Cartesian to Cylindrical ~~~~~~~~~~~~~~~ #
struct CartesiantoCylindricalTransform <: AstroCoordTransformation end

@inline function (::CartesiantoCylindricalTransform)(x::Cartesian, μ::Number)
    return Cylindrical(cart2cylind(params(x), μ))
end

CartesiantoCylindrical = CartesiantoCylindricalTransform()

struct CylindricaltoCartesianTransform <: AstroCoordTransformation end

@inline function (::CylindricaltoCartesianTransform)(x::Cylindrical, μ::Number)
    return Cartesian(cylind2cart(params(x), μ))
end

CylindricaltoCartesian = CylindricaltoCartesianTransform()

Base.inv(::CylindricaltoCartesianTransform) = CartesiantoCylindricalTransform()
Base.inv(::CartesiantoCylindricalTransform) = CylindricaltoCartesianTransform()

# ~~~~~~~~~~~~~~~ Keplerian <=> Cylindrical ~~~~~~~~~~~~~~~ #
KepleriantoCylindrical = CartesiantoCylindrical ∘ KepleriantoCartesian
CylindricaltoKeplerian = CartesiantoKeplerian ∘ CylindricaltoCartesian

# ~~~~~~~~~~~~~~~ USM7 <=> Cylindrical ~~~~~~~~~~~~~~~ #
USM7toCylindrical = CartesiantoCylindrical ∘ USM7toCartesian
CylindricaltoUSM7 = CartesiantoUSM7 ∘ CylindricaltoCartesian

# ~~~~~~~~~~~~~~~ USM6 <=> Cylindrical ~~~~~~~~~~~~~~~ #
USM6toCylindrical = CartesiantoCylindrical ∘ USM6toCartesian
CylindricaltoUSM6 = CartesiantoUSM6 ∘ CylindricaltoCartesian

# ~~~~~~~~~~~~~~~ USMEM <=> Cylindrical ~~~~~~~~~~~~~~~ #
USMEMtoCylindrical = CartesiantoCylindrical ∘ USMEMtoCartesian
CylindricaltoUSMEM = CartesiantoUSMEM ∘ CylindricaltoCartesian

# ~~~~~~~~~~~~~~~ Milankovich <=> Cylindrical ~~~~~~~~~~~~~~~ #
MilankovichtoCylindrical = CartesiantoCylindrical ∘ MilankovichtoCartesian
CylindricaltoMilankovich = CartesiantoMilankovich ∘ CylindricaltoCartesian

# ~~~~~~~~~~~~~~~ Modified Equinoctial <=> Cylindrical ~~~~~~~~~~~~~~~ #
ModifiedEquinoctialtoCylindrical = CartesiantoCylindrical ∘ ModifiedEquinoctialtoCartesian
CylindricaltoModifiedEquinoctial = CartesiantoModifiedEquinoctial ∘ CylindricaltoCartesian

# ~~~~~~~~~~~~~~~ Cartesian to Spherical ~~~~~~~~~~~~~~~ #
struct CartesiantoSphericalTransform <: AstroCoordTransformation end

@inline function (::CartesiantoSphericalTransform)(x::Cartesian, μ::Number)
    return Spherical(cart2sphere(params(x), μ))
end

CartesiantoSpherical = CartesiantoSphericalTransform()

struct SphericaltoCartesianTransform <: AstroCoordTransformation end

@inline function (::SphericaltoCartesianTransform)(x::Spherical, μ::Number)
    return Cartesian(sphere2cart(params(x), μ))
end

SphericaltoCartesian = SphericaltoCartesianTransform()

Base.inv(::SphericaltoCartesianTransform) = CartesiantoSphericalTransform()
Base.inv(::CartesiantoSphericalTransform) = SphericaltoCartesianTransform()

# ~~~~~~~~~~~~~~~ Keplerian <=> Spherical ~~~~~~~~~~~~~~~ #
KepleriantoSpherical = CartesiantoSpherical ∘ KepleriantoCartesian
SphericaltoKeplerian = CartesiantoKeplerian ∘ SphericaltoCartesian

# ~~~~~~~~~~~~~~~ USM7 <=> Spherical ~~~~~~~~~~~~~~~ #
USM7toSpherical = CartesiantoSpherical ∘ USM7toCartesian
SphericaltoUSM7 = CartesiantoUSM7 ∘ SphericaltoCartesian

# ~~~~~~~~~~~~~~~ USM6 <=> Spherical ~~~~~~~~~~~~~~~ #
USM6toSpherical = CartesiantoSpherical ∘ USM6toCartesian
SphericaltoUSM6 = CartesiantoUSM6 ∘ SphericaltoCartesian

# ~~~~~~~~~~~~~~~ USMEM <=> Spherical ~~~~~~~~~~~~~~~ #
USMEMtoSpherical = CartesiantoSpherical ∘ USMEMtoCartesian
SphericaltoUSMEM = CartesiantoUSMEM ∘ SphericaltoCartesian

# ~~~~~~~~~~~~~~~ Milankovich <=> Spherical ~~~~~~~~~~~~~~~ #
MilankovichtoSpherical = CartesiantoSpherical ∘ MilankovichtoCartesian
SphericaltoMilankovich = CartesiantoMilankovich ∘ SphericaltoCartesian

# ~~~~~~~~~~~~~~~ Modified Equinoctial <=> Spherical ~~~~~~~~~~~~~~~ #
ModifiedEquinoctialtoSpherical = CartesiantoSpherical ∘ ModifiedEquinoctialtoCartesian
SphericaltoModifiedEquinoctial = CartesiantoModifiedEquinoctial ∘ SphericaltoCartesian

# ~~~~~~~~~~~~~~~ Cylindrical <=> Spherical ~~~~~~~~~~~~~~~ #
CylindricaltoSpherical = CartesiantoSpherical ∘ CylindricaltoCartesian
SphericaltoCylindrical = CartesiantoCylindrical ∘ SphericaltoCartesian

# ~~~~~~~~~~~~~~~ Cartesian to Delaunay ~~~~~~~~~~~~~~~ #
struct CartesiantoDelaunayTransform <: AstroCoordTransformation end

@inline function (::CartesiantoDelaunayTransform)(x::Cartesian, μ::Number)
    return Delaunay(cart2delaunay(params(x), μ))
end

CartesiantoDelaunay = CartesiantoDelaunayTransform()

struct DelaunaytoCartesianTransform <: AstroCoordTransformation end

@inline function (::DelaunaytoCartesianTransform)(x::Delaunay, μ::Number)
    return Cartesian(delaunay2cart(params(x), μ))
end

DelaunaytoCartesian = DelaunaytoCartesianTransform()

Base.inv(::DelaunaytoCartesianTransform) = CartesiantoDelaunayTransform()
Base.inv(::CartesiantoDelaunayTransform) = DelaunaytoCartesianTransform()

# ~~~~~~~~~~~~~~~ Keplerian <=> Delaunay ~~~~~~~~~~~~~~~ #
KepleriantoDelaunay = CartesiantoDelaunay ∘ KepleriantoCartesian
DelaunaytoKeplerian = CartesiantoKeplerian ∘ DelaunaytoCartesian

# ~~~~~~~~~~~~~~~ USM7 <=> Delaunay ~~~~~~~~~~~~~~~ #
USM7toDelaunay = CartesiantoDelaunay ∘ USM7toCartesian
DelaunaytoUSM7 = CartesiantoUSM7 ∘ DelaunaytoCartesian

# ~~~~~~~~~~~~~~~ USM6 <=> Delaunay ~~~~~~~~~~~~~~~ #
USM6toDelaunay = CartesiantoDelaunay ∘ USM6toCartesian
DelaunaytoUSM6 = CartesiantoUSM6 ∘ DelaunaytoCartesian

# ~~~~~~~~~~~~~~~ USMEM <=> Delaunay ~~~~~~~~~~~~~~~ #
USMEMtoDelaunay = CartesiantoDelaunay ∘ USMEMtoCartesian
DelaunaytoUSMEM = CartesiantoUSMEM ∘ DelaunaytoCartesian

# ~~~~~~~~~~~~~~~ Milankovich <=> Delaunay ~~~~~~~~~~~~~~~ #
MilankovichtoDelaunay = CartesiantoDelaunay ∘ MilankovichtoCartesian
DelaunaytoMilankovich = CartesiantoMilankovich ∘ DelaunaytoCartesian

# ~~~~~~~~~~~~~~~ Modified Equinoctial <=> Delaunay ~~~~~~~~~~~~~~~ #
ModifiedEquinoctialtoDelaunay = CartesiantoDelaunay ∘ ModifiedEquinoctialtoCartesian
DelaunaytoModifiedEquinoctial = CartesiantoModifiedEquinoctial ∘ DelaunaytoCartesian

# ~~~~~~~~~~~~~~~ Cylindrical <=> Delaunay ~~~~~~~~~~~~~~~ #
CylindricaltoDelaunay = CartesiantoDelaunay ∘ CylindricaltoCartesian
DelaunaytoCylindrical = CartesiantoCylindrical ∘ DelaunaytoCartesian

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~ Additional Constructors ~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#TODO: NOT MAINTAINABLE
#TODO: CAN I DO THIS WITH METAPROGRAMMING?
# ~~~~~~~~~~~~~~~ Cartesian ~~~~~~~~~~~~~~~ #
Cartesian(X::Cartesian, μ::Number) = X
Cartesian(X::Keplerian, μ::Number) = KepleriantoCartesian(X, μ)
Cartesian(X::USM7, μ::Number) = USM7toCartesian(X, μ)
Cartesian(X::USM6, μ::Number) = USM6toCartesian(X, μ)
Cartesian(X::USMEM, μ::Number) = USMEMtoCartesian(X, μ)
Cartesian(X::Milankovich, μ::Number) = MilankovichtoCartesian(X, μ)
Cartesian(X::ModEq, μ::Number) = ModifiedEquinoctialtoCartesian(X, μ)
Cartesian(X::Cylindrical, μ::Number) = CylindricaltoCartesian(X, μ)
Cartesian(X::Spherical, μ::Number) = SphericaltoCartesian(X, μ)
Cartesian(X::Delaunay, μ::Number) = DelaunaytoCartesian(X, μ)

# ~~~~~~~~~~~~~~~ Keplerian ~~~~~~~~~~~~~~~ #
Keplerian(X::Keplerian, μ::Number) = X
Keplerian(X::Cartesian, μ::Number) = CartesiantoKeplerian(X, μ)
Keplerian(X::USM7, μ::Number) = USM7toKeplerian(X, μ)
Keplerian(X::USM6, μ::Number) = USM6toKeplerian(X, μ)
Keplerian(X::USMEM, μ::Number) = USMEMtoKeplerian(X, μ)
Keplerian(X::Milankovich, μ::Number) = MilankovichtoKeplerian(X, μ)
Keplerian(X::ModEq, μ::Number) = ModifiedEquinoctialtoKeplerian(X, μ)
Keplerian(X::Cylindrical, μ::Number) = CylindricaltoKeplerian(X, μ)
Keplerian(X::Spherical, μ::Number) = SphericaltoKeplerian(X, μ)
Keplerian(X::Delaunay, μ::Number) = DelaunaytoKeplerian(X, μ)

# ~~~~~~~~~~~~~~~ USM7 ~~~~~~~~~~~~~~~ #
USM7(X::USM7, μ::Number) = X
USM7(X::Cartesian, μ::Number) = CartesiantoUSM7(X, μ)
USM7(X::Keplerian, μ::Number) = KepleriantoUSM7(X, μ)
USM7(X::USM6, μ::Number) = USM6toUSM7(X, μ)
USM7(X::USMEM, μ::Number) = USMEMtoUSM7(X, μ)
USM7(X::Milankovich, μ::Number) = MilankovichtoUSM7(X, μ)
USM7(X::ModEq, μ::Number) = ModifiedEquinoctialtoUSM7(X, μ)
USM7(X::Cylindrical, μ::Number) = CylindricaltoUSM7(X, μ)
USM7(X::Spherical, μ::Number) = SphericaltoUSM7(X, μ)
USM7(X::Delaunay, μ::Number) = DelaunaytoUSM7(X, μ)

# ~~~~~~~~~~~~~~~ USM6 ~~~~~~~~~~~~~~~ #
USM6(X::USM6, μ::Number) = X
USM6(X::Cartesian, μ::Number) = CartesiantoUSM6(X, μ)
USM6(X::Keplerian, μ::Number) = KepleriantoUSM6(X, μ)
USM6(X::USM7, μ::Number) = USM7toUSM6(X, μ)
USM6(X::USMEM, μ::Number) = USMEMtoUSM6(X, μ)
USM6(X::Milankovich, μ::Number) = MilankovichtoUSM6(X, μ)
USM6(X::ModEq, μ::Number) = ModifiedEquinoctialtoUSM6(X, μ)
USM6(X::Cylindrical, μ::Number) = CylindricaltoUSM6(X, μ)
USM6(X::Spherical, μ::Number) = SphericaltoUSM6(X, μ)
USM6(X::Delaunay, μ::Number) = DelaunaytoUSM6(X, μ)

# ~~~~~~~~~~~~~~~ USMEM ~~~~~~~~~~~~~~~ #
USMEM(X::USMEM, μ::Number) = X
USMEM(X::Cartesian, μ::Number) = CartesiantoUSMEM(X, μ)
USMEM(X::Keplerian, μ::Number) = KepleriantoUSMEM(X, μ)
USMEM(X::USM7, μ::Number) = USM7toUSMEM(X, μ)
USMEM(X::USM6, μ::Number) = USM6toUSMEM(X, μ)
USMEM(X::Milankovich, μ::Number) = MilankovichtoUSMEM(X, μ)
USMEM(X::ModEq, μ::Number) = ModifiedEquinoctialtoUSMEM(X, μ)
USMEM(X::Cylindrical, μ::Number) = CylindricaltoUSMEM(X, μ)
USMEM(X::Spherical, μ::Number) = SphericaltoUSMEM(X, μ)
USMEM(X::Delaunay, μ::Number) = DelaunaytoUSMEM(X, μ)

# ~~~~~~~~~~~~~~~ Milankovich ~~~~~~~~~~~~~~~ #
Milankovich(X::Milankovich, μ::Number) = X
Milankovich(X::Cartesian, μ::Number) = CartesiantoMilankovich(X, μ)
Milankovich(X::Keplerian, μ::Number) = KepleriantoMilankovich(X, μ)
Milankovich(X::USM7, μ::Number) = USM7toMilankovich(X, μ)
Milankovich(X::USM6, μ::Number) = USM6toMilankovich(X, μ)
Milankovich(X::USMEM, μ::Number) = USMEMtoMilankovich(X, μ)
Milankovich(X::ModEq, μ::Number) = ModifiedEquinoctialtoMilankovich(X, μ)
Milankovich(X::Cylindrical, μ::Number) = CylindricaltoMilankovich(X, μ)
Milankovich(X::Spherical, μ::Number) = SphericaltoMilankovich(X, μ)
Milankovich(X::Delaunay, μ::Number) = DelaunaytoMilankovich(X, μ)

# ~~~~~~~~~~~~~~~ ModEq ~~~~~~~~~~~~~~~ #
ModEq(X::ModEq, μ::Number) = X
ModEq(X::Cartesian, μ::Number) = CartesiantoModifiedEquinoctial(X, μ)
ModEq(X::Keplerian, μ::Number) = KepleriantoModifiedEquinoctial(X, μ)
ModEq(X::USM7, μ::Number) = USM7toModifiedEquinoctial(X, μ)
ModEq(X::USM6, μ::Number) = USM6toModifiedEquinoctial(X, μ)
ModEq(X::USMEM, μ::Number) = USMEMtoModifiedEquinoctial(X, μ)
ModEq(X::Milankovich, μ::Number) = MilankovichtoModifiedEquinoctial(X, μ)
ModEq(X::Cylindrical, μ::Number) = CylindricaltoModifiedEquinoctial(X, μ)
ModEq(X::Spherical, μ::Number) = SphericaltoModifiedEquinoctial(X, μ)
ModEq(X::Delaunay, μ::Number) = DelaunaytoModifiedEquinoctial(X, μ)

# ~~~~~~~~~~~~~~~ Cylindrical ~~~~~~~~~~~~~~~ #
Cylindrical(X::Cylindrical, μ::Number) = X
Cylindrical(X::Cartesian, μ::Number) = CartesiantoCylindrical(X, μ)
Cylindrical(X::Keplerian, μ::Number) = KepleriantoCylindrical(X, μ)
Cylindrical(X::USM7, μ::Number) = USM7toCylindrical(X, μ)
Cylindrical(X::USM6, μ::Number) = USM6toCylindrical(X, μ)
Cylindrical(X::USMEM, μ::Number) = USMEMtoCylindrical(X, μ)
Cylindrical(X::Milankovich, μ::Number) = MilankovichtoCylindrical(X, μ)
Cylindrical(X::ModEq, μ::Number) = ModifiedEquinoctialtoCylindrical(X, μ)
Cylindrical(X::Spherical, μ::Number) = SphericaltoCylindrical(X, μ)
Cylindrical(X::Delaunay, μ::Number) = DelaunaytoCylindrical(X, μ)

# ~~~~~~~~~~~~~~~ Spherical ~~~~~~~~~~~~~~~ #
Spherical(X::Spherical, μ::Number) = X
Spherical(X::Cartesian, μ::Number) = CartesiantoSpherical(X, μ)
Spherical(X::Keplerian, μ::Number) = KepleriantoSpherical(X, μ)
Spherical(X::USM7, μ::Number) = USM7toSpherical(X, μ)
Spherical(X::USM6, μ::Number) = USM6toSpherical(X, μ)
Spherical(X::USMEM, μ::Number) = USMEMtoSpherical(X, μ)
Spherical(X::Milankovich, μ::Number) = MilankovichtoSpherical(X, μ)
Spherical(X::ModEq, μ::Number) = ModifiedEquinoctialtoSpherical(X, μ)
Spherical(X::Cylindrical, μ::Number) = CylindricaltoSpherical(X, μ)
Spherical(X::Delaunay, μ::Number) = DelaunaytoSpherical(X, μ)

# ~~~~~~~~~~~~~~~ Delaunay ~~~~~~~~~~~~~~~ #
Delaunay(X::Delaunay, μ::Number) = X
Delaunay(X::Cartesian, μ::Number) = CartesiantoDelaunay(X, μ)
Delaunay(X::Keplerian, μ::Number) = KepleriantoDelaunay(X, μ)
Delaunay(X::USM7, μ::Number) = USM7toDelaunay(X, μ)
Delaunay(X::USM6, μ::Number) = USM6toDelaunay(X, μ)
Delaunay(X::USMEM, μ::Number) = USMEMtoDelaunay(X, μ)
Delaunay(X::Milankovich, μ::Number) = MilankovichtoDelaunay(X, μ)
Delaunay(X::ModEq, μ::Number) = ModifiedEquinoctialtoDelaunay(X, μ)
Delaunay(X::Cylindrical, μ::Number) = CylindricaltoDelaunay(X, μ)
Delaunay(X::Spherical, μ::Number) = SphericaltoDelaunay(X, μ)
