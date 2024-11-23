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
struct CartesianToKeplerianTransform <: AstroCoordTransformation end

@inline function (::CartesianToKeplerianTransform)(
    x::Cartesian{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Keplerian{RT}(cart2koe(params(x), μ))
end

const CartesianToKeplerian = CartesianToKeplerianTransform()

struct KeplerianToCartesianTransform <: AstroCoordTransformation end

@inline function (::KeplerianToCartesianTransform)(
    x::Keplerian{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Cartesian{RT}(koe2cart(params(x), μ))
end

const KeplerianToCartesian = KeplerianToCartesianTransform()

Base.inv(::CartesianToKeplerianTransform) = KeplerianToCartesianTransform()
Base.inv(::KeplerianToCartesianTransform) = CartesianToKeplerianTransform()

# ~~~~~~~~~~~~~~~ Keplerian <=> USM7 ~~~~~~~~~~~~~~~ #
struct KeplerianToUSM7Transform <: AstroCoordTransformation end

@inline function (::KeplerianToUSM7Transform)(
    x::Keplerian{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return USM7{RT}(koe2USM7(params(x), μ))
end

const KeplerianToUSM7 = KeplerianToUSM7Transform()

struct USM7ToKeplerianTransform <: AstroCoordTransformation end

@inline function (::USM7ToKeplerianTransform)(x::USM7{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Keplerian{RT}(USM72koe(params(x), μ))
end

const USM7ToKeplerian = USM7ToKeplerianTransform()

Base.inv(::KeplerianToUSM7Transform) = USM7ToKeplerianTransform()
Base.inv(::USM7ToKeplerianTransform) = KeplerianToUSM7Transform()

# ~~~~~~~~~~~~~~~ Cartesian <=> USM7 ~~~~~~~~~~~~~~~ #
const CartesianToUSM7 = KeplerianToUSM7 ∘ CartesianToKeplerian
const USM7ToCartesian = KeplerianToCartesian ∘ USM7ToKeplerian

# ~~~~~~~~~~~~~~~ USM6 <=> USM7 ~~~~~~~~~~~~~~~ #
struct USM6ToUSM7Transform <: AstroCoordTransformation end

@inline function (::USM6ToUSM7Transform)(x::USM6{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return USM7{RT}(USM62USM7(params(x), μ))
end

const USM6ToUSM7 = USM6ToUSM7Transform()

struct USM7ToUSM6Transform <: AstroCoordTransformation end

@inline function (::USM7ToUSM6Transform)(x::USM7{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return USM6{RT}(USM72USM6(params(x), μ))
end

const USM7ToUSM6 = USM7ToUSM6Transform()

Base.inv(::USM6ToUSM7Transform) = USM7ToUSM6Transform()
Base.inv(::USM7ToUSM6Transform) = USM6ToUSM7Transform()

# ~~~~~~~~~~~~~~~ Cartesian <=> USM6 ~~~~~~~~~~~~~~~ #
const CartesianToUSM6 = USM7ToUSM6 ∘ CartesianToUSM7
const USM6ToCartesian = USM7ToCartesian ∘ USM6ToUSM7

# ~~~~~~~~~~~~~~~ Keplerian <=> USM6 ~~~~~~~~~~~~~~~ #
const KeplerianToUSM6 = USM7ToUSM6 ∘ KeplerianToUSM7
const USM6ToKeplerian = USM7ToKeplerian ∘ USM6ToUSM7

# ~~~~~~~~~~~~~~~ USMEM <=> USM7 ~~~~~~~~~~~~~~~ #
struct USMEMToUSM7Transform <: AstroCoordTransformation end

@inline function (::USMEMToUSM7Transform)(x::USMEM{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return USM7{RT}(USMEM2USM7(params(x), μ))
end

const USMEMToUSM7 = USMEMToUSM7Transform()

struct USM7ToUSMEMTransform <: AstroCoordTransformation end

@inline function (::USM7ToUSMEMTransform)(x::USM7{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return USMEM{RT}(USM72USMEM(params(x), μ))
end

const USM7ToUSMEM = USM7ToUSMEMTransform()

Base.inv(::USMEMToUSM7Transform) = USM7ToUSMEMTransform()
Base.inv(::USM7ToUSMEMTransform) = USMEMToUSM7Transform()

# ~~~~~~~~~~~~~~~ Cartesian <=> USMEM ~~~~~~~~~~~~~~~ #
const CartesianToUSMEM = USM7ToUSMEM ∘ CartesianToUSM7
const USMEMToCartesian = USM7ToCartesian ∘ USMEMToUSM7

# ~~~~~~~~~~~~~~~ Keplerian <=> USMEM ~~~~~~~~~~~~~~~ #
const KeplerianToUSMEM = USM7ToUSMEM ∘ KeplerianToUSM7
const USMEMToKeplerian = USM7ToKeplerian ∘ USMEMToUSM7

# ~~~~~~~~~~~~~~~ USM6 <=> USMEM ~~~~~~~~~~~~~~~ #
const USM6ToUSMEM = USM7ToUSMEM ∘ USM6ToUSM7
const USMEMToUSM6 = USM7ToUSM6 ∘ USMEMToUSM7

# ~~~~~~~~~~~~~~~ Cartesian <=> Milankovich ~~~~~~~~~~~~~~~ #
struct CartesianToMilankovichTransform <: AstroCoordTransformation end

@inline function (::CartesianToMilankovichTransform)(
    x::Cartesian{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Milankovich{RT}(cart2Mil(params(x), μ))
end

const CartesianToMilankovich = CartesianToMilankovichTransform()

struct MilankovichToCartesianTransform <: AstroCoordTransformation end

@inline function (::MilankovichToCartesianTransform)(
    x::Milankovich{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Cartesian{RT}(Mil2cart(params(x), μ))
end

const MilankovichToCartesian = MilankovichToCartesianTransform()

Base.inv(::MilankovichToCartesianTransform) = CartesianToMilankovichTransform()
Base.inv(::CartesianToMilankovichTransform) = MilankovichToCartesianTransform()

# ~~~~~~~~~~~~~~~ Keplerian <=> Milankovich ~~~~~~~~~~~~~~~ #
const KeplerianToMilankovich = CartesianToMilankovich ∘ KeplerianToCartesian
const MilankovichToKeplerian = CartesianToKeplerian ∘ MilankovichToCartesian

# ~~~~~~~~~~~~~~~ USM7 <=> Milankovich ~~~~~~~~~~~~~~~ #
const USM7ToMilankovich = CartesianToMilankovich ∘ USM7ToCartesian
const MilankovichToUSM7 = CartesianToUSM7 ∘ MilankovichToCartesian

# ~~~~~~~~~~~~~~~ USM6 <=> Milankovich ~~~~~~~~~~~~~~~ #
const USM6ToMilankovich = CartesianToMilankovich ∘ USM6ToCartesian
const MilankovichToUSM6 = CartesianToUSM6 ∘ MilankovichToCartesian

# ~~~~~~~~~~~~~~~ USMEM <=> Milankovich ~~~~~~~~~~~~~~~ #
const USMEMToMilankovich = CartesianToMilankovich ∘ USMEMToCartesian
const MilankovichToUSMEM = CartesianToUSMEM ∘ MilankovichToCartesian

# ~~~~~~~~~~~~~~~ Keplerian To ModifiedEquinoctial ~~~~~~~~~~~~~~~ #
struct KeplerianToModifiedEquinoctialTransform <: AstroCoordTransformation end

@inline function (::KeplerianToModifiedEquinoctialTransform)(
    x::Keplerian{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return ModEq{RT}(koe2ModEq(params(x), μ))
end

const KeplerianToModifiedEquinoctial = KeplerianToModifiedEquinoctialTransform()

struct ModifiedEquinoctialToKeplerianTransform <: AstroCoordTransformation end

@inline function (::ModifiedEquinoctialToKeplerianTransform)(
    x::ModEq{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Keplerian{RT}(ModEq2koe(params(x), μ))
end

const ModifiedEquinoctialToKeplerian = ModifiedEquinoctialToKeplerianTransform()

function Base.inv(::KeplerianToModifiedEquinoctialTransform)
    return ModifiedEquinoctialToKeplerianTransform()
end
function Base.inv(::ModifiedEquinoctialToKeplerianTransform)
    return KeplerianToModifiedEquinoctialTransform()
end

# ~~~~~~~~~~~~~~~ Cartesian <=> Modified Equinoctial ~~~~~~~~~~~~~~~ #
const CartesianToModifiedEquinoctial = KeplerianToModifiedEquinoctial ∘ CartesianToKeplerian
const ModifiedEquinoctialToCartesian = KeplerianToCartesian ∘ ModifiedEquinoctialToKeplerian

# ~~~~~~~~~~~~~~~ USM7 <=> Modified Equinoctial ~~~~~~~~~~~~~~~ #
const USM7ToModifiedEquinoctial = KeplerianToModifiedEquinoctial ∘ USM7ToKeplerian
const ModifiedEquinoctialToUSM7 = KeplerianToUSM7 ∘ ModifiedEquinoctialToKeplerian

# ~~~~~~~~~~~~~~~ USM6 <=> Modified Equinoctial ~~~~~~~~~~~~~~~ #
const USM6ToModifiedEquinoctial = KeplerianToModifiedEquinoctial ∘ USM6ToKeplerian
const ModifiedEquinoctialToUSM6 = KeplerianToUSM6 ∘ ModifiedEquinoctialToKeplerian

# ~~~~~~~~~~~~~~~ USMEM <=> Modified Equinoctial ~~~~~~~~~~~~~~~ #
const USMEMToModifiedEquinoctial = KeplerianToModifiedEquinoctial ∘ USMEMToKeplerian
const ModifiedEquinoctialToUSMEM = KeplerianToUSMEM ∘ ModifiedEquinoctialToKeplerian

# ~~~~~~~~~~~~~~~ Milankovich <=> Modified Equinoctial ~~~~~~~~~~~~~~~ #
const MilankovichToModifiedEquinoctial =
    KeplerianToModifiedEquinoctial ∘ MilankovichToKeplerian
const ModifiedEquinoctialToMilankovich =
    KeplerianToMilankovich ∘ ModifiedEquinoctialToKeplerian

# ~~~~~~~~~~~~~~~ Cartesian To Cylindrical ~~~~~~~~~~~~~~~ #
struct CartesianToCylindricalTransform <: AstroCoordTransformation end

@inline function (::CartesianToCylindricalTransform)(
    x::Cartesian{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Cylindrical{RT}(cart2cylind(params(x), μ))
end

const CartesianToCylindrical = CartesianToCylindricalTransform()

struct CylindricalToCartesianTransform <: AstroCoordTransformation end

@inline function (::CylindricalToCartesianTransform)(
    x::Cylindrical{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Cartesian{RT}(cylind2cart(params(x), μ))
end

const CylindricalToCartesian = CylindricalToCartesianTransform()

Base.inv(::CylindricalToCartesianTransform) = CartesianToCylindricalTransform()
Base.inv(::CartesianToCylindricalTransform) = CylindricalToCartesianTransform()

# ~~~~~~~~~~~~~~~ Keplerian <=> Cylindrical ~~~~~~~~~~~~~~~ #
const KeplerianToCylindrical = CartesianToCylindrical ∘ KeplerianToCartesian
const CylindricalToKeplerian = CartesianToKeplerian ∘ CylindricalToCartesian

# ~~~~~~~~~~~~~~~ USM7 <=> Cylindrical ~~~~~~~~~~~~~~~ #
const USM7ToCylindrical = CartesianToCylindrical ∘ USM7ToCartesian
const CylindricalToUSM7 = CartesianToUSM7 ∘ CylindricalToCartesian

# ~~~~~~~~~~~~~~~ USM6 <=> Cylindrical ~~~~~~~~~~~~~~~ #
const USM6ToCylindrical = CartesianToCylindrical ∘ USM6ToCartesian
const CylindricalToUSM6 = CartesianToUSM6 ∘ CylindricalToCartesian

# ~~~~~~~~~~~~~~~ USMEM <=> Cylindrical ~~~~~~~~~~~~~~~ #
const USMEMToCylindrical = CartesianToCylindrical ∘ USMEMToCartesian
const CylindricalToUSMEM = CartesianToUSMEM ∘ CylindricalToCartesian

# ~~~~~~~~~~~~~~~ Milankovich <=> Cylindrical ~~~~~~~~~~~~~~~ #
const MilankovichToCylindrical = CartesianToCylindrical ∘ MilankovichToCartesian
const CylindricalToMilankovich = CartesianToMilankovich ∘ CylindricalToCartesian

# ~~~~~~~~~~~~~~~ Modified Equinoctial <=> Cylindrical ~~~~~~~~~~~~~~~ #
const ModifiedEquinoctialToCylindrical =
    CartesianToCylindrical ∘ ModifiedEquinoctialToCartesian
const CylindricalToModifiedEquinoctial =
    CartesianToModifiedEquinoctial ∘ CylindricalToCartesian

# ~~~~~~~~~~~~~~~ Cartesian To Spherical ~~~~~~~~~~~~~~~ #
struct CartesianToSphericalTransform <: AstroCoordTransformation end

@inline function (::CartesianToSphericalTransform)(
    x::Cartesian{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Spherical{RT}(cart2sphere(params(x), μ))
end

const CartesianToSpherical = CartesianToSphericalTransform()

struct SphericalToCartesianTransform <: AstroCoordTransformation end

@inline function (::SphericalToCartesianTransform)(
    x::Spherical{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Cartesian{RT}(sphere2cart(params(x), μ))
end

const SphericalToCartesian = SphericalToCartesianTransform()

Base.inv(::SphericalToCartesianTransform) = CartesianToSphericalTransform()
Base.inv(::CartesianToSphericalTransform) = SphericalToCartesianTransform()

# ~~~~~~~~~~~~~~~ Keplerian <=> Spherical ~~~~~~~~~~~~~~~ #
const KeplerianToSpherical = CartesianToSpherical ∘ KeplerianToCartesian
const SphericalToKeplerian = CartesianToKeplerian ∘ SphericalToCartesian

# ~~~~~~~~~~~~~~~ USM7 <=> Spherical ~~~~~~~~~~~~~~~ #
const USM7ToSpherical = CartesianToSpherical ∘ USM7ToCartesian
const SphericalToUSM7 = CartesianToUSM7 ∘ SphericalToCartesian

# ~~~~~~~~~~~~~~~ USM6 <=> Spherical ~~~~~~~~~~~~~~~ #
const USM6ToSpherical = CartesianToSpherical ∘ USM6ToCartesian
const SphericalToUSM6 = CartesianToUSM6 ∘ SphericalToCartesian

# ~~~~~~~~~~~~~~~ USMEM <=> Spherical ~~~~~~~~~~~~~~~ #
const USMEMToSpherical = CartesianToSpherical ∘ USMEMToCartesian
const SphericalToUSMEM = CartesianToUSMEM ∘ SphericalToCartesian

# ~~~~~~~~~~~~~~~ Milankovich <=> Spherical ~~~~~~~~~~~~~~~ #
const MilankovichToSpherical = CartesianToSpherical ∘ MilankovichToCartesian
const SphericalToMilankovich = CartesianToMilankovich ∘ SphericalToCartesian

# ~~~~~~~~~~~~~~~ Modified Equinoctial <=> Spherical ~~~~~~~~~~~~~~~ #
const ModifiedEquinoctialToSpherical = CartesianToSpherical ∘ ModifiedEquinoctialToCartesian
const SphericalToModifiedEquinoctial = CartesianToModifiedEquinoctial ∘ SphericalToCartesian

# ~~~~~~~~~~~~~~~ Cylindrical <=> Spherical ~~~~~~~~~~~~~~~ #
const CylindricalToSpherical = CartesianToSpherical ∘ CylindricalToCartesian
const SphericalToCylindrical = CartesianToCylindrical ∘ SphericalToCartesian

# ~~~~~~~~~~~~~~~ Cartesian To Delaunay ~~~~~~~~~~~~~~~ #
struct CartesianToDelaunayTransform <: AstroCoordTransformation end

@inline function (::CartesianToDelaunayTransform)(
    x::Cartesian{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Delaunay{RT}(cart2delaunay(params(x), μ))
end

const CartesianToDelaunay = CartesianToDelaunayTransform()

struct DelaunayToCartesianTransform <: AstroCoordTransformation end

@inline function (::DelaunayToCartesianTransform)(
    x::Delaunay{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Cartesian{RT}(delaunay2cart(params(x), μ))
end

const DelaunayToCartesian = DelaunayToCartesianTransform()

Base.inv(::DelaunayToCartesianTransform) = CartesianToDelaunayTransform()
Base.inv(::CartesianToDelaunayTransform) = DelaunayToCartesianTransform()

# ~~~~~~~~~~~~~~~ Keplerian <=> Delaunay ~~~~~~~~~~~~~~~ #
const KeplerianToDelaunay = CartesianToDelaunay ∘ KeplerianToCartesian
const DelaunayToKeplerian = CartesianToKeplerian ∘ DelaunayToCartesian

# ~~~~~~~~~~~~~~~ USM7 <=> Delaunay ~~~~~~~~~~~~~~~ #
const USM7ToDelaunay = CartesianToDelaunay ∘ USM7ToCartesian
const DelaunayToUSM7 = CartesianToUSM7 ∘ DelaunayToCartesian

# ~~~~~~~~~~~~~~~ USM6 <=> Delaunay ~~~~~~~~~~~~~~~ #
const USM6ToDelaunay = CartesianToDelaunay ∘ USM6ToCartesian
const DelaunayToUSM6 = CartesianToUSM6 ∘ DelaunayToCartesian

# ~~~~~~~~~~~~~~~ USMEM <=> Delaunay ~~~~~~~~~~~~~~~ #
const USMEMToDelaunay = CartesianToDelaunay ∘ USMEMToCartesian
const DelaunayToUSMEM = CartesianToUSMEM ∘ DelaunayToCartesian

# ~~~~~~~~~~~~~~~ Milankovich <=> Delaunay ~~~~~~~~~~~~~~~ #
const MilankovichToDelaunay = CartesianToDelaunay ∘ MilankovichToCartesian
const DelaunayToMilankovich = CartesianToMilankovich ∘ DelaunayToCartesian

# ~~~~~~~~~~~~~~~ Modified Equinoctial <=> Delaunay ~~~~~~~~~~~~~~~ #
const ModifiedEquinoctialToDelaunay = CartesianToDelaunay ∘ ModifiedEquinoctialToCartesian
const DelaunayToModifiedEquinoctial = CartesianToModifiedEquinoctial ∘ DelaunayToCartesian

# ~~~~~~~~~~~~~~~ Cylindrical <=> Delaunay ~~~~~~~~~~~~~~~ #
const CylindricalToDelaunay = CartesianToDelaunay ∘ CylindricalToCartesian
const DelaunayToCylindrical = CartesianToCylindrical ∘ DelaunayToCartesian

# ~~~~~~~~~~~~~~~ Spherical <=> Delaunay ~~~~~~~~~~~~~~~ #
const SphericalToDelaunay = CartesianToDelaunay ∘ SphericalToCartesian
const DelaunayToSpherical = CartesianToSpherical ∘ DelaunayToCartesian

# ~~~~~~~~~~~~~~~ Cartesian To J2 Modified Equinoctial ~~~~~~~~~~~~~~~ #
struct CartesianToJ2EqOETransform <: AstroCoordTransformation end

@inline function (::CartesianToJ2EqOETransform)(
    x::Cartesian{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return J2EqOE{RT}(cart2J2EqOE(params(x), μ))
end

const CartesianToJ2EqOE = CartesianToJ2EqOETransform()

struct J2EqOEToCartesianTransform <: AstroCoordTransformation end

@inline function (::J2EqOEToCartesianTransform)(
    x::J2EqOE{T}, μ::V
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    return Cartesian{RT}(J2EqOE2cart(params(x), μ))
end

const J2EqOEToCartesian = J2EqOEToCartesianTransform()

Base.inv(::J2EqOEToCartesianTransform) = CartesianTo2EqOETransform()
Base.inv(::CartesianToJ2EqOETransform) = J2EqOEToCartesianTransform()

# ~~~~~~~~~~~~~~~ Keplerian <=> J2EqOE ~~~~~~~~~~~~~~~ #
const KeplerianToJ2EqOE = CartesianToJ2EqOE ∘ KeplerianToCartesian
const J2EqOEToKeplerian = CartesianToKeplerian ∘ J2EqOEToCartesian

# ~~~~~~~~~~~~~~~ USM7 <=> J2EqOE ~~~~~~~~~~~~~~~ #
const USM7ToJ2EqOE = CartesianToJ2EqOE ∘ USM7ToCartesian
const J2EqOEToUSM7 = CartesianToUSM7 ∘ J2EqOEToCartesian

# ~~~~~~~~~~~~~~~ USM6 <=> J2EqOE ~~~~~~~~~~~~~~~ #
const USM6ToJ2EqOE = CartesianToJ2EqOE ∘ USM6ToCartesian
const J2EqOEToUSM6 = CartesianToUSM6 ∘ J2EqOEToCartesian

# ~~~~~~~~~~~~~~~ USMEM <=> J2EqOE ~~~~~~~~~~~~~~~ #
const USMEMToJ2EqOE = CartesianToJ2EqOE ∘ USMEMToCartesian
const J2EqOEToUSMEM = CartesianToUSMEM ∘ J2EqOEToCartesian

# ~~~~~~~~~~~~~~~ Milankovich <=> J2EqOE ~~~~~~~~~~~~~~~ #
const MilankovichToJ2EqOE = CartesianToJ2EqOE ∘ MilankovichToCartesian
const J2EqOEToMilankovich = CartesianToMilankovich ∘ J2EqOEToCartesian

# ~~~~~~~~~~~~~~~ Modified Equinoctial <=> J2EqOE ~~~~~~~~~~~~~~~ #
const ModifiedEquinoctialToJ2EqOE = CartesianToJ2EqOE ∘ ModifiedEquinoctialToCartesian
const J2EqOEToModifiedEquinoctial = CartesianToModifiedEquinoctial ∘ J2EqOEToCartesian

# ~~~~~~~~~~~~~~~ Cylindrical <=> J2EqOE ~~~~~~~~~~~~~~~ #
const CylindricalToJ2EqOE = CartesianToJ2EqOE ∘ CylindricalToCartesian
const J2EqOEToCylindrical = CartesianToCylindrical ∘ J2EqOEToCartesian

# ~~~~~~~~~~~~~~~ Spherical <=> J2EqOE ~~~~~~~~~~~~~~~ #
const SphericalToJ2EqOE = CartesianToJ2EqOE ∘ SphericalToCartesian
const J2EqOEToSpherical = CartesianToSpherical ∘ J2EqOEToCartesian

# ~~~~~~~~~~~~~~~ Delaunay <=> J2EqOE ~~~~~~~~~~~~~~~ #
const DelaunayToJ2EqOE = CartesianToJ2EqOE ∘ DelaunayToCartesian
const J2EqOEToDelaunay = CartesianToDelaunay ∘ J2EqOEToCartesian

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~ Additional Constructors ~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#TODO: NOT MAINTAINABLE
#TODO: CAN I DO THIS WITH METAPROGRAMMING?
# ~~~~~~~~~~~~~~~ Cartesian ~~~~~~~~~~~~~~~ #
Cartesian(X::Cartesian{T}, μ::Number) where {T<:Number} = X
Cartesian(X::Keplerian{T}, μ::Number) where {T<:Number} = KeplerianToCartesian(X, μ)
Cartesian(X::USM7{T}, μ::Number) where {T<:Number} = USM7ToCartesian(X, μ)
Cartesian(X::USM6{T}, μ::Number) where {T<:Number} = USM6ToCartesian(X, μ)
Cartesian(X::USMEM{T}, μ::Number) where {T<:Number} = USMEMToCartesian(X, μ)
Cartesian(X::Milankovich{T}, μ::Number) where {T<:Number} = MilankovichToCartesian(X, μ)
Cartesian(X::ModEq{T}, μ::Number) where {T<:Number} = ModifiedEquinoctialToCartesian(X, μ)
Cartesian(X::Cylindrical{T}, μ::Number) where {T<:Number} = CylindricalToCartesian(X, μ)
Cartesian(X::Spherical{T}, μ::Number) where {T<:Number} = SphericalToCartesian(X, μ)
Cartesian(X::Delaunay{T}, μ::Number) where {T<:Number} = DelaunayToCartesian(X, μ)
Cartesian(X::J2EqOE{T}, μ::Number) where {T<:Number} = J2EqOEToCartesian(X, μ)

# ~~~~~~~~~~~~~~~ Keplerian ~~~~~~~~~~~~~~~ #
Keplerian(X::Keplerian{T}, μ::Number) where {T<:Number} = X
Keplerian(X::Cartesian{T}, μ::Number) where {T<:Number} = CartesianToKeplerian(X, μ)
Keplerian(X::USM7{T}, μ::Number) where {T<:Number} = USM7ToKeplerian(X, μ)
Keplerian(X::USM6{T}, μ::Number) where {T<:Number} = USM6ToKeplerian(X, μ)
Keplerian(X::USMEM{T}, μ::Number) where {T<:Number} = USMEMToKeplerian(X, μ)
Keplerian(X::Milankovich{T}, μ::Number) where {T<:Number} = MilankovichToKeplerian(X, μ)
Keplerian(X::ModEq{T}, μ::Number) where {T<:Number} = ModifiedEquinoctialToKeplerian(X, μ)
Keplerian(X::Cylindrical{T}, μ::Number) where {T<:Number} = CylindricalToKeplerian(X, μ)
Keplerian(X::Spherical{T}, μ::Number) where {T<:Number} = SphericalToKeplerian(X, μ)
Keplerian(X::Delaunay{T}, μ::Number) where {T<:Number} = DelaunayToKeplerian(X, μ)
Keplerian(X::J2EqOE{T}, μ::Number) where {T<:Number} = J2EqOEToKeplerian(X, μ)

# ~~~~~~~~~~~~~~~ USM7 ~~~~~~~~~~~~~~~ #
USM7(X::USM7{T}, μ::Number) where {T<:Number} = X
USM7(X::Cartesian{T}, μ::Number) where {T<:Number} = CartesianToUSM7(X, μ)
USM7(X::Keplerian{T}, μ::Number) where {T<:Number} = KeplerianToUSM7(X, μ)
USM7(X::USM6{T}, μ::Number) where {T<:Number} = USM6ToUSM7(X, μ)
USM7(X::USMEM{T}, μ::Number) where {T<:Number} = USMEMToUSM7(X, μ)
USM7(X::Milankovich{T}, μ::Number) where {T<:Number} = MilankovichToUSM7(X, μ)
USM7(X::ModEq{T}, μ::Number) where {T<:Number} = ModifiedEquinoctialToUSM7(X, μ)
USM7(X::Cylindrical{T}, μ::Number) where {T<:Number} = CylindricalToUSM7(X, μ)
USM7(X::Spherical{T}, μ::Number) where {T<:Number} = SphericalToUSM7(X, μ)
USM7(X::Delaunay{T}, μ::Number) where {T<:Number} = DelaunayToUSM7(X, μ)
USM7(X::J2EqOE{T}, μ::Number) where {T<:Number} = J2EqOEToUSM7(X, μ)

# ~~~~~~~~~~~~~~~ USM6 ~~~~~~~~~~~~~~~ #
USM6(X::USM6{T}, μ::Number) where {T<:Number} = X
USM6(X::Cartesian{T}, μ::Number) where {T<:Number} = CartesianToUSM6(X, μ)
USM6(X::Keplerian{T}, μ::Number) where {T<:Number} = KeplerianToUSM6(X, μ)
USM6(X::USM7{T}, μ::Number) where {T<:Number} = USM7ToUSM6(X, μ)
USM6(X::USMEM{T}, μ::Number) where {T<:Number} = USMEMToUSM6(X, μ)
USM6(X::Milankovich{T}, μ::Number) where {T<:Number} = MilankovichToUSM6(X, μ)
USM6(X::ModEq{T}, μ::Number) where {T<:Number} = ModifiedEquinoctialToUSM6(X, μ)
USM6(X::Cylindrical{T}, μ::Number) where {T<:Number} = CylindricalToUSM6(X, μ)
USM6(X::Spherical{T}, μ::Number) where {T<:Number} = SphericalToUSM6(X, μ)
USM6(X::Delaunay{T}, μ::Number) where {T<:Number} = DelaunayToUSM6(X, μ)
USM6(X::J2EqOE{T}, μ::Number) where {T<:Number} = J2EqOEToUSM6(X, μ)

# ~~~~~~~~~~~~~~~ USMEM ~~~~~~~~~~~~~~~ #
USMEM(X::USMEM{T}, μ::Number) where {T<:Number} = X
USMEM(X::Cartesian{T}, μ::Number) where {T<:Number} = CartesianToUSMEM(X, μ)
USMEM(X::Keplerian{T}, μ::Number) where {T<:Number} = KeplerianToUSMEM(X, μ)
USMEM(X::USM7{T}, μ::Number) where {T<:Number} = USM7ToUSMEM(X, μ)
USMEM(X::USM6{T}, μ::Number) where {T<:Number} = USM6ToUSMEM(X, μ)
USMEM(X::Milankovich{T}, μ::Number) where {T<:Number} = MilankovichToUSMEM(X, μ)
USMEM(X::ModEq{T}, μ::Number) where {T<:Number} = ModifiedEquinoctialToUSMEM(X, μ)
USMEM(X::Cylindrical{T}, μ::Number) where {T<:Number} = CylindricalToUSMEM(X, μ)
USMEM(X::Spherical{T}, μ::Number) where {T<:Number} = SphericalToUSMEM(X, μ)
USMEM(X::Delaunay{T}, μ::Number) where {T<:Number} = DelaunayToUSMEM(X, μ)
USMEM(X::J2EqOE{T}, μ::Number) where {T<:Number} = J2EqOEToUSMEM(X, μ)

# ~~~~~~~~~~~~~~~ Milankovich ~~~~~~~~~~~~~~~ #
Milankovich(X::Milankovich{T}, μ::Number) where {T<:Number} = X
Milankovich(X::Cartesian{T}, μ::Number) where {T<:Number} = CartesianToMilankovich(X, μ)
Milankovich(X::Keplerian{T}, μ::Number) where {T<:Number} = KeplerianToMilankovich(X, μ)
Milankovich(X::USM7{T}, μ::Number) where {T<:Number} = USM7ToMilankovich(X, μ)
Milankovich(X::USM6{T}, μ::Number) where {T<:Number} = USM6ToMilankovich(X, μ)
Milankovich(X::USMEM{T}, μ::Number) where {T<:Number} = USMEMToMilankovich(X, μ)
function Milankovich(X::ModEq{T}, μ::Number) where {T<:Number}
    return ModifiedEquinoctialToMilankovich(X, μ)
end
Milankovich(X::Cylindrical{T}, μ::Number) where {T<:Number} = CylindricalToMilankovich(X, μ)
Milankovich(X::Spherical{T}, μ::Number) where {T<:Number} = SphericalToMilankovich(X, μ)
Milankovich(X::Delaunay{T}, μ::Number) where {T<:Number} = DelaunayToMilankovich(X, μ)
Milankovich(X::J2EqOE{T}, μ::Number) where {T<:Number} = J2EqOEToMilankovich(X, μ)

# ~~~~~~~~~~~~~~~ ModEq ~~~~~~~~~~~~~~~ #
ModEq(X::ModEq{T}, μ::Number) where {T<:Number} = X
ModEq(X::Cartesian{T}, μ::Number) where {T<:Number} = CartesianToModifiedEquinoctial(X, μ)
ModEq(X::Keplerian{T}, μ::Number) where {T<:Number} = KeplerianToModifiedEquinoctial(X, μ)
ModEq(X::USM7{T}, μ::Number) where {T<:Number} = USM7ToModifiedEquinoctial(X, μ)
ModEq(X::USM6{T}, μ::Number) where {T<:Number} = USM6ToModifiedEquinoctial(X, μ)
ModEq(X::USMEM{T}, μ::Number) where {T<:Number} = USMEMToModifiedEquinoctial(X, μ)
function ModEq(X::Milankovich{T}, μ::Number) where {T<:Number}
    return MilankovichToModifiedEquinoctial(X, μ)
end
function ModEq(X::Cylindrical{T}, μ::Number) where {T<:Number}
    return CylindricalToModifiedEquinoctial(X, μ)
end
ModEq(X::Spherical{T}, μ::Number) where {T<:Number} = SphericalToModifiedEquinoctial(X, μ)
ModEq(X::Delaunay{T}, μ::Number) where {T<:Number} = DelaunayToModifiedEquinoctial(X, μ)
ModEq(X::J2EqOE{T}, μ::Number) where {T<:Number} = J2EqOEToModEq(X, μ)

# ~~~~~~~~~~~~~~~ Cylindrical ~~~~~~~~~~~~~~~ #
Cylindrical(X::Cylindrical{T}, μ::Number) where {T<:Number} = X
Cylindrical(X::Cartesian{T}, μ::Number) where {T<:Number} = CartesianToCylindrical(X, μ)
Cylindrical(X::Keplerian{T}, μ::Number) where {T<:Number} = KeplerianToCylindrical(X, μ)
Cylindrical(X::USM7{T}, μ::Number) where {T<:Number} = USM7ToCylindrical(X, μ)
Cylindrical(X::USM6{T}, μ::Number) where {T<:Number} = USM6ToCylindrical(X, μ)
Cylindrical(X::USMEM{T}, μ::Number) where {T<:Number} = USMEMToCylindrical(X, μ)
Cylindrical(X::Milankovich{T}, μ::Number) where {T<:Number} = MilankovichToCylindrical(X, μ)
function Cylindrical(X::ModEq{T}, μ::Number) where {T<:Number}
    return ModifiedEquinoctialToCylindrical(X, μ)
end
Cylindrical(X::Spherical{T}, μ::Number) where {T<:Number} = SphericalToCylindrical(X, μ)
Cylindrical(X::Delaunay{T}, μ::Number) where {T<:Number} = DelaunayToCylindrical(X, μ)
Cylindrical(X::J2EqOE{T}, μ::Number) where {T<:Number} = J2EqOEToCylindrical(X, μ)

# ~~~~~~~~~~~~~~~ Spherical ~~~~~~~~~~~~~~~ #
Spherical(X::Spherical{T}, μ::Number) where {T<:Number} = X
Spherical(X::Cartesian{T}, μ::Number) where {T<:Number} = CartesianToSpherical(X, μ)
Spherical(X::Keplerian{T}, μ::Number) where {T<:Number} = KeplerianToSpherical(X, μ)
Spherical(X::USM7{T}, μ::Number) where {T<:Number} = USM7ToSpherical(X, μ)
Spherical(X::USM6{T}, μ::Number) where {T<:Number} = USM6ToSpherical(X, μ)
Spherical(X::USMEM{T}, μ::Number) where {T<:Number} = USMEMToSpherical(X, μ)
Spherical(X::Milankovich{T}, μ::Number) where {T<:Number} = MilankovichToSpherical(X, μ)
Spherical(X::ModEq{T}, μ::Number) where {T<:Number} = ModifiedEquinoctialToSpherical(X, μ)
Spherical(X::Cylindrical{T}, μ::Number) where {T<:Number} = CylindricalToSpherical(X, μ)
Spherical(X::Delaunay{T}, μ::Number) where {T<:Number} = DelaunayToSpherical(X, μ)
Spherical(X::J2EqOE{T}, μ::Number) where {T<:Number} = J2EqOEToSpherical(X, μ)

# ~~~~~~~~~~~~~~~ Delaunay ~~~~~~~~~~~~~~~ #
Delaunay(X::Delaunay{T}, μ::Number) where {T<:Number} = X
Delaunay(X::Cartesian{T}, μ::Number) where {T<:Number} = CartesianToDelaunay(X, μ)
Delaunay(X::Keplerian{T}, μ::Number) where {T<:Number} = KeplerianToDelaunay(X, μ)
Delaunay(X::USM7{T}, μ::Number) where {T<:Number} = USM7ToDelaunay(X, μ)
Delaunay(X::USM6{T}, μ::Number) where {T<:Number} = USM6ToDelaunay(X, μ)
Delaunay(X::USMEM{T}, μ::Number) where {T<:Number} = USMEMToDelaunay(X, μ)
Delaunay(X::Milankovich{T}, μ::Number) where {T<:Number} = MilankovichToDelaunay(X, μ)
Delaunay(X::ModEq{T}, μ::Number) where {T<:Number} = ModifiedEquinoctialToDelaunay(X, μ)
Delaunay(X::Cylindrical{T}, μ::Number) where {T<:Number} = CylindricalToDelaunay(X, μ)
Delaunay(X::Spherical{T}, μ::Number) where {T<:Number} = SphericalToDelaunay(X, μ)
Delaunay(X::J2EqOE{T}, μ::Number) where {T<:Number} = J2EqOEToDelaunay(X, μ)

# ~~~~~~~~~~~~~~~ J2EqOE ~~~~~~~~~~~~~~~ #
J2EqOE(X::J2EqOE{T}, μ::Number) where {T<:Number} = X
J2EqOE(X::Cartesian{T}, μ::Number) where {T<:Number} = CartesianToJ2EqOE(X, μ)
J2EqOE(X::Keplerian{T}, μ::Number) where {T<:Number} = KeplerianToJ2EqOE(X, μ)
J2EqOE(X::USM7{T}, μ::Number) where {T<:Number} = USM7ToJ2EqOE(X, μ)
J2EqOE(X::USM6{T}, μ::Number) where {T<:Number} = USM6ToJ2EqOE(X, μ)
J2EqOE(X::USMEM{T}, μ::Number) where {T<:Number} = USMEMToJ2EqOE(X, μ)
J2EqOE(X::Milankovich{T}, μ::Number) where {T<:Number} = MilankovichToJ2EqOE(X, μ)
J2EqOE(X::ModEq{T}, μ::Number) where {T<:Number} = ModifiedEquinoctialToJ2EqOE(X, μ)
J2EqOE(X::Cylindrical{T}, μ::Number) where {T<:Number} = CylindricalToJ2EqOE(X, μ)
J2EqOE(X::Spherical{T}, μ::Number) where {T<:Number} = SphericalToJ2EqOE(X, μ)
J2EqOE(X::Delaunay{T}, μ::Number) where {T<:Number} = DelaunayToJ2EqOE(X, μ)
