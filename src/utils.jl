export angle_between_vectors
"""
    angle_between_vectors(
        v1::AbstractVector{T1}, v2::AbstractVector{T2}
    ) where {T1<:Number,T2<:Number}

Computes the angle between two vectors in a more numerically stable way than dot product.

# Arguments
-`v1::AbstractVector{<:Number}`: The first vector of the computation
-`v2::AbstractVector{<:Number}`: The second vector of the computation

# Returns 
-`angle::Number`: The angle between the two vectors
"""
@inline function angle_between_vectors(
    v1::AbstractVector{T1}, v2::AbstractVector{T2}
) where {T1<:Number,T2<:Number}
    T = promote_type(T1, T2)

    unitv1 = normalize(v1)
    unitv2 = normalize(v2)

    y = unitv1 - unitv2
    x = unitv1 + unitv2

    a = 2.0 * atan(norm(y), norm(x))

    angle::T =
        !(sign(a) == -1.0 || sign(T(π) - a) == -1.0) ? a : (sign(a) == -1.0 ? zero(T) : T(π))

    return angle
end
