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

    unitv1 = v1 ./ √(sum(abs2.(v1)))
    unitv2 = v2 ./ √(sum(abs2.(v2)))

    y = unitv1 - unitv2
    x = unitv1 + unitv2

    a = 2.0 * atan(norm(y), norm(x))

    angle::T =
        !(signbit(a) || signbit(float(T)(π) - a)) ? a : (signbit(a) ? zero(T) : float(T)(π))

    return angle
end
