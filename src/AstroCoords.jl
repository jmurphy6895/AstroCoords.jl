module AstroCoords

using LinearAlgebra
using StaticArraysCore

include("utils.jl")
include("./anomalies.jl")
include("./coordinate_changes.jl")
include("./core_types.jl")

include("./coordinate_sets/cartesian.jl")
include("./coordinate_sets/delauney.jl")
include("./coordinate_sets/keplerian.jl")
include("./coordinate_sets/milankovich.jl")
include("./coordinate_sets/modEq.jl")
include("./coordinate_sets/spherical.jl")
include("./coordinate_sets/usm.jl")

include("./transformations.jl")

include("./quantities.jl")

end
