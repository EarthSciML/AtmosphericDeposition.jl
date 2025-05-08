module AtmosphericDeposition

using DynamicQuantities
using StaticArrays
using ModelingToolkit
using EarthSciMLBase
using DataInterpolations
using ModelingToolkit: t, D

@register_unit ppb 1

include("wesley1989.jl")
include("dry_deposition.jl")
include("wet_deposition.jl")

end
