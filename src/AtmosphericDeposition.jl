module AtmosphericDeposition

using Unitful
using StaticArrays
using ModelingToolkit
using EarthSciMLBase
using IfElse
using DataInterpolations

include("wesley1989.jl")
include("dry_deposition.jl")
include("wet_deposition.jl")

end