module DepositionMTK

using Unitful
using StaticArrays
using ModelingToolkit
using EarthSciMLBase
using IfElse

include("wesley1989.jl")
include("dry_deposition.jl")
include("wet_deposition.jl")

end