using AtmosphericDeposition
using ModelingToolkit
using DifferentialEquations
using EarthSciMLBase
using Unitful

@parameters t [unit = u"s", description="Time"]
model = Wetdeposition(t)