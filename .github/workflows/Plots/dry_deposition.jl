using AtmosphericDeposition
using ModelingToolkit
using DifferentialEquations
using EarthSciMLBase
using Unitful

@parameters t [unit = u"s", description="Time"]

model = DrydepositionG(t)

sys = structural_simplify(get_mtk(model))
tspan = (0.0, 3600*24) 
u0 = [2.0,10.0,5,5,2.34,0.15] 
sol = solve(ODEProblem(sys, u0, tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)

using Plots
plot(sol,xlabel="Time (second)", ylabel="concentration (ppb)",legend=:outerright)

states(sys)
parameters(sys)

@unpack O3 = sys

p1 = [50,0.04,0.44,0,1.2,300,298,0]
p2 = [10,0.04,0.44,0,1.2,300,298,0]
sol1 = solve(ODEProblem(sys, u0, tspan, p1),AutoTsit5(Rosenbrock23()), saveat=10.0)
sol2 = solve(ODEProblem(sys, u0, tspan, p2),AutoTsit5(Rosenbrock23()), saveat=10.0)

plot([sol1[O3],sol2[O3]], label = ["z=50m" "z=10m"], title = "Change of O3 concentration due to dry deposition", xlabel="Time (second)", ylabel="concentration (ppb)")