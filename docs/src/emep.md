# EMEP Wet Deposition Algorithm
## Model
This is an implementation of a box model used to calculate wet deposition based on formulas at [Unified EMEP Model documentation, chapter 9](https://www.emep.int/publ/reports/2003/emep_report_1_part1_2003.pdf).

## Running the model
```@example 1
using AtmosphericDeposition
using ModelingToolkit
using DifferentialEquations
using EarthSciMLBase
using Unitful

@parameters t [unit = u"s", description="Time"]
model = Wetdeposition(t)
```

Before running any simulations with the model we need to convert it into a system of differential equations.
```@example 1
sys = structural_simplify(get_mtk(model))
tspan = (0.0, 3600*24)
u0 = [2.0,10.0,5,1400,275,50,0.15]  # initial concentration of SO₂, O₃, NO₂, CH₄, CO, DMS, ISOP
prob = ODEProblem(sys, u0, tspan, [])
sol = solve(prob,AutoTsit5(Rosenbrock23()), saveat=10.0) # default parameters
```
which we can plot as
```@example 1
using Plots
plot(sol, xlabel="Time (second)", ylabel="concentration (ppb)", legend=:outerright)
```

## Parameters
The parameters in the model are:
```@example 1
parameters(sys) # [cloudFrac, qrain, ρ_air, Δz]
```
where ```cloudFrac``` is fraction of grid cell covered by clouds, ```qrain``` is rain mixing ratio, ```ρ_air``` is air density [kg/m3], and ```Δz``` is fall distance [m].

Let's run some simulation with different value for parameter ```cloudFrac```. 
```@example 1
@unpack O3 = sys

p1 = [0.3,0.5,1.204,200]
p2 = [0.6,0.5,1.204,200]
sol1 = solve(ODEProblem(sys, u0, tspan, p1),AutoTsit5(Rosenbrock23()), saveat=10.0)
sol2 = solve(ODEProblem(sys, u0, tspan, p2),AutoTsit5(Rosenbrock23()), saveat=10.0)

plot([sol1[O3],sol2[O3]], label = ["cloudFrac=0.3" "cloudFrac=0.6"], title = "Change of O3 concentration due to wet deposition", xlabel="Time (second)", ylabel="concentration (ppb)")
```
From the plot we could see that with larger cloud fraction, the wet deposition rate increases. 

Let's run some simulation with different value for parameter ```qrain``` 
```@example 1
p3 = [0.5,0.3,1.204,200]
p4 = [0.5,0.6,1.204,200]
sol3 = solve(ODEProblem(sys, u0, tspan, p3),AutoTsit5(Rosenbrock23()), saveat=10.0)
sol4 = solve(ODEProblem(sys, u0, tspan, p4),AutoTsit5(Rosenbrock23()), saveat=10.0)

plot([sol3[O3],sol4[O3]], label = ["cloudFrac=0.3" "cloudFrac=0.6"], title = "Change of O3 concentration due to wet deposition", xlabel="Time (second)", ylabel="concentration (ppb)")
```
The graph indicates that an increase in the rain mixing ratio leads to a corresponding rise in the rate of wet deposition.

