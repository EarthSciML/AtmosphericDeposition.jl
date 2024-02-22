# Dry Deposition
## Model
This is an implementation of a box model used to calculate changes in gas species concentration due to dry deposition.

## Running the model
Here's an example of how concentration of different species, such as SO<sub>2</sub>, O<sub>3</sub>, NO<sub>2</sub>, NO, H<sub>2</sub>O<sub>2</sub> and CH<sub>2</sub>O change due to dry deposition. 

We can create an instance of the model in the following manner:
```@example 1
using AtmosphericDeposition
using ModelingToolkit
using DifferentialEquations
using EarthSciMLBase
using Unitful

@parameters t [unit = u"s", description="Time"]

model = DrydepositionG(t)
```
Before running any simulations with the model we need to convert it into a system of differential equations.
```@example 1
sys = structural_simplify(get_mtk(model))
tspan = (0.0, 3600*24) 
u0 = [2.0,10.0,5,5,2.34,0.15] 
sol = solve(ODEProblem(sys, u0, tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)
```
which we can plot as
```@example 1
using Plots
plot(sol, xlabel="Time (second)", ylabel="concentration (ppb)", legend=:outerright)
```