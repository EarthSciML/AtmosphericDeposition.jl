# Dry Deposition
## Model
This is an implementation of a box model used to calculate changes in gas species concentration due to dry deposition.

## Running the model
Here's an example of how concentration of different species, such as SO₂, O₃, NO₂, NO, H₂O₂ and CH₂O change due to dry deposition. 

We can create an instance of the model in the following manner:
```julia
using AtmosphericDeposition
using ModelingToolkit
using DifferentialEquations
using EarthSciMLBase
using Unitful

@parameters t [unit = u"s", description="Time"]
model = DrydepositionG(t)
```
Before running any simulations with the model we need to convert it into a system of differential equations.
```julia
sys = structural_simplify(get_mtk(model))
tspan = (0.0, 3600*24) 
u0 = [2.0,10.0,5,5,2.34,0.15] # initial concentrations of SO₂, O₃, NO₂, NO, H₂O₂, CH₂O
sol = solve(ODEProblem(sys, u0, tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0) # default parameters
```

```@setup 1
using AtmosphericDeposition
using ModelingToolkit
using DifferentialEquations
using EarthSciMLBase
using Unitful

@parameters t [unit = u"s", description="Time"]
model = DrydepositionG(t)

sys = structural_simplify(get_mtk(model))
tspan = (0.0, 3600*24) 
u0 = [2.0,10.0,5,5,2.34,0.15] # initial concentrations of SO₂, O₃, NO₂, NO, H₂O₂, CH₂O
sol = solve(ODEProblem(sys, u0, tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0) # default parameters
```

which we can plot as
```@example 1
using Plots
plot(sol, xlabel="Time (second)", ylabel="concentration (ppb)", legend=:outerright)
```

## Parameters
The parameters in the model are:
```julia
parameters(sys) # [z, z₀, u_star, L, ρA, G, T, θ]
```
where ```z``` is the top of the surface layer [m], ```z₀``` is the roughness length [m], ```u_star``` is friction velocity [m/s], and ```L``` is Monin-Obukhov length [m], ```ρA``` is air density [kg/m3], ```T``` is surface air temperature [K], ```G``` is solar irradiation [W m-2], ```Θ``` is the slope of the local terrain [radians].

Let's run some simulation with different value for parameter ```z```. 
```@example 1
@unpack O3 = sys

p1 = [50,0.04,0.44,0,1.2,300,298,0]
p2 = [10,0.04,0.44,0,1.2,300,298,0]
sol1 = solve(ODEProblem(sys, u0, tspan, p1),AutoTsit5(Rosenbrock23()), saveat=10.0)
sol2 = solve(ODEProblem(sys, u0, tspan, p2),AutoTsit5(Rosenbrock23()), saveat=10.0)

plot([sol1[O3],sol2[O3]], label = ["z=50m" "z=10m"], title = "Change of O3 concentration due to dry deposition", xlabel="Time (second)", ylabel="concentration (ppb)")
```
