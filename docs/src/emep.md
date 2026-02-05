# EMEP Wet Deposition Algorithm

## Model

This is an implementation of a box model used to calculate wet deposition based on formulas at [Unified EMEP Model documentation, chapter 9](https://www.emep.int/publ/reports/2003/emep_report_1_part1_2003.pdf).

## Running the model

```julia
using AtmosphericDeposition
using ModelingToolkit
using DifferentialEquations
using EarthSciMLBase
using DynamicQuantities
using ModelingToolkit: t

model = WetDeposition()
```

Before running any simulations with the model we need to convert it into a system of differential equations.

```julia
sys = structural_simplify(model)
tspan = (0.0, 3600*24)
u0 = [2.0, 10.0, 5, 1400, 275, 50, 0.15, 2.34, 10, 0.15]  # initial concentration of SO₂, O₃, NO₂, CH₄, CO, DMS, ISOP, H₂O₂, HNO₃, CH₂O
prob = ODEProblem(sys, u0, tspan, [])
sol = solve(prob, AutoTsit5(Rosenbrock23()), saveat = 10.0) # default parameters
```

```@setup 1
using AtmosphericDeposition
using ModelingToolkit
using DifferentialEquations
using EarthSciMLBase
using DynamicQuantities
using ModelingToolkit:t

model = WetDeposition()

sys = structural_simplify(model)
tspan = (0.0, 3600*24)
prob = ODEProblem(sys, [], tspan, []) # default initial concentration of SO₂, O₃, NO₂, H₂O₂, HNO₃, CH₂O
sol = solve(prob,AutoTsit5(Rosenbrock23()), saveat=10.0) # default parameters
```

which we can plot as

```@example 1
using Plots
plot(sol, xlabel = "Time (second)", ylabel = "concentration (ppb)", legend = :outerright)
```

## Parameters

The parameters in the model are:

```julia @example 1
parameters(sys) #[ρ_air, qrain, Δz, cloudFrac]
```

where `cloudFrac` is fraction of grid cell covered by clouds, `qrain` is rain mixing ratio, `ρ_air` is air density [kg/m3], and `Δz` is fall distance [m].

Let's run some simulation with different value for parameter `cloudFrac`.

```@example 1
@unpack O3, cloudFrac, qrain = sys

p1 = [cloudFrac=>0.3]
p2 = [cloudFrac=>0.6]
sol1 = solve(ODEProblem(sys, [], tspan, p1), AutoTsit5(Rosenbrock23()), saveat = 10.0)
sol2 = solve(ODEProblem(sys, [], tspan, p2), AutoTsit5(Rosenbrock23()), saveat = 10.0)

plot([sol1[O3], sol2[O3]], label = ["cloudFrac=0.3" "cloudFrac=0.6"],
    title = "Change of O3 concentration due to wet deposition",
    xlabel = "Time (second)", ylabel = "concentration (ppb)")
```

From the plot we could see that with larger cloud fraction, the wet deposition rate increases.

Let's run some simulation with different value for parameter `qrain`

```@example 1
p3 = [qrain=>0.3]
p4 = [qrain=>0.6]
sol3 = solve(ODEProblem(sys, [], tspan, p3), AutoTsit5(Rosenbrock23()), saveat = 10.0)
sol4 = solve(ODEProblem(sys, [], tspan, p4), AutoTsit5(Rosenbrock23()), saveat = 10.0)

plot([sol3[O3], sol4[O3]], label = ["qrain=0.3" "qrain=0.6"],
    title = "Change of O3 concentration due to wet deposition",
    xlabel = "Time (second)", ylabel = "concentration (ppb)")
```

The graph indicates that an increase in the rain mixing ratio leads to a corresponding rise in the rate of wet deposition.
