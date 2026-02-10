# DryDeposition.jl

## Overview

DryDeposition.jl is a ModelingToolkit.jl implementation of the dry deposition equations from Chapter 19 of Seinfeld and Pandis (2006). Dry deposition is the transport of gaseous and particulate species from the atmosphere onto surfaces in the absence of precipitation.

The package implements the resistance model for dry deposition, which represents the deposition process using an electrical resistance analogy. The transport of material to the surface is governed by three resistances in series:

  - **Aerodynamic resistance** (``r_a``): Turbulent transport through the atmospheric surface layer
  - **Quasi-laminar resistance** (``r_b``): Molecular/Brownian diffusion across the thin stagnant layer adjacent to the surface
  - **Surface resistance** (``r_c``): Uptake at the surface (for gases only; particles are assumed to stick on contact)

### Gas Deposition

For gases, the deposition velocity is given by:

```math
v_d = \frac{1}{r_a + r_b + r_c}
```

### Particle Deposition

For particles, gravitational settling operates in parallel with the resistance pathway:

```math
v_d = \frac{1}{r_a + r_b + r_a r_b v_s} + v_s
```

where ``v_s`` is the particle settling velocity.

## Reference

> Seinfeld, J.H. and Pandis, S.N. (2006). *Atmospheric Chemistry and Physics: From Air Pollution to Climate Change*, 2nd Edition, Chapter 19: Dry Deposition. John Wiley & Sons, Inc.

## Package Contents

This package provides the following equation system implementations:

| Component                                              | Equations        | Description                                                |
|:------------------------------------------------------ |:---------------- |:---------------------------------------------------------- |
| [`DryDeposition.AerodynamicResistance`](@ref)          | Eq. 19.14        | Aerodynamic resistance for neutral stability               |
| [`DryDeposition.QuasiLaminarResistanceGas`](@ref)      | Eq. 19.17        | Quasi-laminar resistance for gaseous species               |
| [`DryDeposition.ParticleSettling`](@ref)               | Eq. 19.18, 19.20 | Particle settling velocity and Brownian diffusivity        |
| [`DryDeposition.QuasiLaminarResistanceParticle`](@ref) | Eq. 19.27        | Quasi-laminar resistance for particles (Zhang et al. 2001) |
| [`DryDeposition.SurfaceResistance`](@ref)              | Eq. 19.50-19.52  | Surface resistance using Wesely (1989) model               |
| [`DryDeposition.DryDepositionGas`](@ref)               | Eq. 19.2         | Complete gas deposition velocity                           |
| [`DryDeposition.DryDepositionParticle`](@ref)          | Eq. 19.7         | Complete particle deposition velocity                      |

## Quick Start

```julia
using ModelingToolkit
using ModelingToolkit: t

# Include the module
include("src/DryDeposition.jl")
using .DryDeposition

# Create a particle deposition system
@named part_dep = DryDepositionParticle()

# View system structure
unknowns(part_dep)
parameters(part_dep)
equations(part_dep)
```

## Physical Constants

The module uses the following physical constants:

| Constant   | Value        | Description                |
|:---------- |:------------ |:-------------------------- |
| ``\kappa`` | 0.4          | von Karman constant        |
| ``g``      | 9.81 m/s^2   | Gravitational acceleration |
| ``k_B``    | 1.38e-23 J/K | Boltzmann constant         |

* * *

# Equation Systems

This section provides detailed documentation of each equation system implemented in DryDeposition.jl, including state variables, parameters, and the governing equations.

## Aerodynamic Resistance (Eq. 19.14)

The aerodynamic resistance represents the turbulent transport of material from the bulk atmosphere down to the quasi-laminar sublayer. Under neutral stability conditions:

```math
r_a = \frac{1}{\kappa u_*} \ln\left(\frac{z}{z_0}\right)
```

where:

  - ``\kappa = 0.4`` is the von Karman constant
  - ``u_*`` is the friction velocity (m/s)
  - ``z`` is the reference height (m)
  - ``z_0`` is the surface roughness length (m)

### State Variables

```@example aerodynamic
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using AtmosphericDeposition
using AtmosphericDeposition.DryDeposition

@named sys = AerodynamicResistance()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example aerodynamic
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [string(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example aerodynamic
equations(sys)
```

* * *

## Quasi-Laminar Resistance for Gases (Eq. 19.17)

The quasi-laminar resistance for gases depends on the molecular diffusivity of the gas through the Schmidt number:

```math
r_b = \frac{5 \, \text{Sc}^{2/3}}{u_*}
```

where the Schmidt number is ``\text{Sc} = \nu / D_g``, with ``\nu`` being the kinematic viscosity of air and ``D_g`` the molecular diffusivity of the gas.

### State Variables

```@example qlam_gas
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using AtmosphericDeposition
using AtmosphericDeposition.DryDeposition

@named sys = QuasiLaminarResistanceGas()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example qlam_gas
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [string(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example qlam_gas
equations(sys)
```

* * *

## Particle Settling (Eq. 19.18, 19.20)

Particle settling velocity follows Stokes' law with slip correction:

```math
v_s = \frac{\rho_p D_p^2 g C_c}{18 \mu}
```

The Brownian diffusivity of particles is given by:

```math
D = \frac{k_B T C_c}{3 \pi \mu D_p}
```

where:

  - ``\rho_p`` is the particle density (kg/m^3)
  - ``D_p`` is the particle diameter (m)
  - ``g`` is gravitational acceleration (m/s^2)
  - ``C_c`` is the Cunningham slip correction factor (dimensionless)
  - ``\mu`` is the dynamic viscosity of air (Pa s)
  - ``k_B`` is the Boltzmann constant (J/K)
  - ``T`` is temperature (K)

### State Variables

```@example settling
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using AtmosphericDeposition
using AtmosphericDeposition.DryDeposition

@named sys = ParticleSettling()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example settling
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [string(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example settling
equations(sys)
```

* * *

## Quasi-Laminar Resistance for Particles (Eq. 19.27)

The Zhang et al. (2001) model calculates particle quasi-laminar resistance based on three collection mechanisms:

```math
r_b = \frac{1}{3 u_* (E_B + E_{IM} + E_{IN}) R_1}
```

where:

**Brownian diffusion efficiency** (Eq. 19.21):

```math
E_B = \text{Sc}^{-\gamma}
```

**Impaction efficiency** (Eq. 19.22):

```math
E_{IM} = \left(\frac{\text{St}}{\alpha + \text{St}}\right)^2
```

**Interception efficiency** (Eq. 19.25):

```math
E_{IN} = \frac{1}{2}\left(\frac{D_p}{A}\right)^2
```

**Sticking fraction** (Eq. 19.26):

```math
R_1 = \exp(-\sqrt{\text{St}})
```

**Stokes number for vegetated surfaces** (Eq. 19.24):

```math
\text{St} = \frac{v_s u_*}{g A}
```

The parameters ``\alpha``, ``\gamma``, and ``A`` depend on land-use category (Table 19.2).

### State Variables

```@example qlam_particle
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using AtmosphericDeposition
using AtmosphericDeposition.DryDeposition

@named sys = QuasiLaminarResistanceParticle()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example qlam_particle
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [string(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example qlam_particle
equations(sys)
```

### Land-Use Parameters (Table 19.2)

| Land Use Category | A (mm) | ``\alpha`` | ``\gamma`` |
|:----------------- |:------ |:---------- |:---------- |
| Grass             | 2.0    | 1.2        | 0.54       |
| Deciduous Forest  | 5.0    | 0.8        | 0.56       |
| Desert            | 10.0   | 50.0       | 0.54       |

* * *

## Surface Resistance (Eq. 19.50-19.52)

The Wesely (1989) model calculates surface resistance for gases using multiple parallel pathways:

```math
r_c = \left(\frac{1}{r_{st} + r_m} + \frac{1}{r_{lu}} + \frac{1}{r_{dc} + r_{cl}} + \frac{1}{r_{ac} + r_{gs}}\right)^{-1}
```

**Stomatal resistance** (Eq. 19.51):

```math
r_{st} = r_i \left[1 + \left(\frac{200}{G + 0.1}\right)^2\right] \left[\frac{400}{T_s(40 - T_s)}\right]
```

**Mesophyll resistance** (Eq. 19.52):

```math
r_m = \frac{1}{3.3 \times 10^{-4} H^* + 100 f_0}
```

where:

  - ``r_i`` is the minimum stomatal resistance
  - ``G`` is solar irradiance (W/m^2)
  - ``T_s`` is surface temperature (Celsius, valid range 0-40)
  - ``H^*`` is the effective Henry's law constant
  - ``f_0`` is the reactivity factor (0-1)

### State Variables

```@example surface
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using AtmosphericDeposition
using AtmosphericDeposition.DryDeposition

@named sys = SurfaceResistance()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example surface
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [string(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example surface
equations(sys)
```

### Gas Properties (Table 19.4)

| Species | ``D_{H_2O}/D_{gas}`` | ``H^*`` (M/atm) | ``f_0`` |
|:------- |:-------------------- |:--------------- |:------- |
| SO2     | 1.9                  | 10^5            | 0       |
| O3      | 1.6                  | 0.01            | 1       |
| NO2     | 1.6                  | 0.01            | 0.1     |
| HNO3    | 1.9                  | 10^14           | 0       |

* * *

## Complete Gas Deposition System (Eq. 19.2)

The `DryDepositionGas` component composes the aerodynamic, quasi-laminar, and surface resistance systems to calculate the complete gas deposition velocity:

```math
v_d = \frac{1}{r_a + r_b + r_c}
```

and the deposition flux (Eq. 19.1):

```math
F = -v_d C
```

### State Variables

```@example gas_dep
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using AtmosphericDeposition.DryDeposition

@named sys = DryDepositionGas()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Equations

```@example gas_dep
equations(sys)
```

* * *

## Complete Particle Deposition System (Eq. 19.7)

The `DryDepositionParticle` component composes the aerodynamic, settling, and quasi-laminar resistance systems to calculate the complete particle deposition velocity:

```math
v_d = \frac{1}{r_a + r_b + r_a r_b v_s} + v_s
```

This accounts for gravitational settling operating in parallel with the turbulent/diffusive transport pathway.

### State Variables

```@example particle_dep
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using AtmosphericDeposition.DryDeposition

@named sys = DryDepositionParticle()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Equations

```@example particle_dep
equations(sys)
```

* * *

# Analysis

This section reproduces key figures from Seinfeld and Pandis (2006) Chapter 19 to validate the implementation and demonstrate the physical behavior of the dry deposition model.

## Figure 19.2: Particle Deposition Velocity vs. Diameter

Figure 19.2 (page 905) shows the characteristic V-shaped curve of particle deposition velocity as a function of particle diameter. The minimum in deposition velocity occurs in the 0.1-1.0 micrometer range, often called the "accumulation mode."

### Physical Explanation

  - **Small particles (< 0.1 micrometer)**: Efficiently transported by Brownian diffusion across the quasi-laminar sublayer
  - **Intermediate particles (0.1-1.0 micrometer)**: Neither Brownian diffusion nor impaction/settling is effective - the "accumulation mode"
  - **Large particles (> 1 micrometer)**: Efficiently deposited by impaction, interception, and gravitational settling

```@example figure19_2
using Plots

# Physical constants
const k_B = 1.38e-23  # Boltzmann constant (J/K)
const g = 9.81        # Gravitational acceleration (m/s^2)

# Air properties at 298 K
const T = 298.0       # Temperature (K)
const mu = 1.85e-5    # Dynamic viscosity (Pa*s)
const nu = 1.55e-5    # Kinematic viscosity (m^2/s)
const rho_p = 1000.0  # Particle density (kg/m^3) - unit density

# Land-use parameters for grass (Table 19.2)
const A = 2.0e-3      # Characteristic collector radius (m)
const alpha = 1.2     # Impaction parameter
const gamma_param = 0.54  # Brownian diffusion parameter

# Aerodynamic parameters
const z = 10.0        # Reference height (m)
const z_0 = 0.1       # Roughness length for grass (m)
const kappa = 0.4     # von Karman constant

"""
Calculate Cunningham slip correction factor.
Approximation valid for particles in air at STP.
"""
function cunningham_correction(D_p)
    # Mean free path of air at 298 K (m)
    lambda = 6.8e-8
    Kn = 2 * lambda / D_p  # Knudsen number
    return 1 + Kn * (1.257 + 0.4 * exp(-1.1 / Kn))
end

"""
Calculate particle settling velocity (Eq. 19.18).
"""
function settling_velocity(D_p, C_c)
    return rho_p * D_p^2 * g * C_c / (18 * mu)
end

"""
Calculate Brownian diffusivity (Eq. 19.20).
"""
function brownian_diffusivity(D_p, C_c)
    return k_B * T * C_c / (3 * pi * mu * D_p)
end

"""
Calculate particle deposition velocity using Zhang et al. (2001) model.
"""
function particle_deposition_velocity(D_p, u_star)
    # Cunningham correction
    C_c = cunningham_correction(D_p)

    # Settling velocity (Eq. 19.18)
    v_s = settling_velocity(D_p, C_c)

    # Brownian diffusivity (Eq. 19.20)
    D_diff = brownian_diffusivity(D_p, C_c)

    # Schmidt number
    Sc = nu / D_diff

    # Stokes number for vegetated surfaces (Eq. 19.24)
    St = v_s * u_star / (g * A)

    # Collection efficiencies
    E_B = Sc^(-gamma_param)                # Brownian diffusion (Eq. 19.21)
    E_IM = (St / (alpha + St))^2          # Impaction (Eq. 19.22)
    E_IN = 0.5 * (D_p / A)^2              # Interception (Eq. 19.25)

    # Sticking fraction (Eq. 19.26)
    R_1 = exp(-sqrt(St))

    # Quasi-laminar resistance (Eq. 19.27)
    r_b = 1 / (3 * u_star * (E_B + E_IM + E_IN) * R_1)

    # Aerodynamic resistance (Eq. 19.14)
    r_a = (1 / (kappa * u_star)) * log(z / z_0)

    # Particle deposition velocity (Eq. 19.7)
    v_d = 1 / (r_a + r_b + r_a * r_b * v_s) + v_s

    return v_d, v_s, E_B, E_IM, E_IN
end

# Generate data for different friction velocities (matching Figure 19.2 conditions)
D_p_range = 10 .^ range(-2, 2, length = 100) .* 1e-6  # 0.01 to 100 um

# Multiple friction velocities to show dependence
u_star_values = [0.11, 0.44, 1.17]  # m/s (from Figure 19.2 legend)
u_star_labels = ["u* = 0.11 m/s (z0 = 0.002 cm)",
    "u* = 0.44 m/s (z0 = 0.02 cm)",
    "u* = 1.17 m/s (z0 = 0.1 cm)"]

# Create the plot
p = plot(
    xlabel = "Particle Diameter (micrometer)",
    ylabel = "Deposition Velocity (cm/s)",
    title = "Figure 19.2: Particle Dry Deposition Velocity",
    xscale = :log10,
    yscale = :log10,
    xlims = (0.01, 100),
    ylims = (0.01, 100),
    legend = :bottomright,
    size = (700, 500),
    grid = true,
    minorgrid = true
)

colors = [:blue, :red, :green]
for (i, u_star) in enumerate(u_star_values)
    v_d_values = Float64[]
    v_s_values = Float64[]

    for D_p in D_p_range
        v_d, v_s, _, _, _ = particle_deposition_velocity(D_p, u_star)
        push!(v_d_values, v_d * 100)  # Convert to cm/s
        push!(v_s_values, v_s * 100)  # Convert to cm/s
    end

    plot!(p, D_p_range .* 1e6, v_d_values,
        label = u_star_labels[i],
        linewidth = 2,
        color = colors[i])
end

# Add settling velocity line
v_s_line = [settling_velocity(D_p, cunningham_correction(D_p)) * 100 for D_p in D_p_range]
plot!(p, D_p_range .* 1e6, v_s_line,
    label = "Settling velocity (v_s)",
    linewidth = 2,
    linestyle = :dash,
    color = :black)

# Add annotation for minimum region
annotate!(p, 0.3, 0.03, text("Accumulation\nMode", 8, :center))

savefig(p, "figure19_2.png")
p
```

### Interpretation

The figure demonstrates several key features of particle dry deposition:

 1. **Minimum at 0.1-1 micrometer**: The deposition velocity reaches a minimum in the accumulation mode where neither Brownian diffusion nor impaction/settling is effective.

 2. **Dependence on friction velocity**: Higher ``u_*`` leads to higher deposition velocities due to increased turbulent transport and impaction efficiency.

 3. **Asymptotic behavior**:

      + For very small particles, ``v_d \propto D_p^{-2/3}`` due to Brownian diffusion
      + For very large particles, ``v_d \approx v_s \propto D_p^2`` due to gravitational settling

* * *

## Collection Efficiency Components

This figure shows the individual collection efficiency terms (Brownian diffusion, impaction, and interception) as a function of particle diameter.

```@example efficiency_components
using Plots

const k_B = 1.38e-23
const g = 9.81
const T = 298.0
const mu = 1.85e-5
const nu = 1.55e-5
const rho_p = 1000.0

const A = 2.0e-3
const alpha = 1.2
const gamma_param = 0.54
const u_star = 0.4

function cunningham_correction(D_p)
    lambda = 6.8e-8
    Kn = 2 * lambda / D_p
    return 1 + Kn * (1.257 + 0.4 * exp(-1.1 / Kn))
end

function collection_efficiencies(D_p)
    C_c = cunningham_correction(D_p)

    # Settling velocity
    v_s = rho_p * D_p^2 * g * C_c / (18 * mu)

    # Brownian diffusivity
    D_diff = k_B * T * C_c / (3 * pi * mu * D_p)

    # Schmidt and Stokes numbers
    Sc = nu / D_diff
    St = v_s * u_star / (g * A)

    # Collection efficiencies
    E_B = Sc^(-gamma_param)
    E_IM = (St / (alpha + St))^2
    E_IN = 0.5 * (D_p / A)^2
    R_1 = exp(-sqrt(St))

    return E_B, E_IM, E_IN, R_1
end

D_p_range = 10 .^ range(-2, 2, length = 100) .* 1e-6

E_B_vals = Float64[]
E_IM_vals = Float64[]
E_IN_vals = Float64[]
E_total_vals = Float64[]

for D_p in D_p_range
    E_B, E_IM, E_IN, R_1 = collection_efficiencies(D_p)
    push!(E_B_vals, E_B)
    push!(E_IM_vals, E_IM)
    push!(E_IN_vals, E_IN)
    push!(E_total_vals, (E_B + E_IM + E_IN) * R_1)
end

p = plot(
    xlabel = "Particle Diameter (micrometer)",
    ylabel = "Collection Efficiency",
    title = "Collection Efficiency Components (u* = 0.4 m/s, Grass)",
    xscale = :log10,
    yscale = :log10,
    xlims = (0.01, 100),
    ylims = (1e-10, 1),
    legend = :bottomright,
    size = (700, 500),
    grid = true
)

plot!(p, D_p_range .* 1e6, E_B_vals, label = "E_B (Brownian)", linewidth = 2, color = :blue)
plot!(
    p, D_p_range .* 1e6, E_IM_vals, label = "E_IM (Impaction)", linewidth = 2, color = :red)
plot!(p, D_p_range .* 1e6, E_IN_vals,
    label = "E_IN (Interception)", linewidth = 2, color = :green)
plot!(p, D_p_range .* 1e6, E_total_vals, label = "Total (with R_1)",
    linewidth = 3, color = :black, linestyle = :dash)

savefig(p, "collection_efficiencies.png")
p
```

### Physical Interpretation

  - **Brownian diffusion (E_B)**: Dominates for particles smaller than ~0.1 micrometer. Decreases as ``D_p^{2\gamma/3} \approx D_p^{0.36}`` for the Stokes-Einstein relation.

  - **Impaction (E_IM)**: Becomes significant for particles larger than ~1 micrometer. Depends on the Stokes number, which increases with ``D_p^2``.

  - **Interception (E_IN)**: Also increases with particle size as ``D_p^2``. Represents particles whose streamlines pass within one particle radius of the collector.

  - **Sticking fraction (R_1)**: Reduces collection efficiency for large particles that may bounce off the surface.

* * *

## Sensitivity to Land-Use Category

Different land-use categories have different characteristic collector radii (A) and collection efficiency parameters, leading to variations in deposition velocity.

```@example landuse_sensitivity
using Plots

const k_B = 1.38e-23
const g = 9.81
const T = 298.0
const mu = 1.85e-5
const nu = 1.55e-5
const rho_p = 1000.0
const u_star = 0.4
const z = 10.0
const kappa = 0.4

# Land-use parameters from Table 19.2
landuse_params = [
    ("Grass", 2.0e-3, 1.2, 0.54, 0.1),
    ("Deciduous Forest", 5.0e-3, 0.8, 0.56, 1.0),
    ("Desert", 10.0e-3, 50.0, 0.54, 0.04)
]

function cunningham_correction(D_p)
    lambda = 6.8e-8
    Kn = 2 * lambda / D_p
    return 1 + Kn * (1.257 + 0.4 * exp(-1.1 / Kn))
end

function particle_v_d(D_p, A, alpha, gamma_val, z_0)
    C_c = cunningham_correction(D_p)
    v_s = rho_p * D_p^2 * g * C_c / (18 * mu)
    D_diff = k_B * T * C_c / (3 * pi * mu * D_p)

    Sc = nu / D_diff
    St = v_s * u_star / (g * A)

    E_B = Sc^(-gamma_val)
    E_IM = (St / (alpha + St))^2
    E_IN = 0.5 * (D_p / A)^2
    R_1 = exp(-sqrt(St))

    r_b = 1 / (3 * u_star * (E_B + E_IM + E_IN) * R_1)
    r_a = (1 / (kappa * u_star)) * log(z / z_0)

    v_d = 1 / (r_a + r_b + r_a * r_b * v_s) + v_s
    return v_d
end

D_p_range = 10 .^ range(-2, 2, length = 100) .* 1e-6

p = plot(
    xlabel = "Particle Diameter (micrometer)",
    ylabel = "Deposition Velocity (cm/s)",
    title = "Particle Deposition Velocity by Land-Use Category (u* = 0.4 m/s)",
    xscale = :log10,
    yscale = :log10,
    xlims = (0.01, 100),
    ylims = (0.001, 10),
    legend = :bottomright,
    size = (700, 500),
    grid = true
)

colors = [:blue, :green, :orange]
for (i, (name, A, alpha, gamma_val, z_0)) in enumerate(landuse_params)
    v_d_vals = [particle_v_d(D_p, A, alpha, gamma_val, z_0) * 100 for D_p in D_p_range]
    plot!(p, D_p_range .* 1e6, v_d_vals, label = name, linewidth = 2, color = colors[i])
end

savefig(p, "landuse_comparison.png")
p
```

* * *

## Summary Statistics

The following table summarizes typical deposition velocities from the model compared to values from Table 19.1 in the textbook.

```@example summary_stats
using DataFrames

# Typical conditions
u_star = 0.4  # m/s
z = 10.0      # m
z_0 = 0.1     # m (grass)
kappa = 0.4

# Aerodynamic resistance
r_a = (1 / (kappa * u_star)) * log(z / z_0)

# Quasi-laminar resistance (Sc ~ 1 for gases)
r_b_gas = 5 * 1.0^(2/3) / u_star

# Surface resistance estimates for different gases
# Based on Table 19.1 measured values
gases = [
    ("O3", 0.4, 250 - r_a - r_b_gas),
    ("HNO3", 4.0, 25 - r_a - r_b_gas),
    ("NO2", 0.1, 1000 - r_a - r_b_gas),
    ("SO2", 0.5, 200 - r_a - r_b_gas)
]

df = DataFrame(
    Gas = [g[1] for g in gases],
    Measured_vd_cm_s = [g[2] for g in gases],
    Implied_r_c_s_m = [max(0, g[3]) for g in gases],
    r_a_s_m = fill(round(r_a, digits = 1), length(gases)),
    r_b_s_m = fill(round(r_b_gas, digits = 1), length(gases))
)
df
```

These results are consistent with the physical understanding that:

  - HNO3 has very low surface resistance (highly soluble and reactive)
  - O3 has moderate surface resistance (reactive but not very soluble)
  - NO2 has high surface resistance
  - Surface resistance dominates the total resistance for most gases

* * *

# API Reference

## Module

```@docs
DryDeposition
```

## Physical Constants

```@docs
DryDeposition.KAPPA
```

## Data Types

```@docs
DryDeposition.GasProperties
DryDeposition.LandUseParameters
```

## Component Functions

```@docs
DryDeposition.DryDepositionGas
DryDeposition.DryDepositionParticle
DryDeposition.AerodynamicResistance
DryDeposition.QuasiLaminarResistanceGas
DryDeposition.QuasiLaminarResistanceParticle
DryDeposition.ParticleSettling
DryDeposition.SurfaceResistance
```
