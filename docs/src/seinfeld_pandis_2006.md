# Wet Deposition (Seinfeld & Pandis, 2006)

## Overview

This module implements the wet deposition equations from Chapter 20 of Seinfeld and Pandis (2006), covering below-cloud scavenging of gases and particles, and associated scavenging coefficients.

Wet deposition is a major pathway for removing trace gases and aerosol particles from the atmosphere. The implemented equations describe:

- **Gas-phase mass transfer** to falling raindrops via the Sherwood number correlation
- **Irreversible gas scavenging** for highly soluble species (e.g., HNO₃)
- **Reversible gas scavenging** retaining equilibrium (Henry's law) effects
- **Particle scavenging** via the Slinn (1983) semi-empirical collision efficiency
- **Wet deposition flux** and deposition velocity

**Reference**: Seinfeld, J. H. and Pandis, S. N. (2006). *Atmospheric Chemistry and Physics: From Air Pollution to Climate Change*, 2nd ed., Chapter 20: Wet Deposition, pp. 932–979. John Wiley & Sons, Inc.

```@docs
MassTransferCoeff
GasScavengingCoeff
ReversibleGasScavenging
BelowCloudGasScavenging
ParticleCollectionEfficiency
ParticleScavengingCoeff
WetDepositionFlux
```

## Implementation

### Mass Transfer Coefficient (Eq. 20.12)

The mass transfer coefficient ``K_c`` for gas absorption by a spherical raindrop is computed from the Sherwood number correlation:

```math
K_c = \frac{D_g}{D_p}\left[2 + 0.6\,\text{Re}^{1/2}\,\text{Sc}^{1/3}\right]
```

where ``\text{Re} = \rho_\text{air} U_t D_p / \mu_\text{air}`` and ``\text{Sc} = \mu_\text{air} / (\rho_\text{air} D_g)``.

```@example sp2006
using AtmosphericDeposition
using ModelingToolkit
using DataFrames, Symbolics, DynamicQuantities

sys = MassTransferCoeff()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Gas Scavenging Coefficient (Eq. 20.25)

For irreversible uptake by monodisperse drops:

```math
\Lambda = \frac{6\,p_0\,K_c}{D_p\,U_t}
```

```@example sp2006
sys = GasScavengingCoeff()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Particle Collection Efficiency (Eq. 20.53)

The Slinn (1983) semi-empirical correlation computes the collection efficiency ``E`` as a sum of three terms: Brownian diffusion, interception, and inertial impaction.

```@example sp2006
sys = ParticleCollectionEfficiency()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Equations

```@example sp2006
sys = BelowCloudGasScavenging()
equations(sys)
```

## Analysis

### Table 20.1: Estimation of the Scavenging Coefficient Λ for Irreversible Scavenging

**Table 20.1** from Seinfeld & Pandis (2006, p. 941) shows the scavenging coefficient Λ for irreversible gas scavenging in a homogeneous atmosphere with precipitation rate ``p_0 = 1`` mm h⁻¹.

| ``D_p`` (cm) | ``U_t`` (cm s⁻¹) | ``K_c`` (cm s⁻¹) | ``Λ`` (h⁻¹) | ``1/Λ`` (h) |
|:-------------|:-----------------|:-----------------|:------------|:------------|
| 0.001        | 0.3              | 220              | 4.4 × 10⁵   | 2.3 × 10⁻⁶  |
| 0.01         | 26               | 32               | 73.8        | 0.01        |
| 0.1          | 300              | 13               | 0.26        | 3.8         |
| 1.0          | 1000             | 6                | 0.0036      | 278         |

The table demonstrates that the scavenging coefficient depends dramatically on raindrop diameter. Very small drops are very efficient in scavenging soluble gases for two reasons: (1) they fall more slowly, so they have more time in their transit to "clean" the atmosphere; and (2) mass transfer is more efficient for these drops (high ``K_c``). Note that the scavenging rate varies over eight orders of magnitude when the drop diameter increases by three orders. Also note that the scavenging time scale (``1/Λ``) can vary from less than a second to several hours depending on the raindrop size distribution.

### Gas Scavenging Coefficient vs. Raindrop Diameter (Table 20.1)

Validation against Table 20.1 in Seinfeld & Pandis (2006). The scavenging coefficient ``\Lambda`` is shown as a function of raindrop diameter for ``p_0 = 1`` mm/h.

```@example sp2006
using AtmosphericDeposition: mass_transfer_coeff, gas_scavenging_coeff,
    μ_air_sp, ρ_air_sp
using ModelingToolkit
using Plots
using DynamicQuantities
using Symbolics

@parameters D_g_val = 1.26e-5 [unit = u"m^2/s"]
@parameters D_p_val [unit = u"m"]
@parameters U_t_val [unit = u"m/s"]
@parameters p₀_val = 2.778e-7 [unit = u"m/s"]

defaults = [μ_air_sp => 1.72e-5, ρ_air_sp => 1.225]

# Table 20.1 data (Seinfeld & Pandis 2006, p. 941):
# D_p (m), U_t (m/s), K_c from table (m/s), Λ from table (h⁻¹)
# Note: Original table uses cm and cm/s; converted to SI units here
table_data = [
    (1e-5,   0.003,  2.20,  4.4e5),   # D_p = 0.001 cm, U_t = 0.3 cm/s, K_c = 220 cm/s
    (1e-4,   0.26,   0.32,  73.8),    # D_p = 0.01 cm, U_t = 26 cm/s, K_c = 32 cm/s
    (1e-3,   3.00,   0.13,  0.26),    # D_p = 0.1 cm, U_t = 300 cm/s, K_c = 13 cm/s
    (1e-2,  10.00,   0.06,  0.0036),  # D_p = 1.0 cm, U_t = 1000 cm/s, K_c = 6 cm/s
]

Kc_expr = mass_transfer_coeff(D_g_val, D_p_val, U_t_val)
Λ_expr = gas_scavenging_coeff(Kc_expr, U_t_val, D_p_val, p₀_val)

D_p_vals = Float64[]
Λ_computed = Float64[]
Λ_table = Float64[]

for (dp, ut, kc_exp, lam_exp) in table_data
    Λ_val = Symbolics.value(substitute(Λ_expr,
        Dict(D_g_val => 1.26e-5, D_p_val => dp, U_t_val => ut,
            p₀_val => 2.778e-7, defaults...)))
    push!(D_p_vals, dp * 100) # convert to cm
    push!(Λ_computed, Float64(Λ_val) * 3600) # convert s⁻¹ → h⁻¹
    push!(Λ_table, lam_exp)
end

plot(D_p_vals, Λ_computed, yscale=:log10, xscale=:log10,
    marker=:circle, label="Computed",
    xlabel="Drop diameter (cm)", ylabel="Λ (h⁻¹)",
    title="Gas Scavenging Coefficient (Table 20.1)")
scatter!(D_p_vals, Λ_table, marker=:square, label="Table 20.1")
```

### Below-Cloud Gas Concentration Decay (Eq. 20.24)

Exponential decay of gas-phase concentration due to irreversible scavenging at different scavenging coefficients.

```@example sp2006
using OrdinaryDiffEqDefault

sys = BelowCloudGasScavenging()
compiled = mtkcompile(sys)

C_g0 = 1e-5
tspan = (0.0, 3600.0 * 5) # 5 hours

Λ_values = [0.26 / 3600, 1.0 / 3600, 3.3 / 3600] # h⁻¹ → s⁻¹
labels = ["Λ = 0.26 h⁻¹" "Λ = 1.0 h⁻¹" "Λ = 3.3 h⁻¹"]

p = plot(xlabel="Time (hours)", ylabel="C_g / C_g₀",
    title="Below-Cloud Gas Scavenging (Eq. 20.24)")

for (i, Λ_val) in enumerate(Λ_values)
    # Back-calculate p₀ from Λ = 6*p₀*K_c/(D_p*U_t) with K_c=0.13, D_p=1e-3, U_t=3
    p₀_val = Λ_val * 1e-3 * 3.0 / (6 * 0.13)
    prob = ODEProblem(compiled,
        [compiled.C_g => C_g0],
        tspan,
        [compiled.K_c => 0.13, compiled.U_t => 3.0, compiled.D_p => 1e-3,
            compiled.h => 1000.0, compiled.p₀ => p₀_val])
    sol = solve(prob)
    plot!(p, sol.t ./ 3600, sol[compiled.C_g] ./ C_g0,
        label=labels[i])
end
p
```

### Collection Efficiency vs. Particle Size (Figure 20.6)

The Slinn (1983) collision efficiency as a function of collected particle radius, showing the Brownian diffusion regime at small sizes, the Greenfield gap near 1 μm, and the inertial impaction regime at large sizes.

```@example sp2006
using AtmosphericDeposition: slinn_collection_efficiency
using DynamicQuantities
using Symbolics

@parameters D_p_rain [unit = u"m"]
@parameters U_t_rain [unit = u"m/s"]
@parameters d_p_aer [unit = u"m"]
@parameters ρ_p_aer [unit = u"kg/m^3"]
@parameters T_val [unit = u"K"]

E_expr = slinn_collection_efficiency(D_p_rain, U_t_rain, d_p_aer, ρ_p_aer, T_val)

defaults_E = [μ_air_sp => 1.72e-5, ρ_air_sp => 1.225, AtmosphericDeposition.μ_w_sp => 1e-3,
    AtmosphericDeposition.ρ_w_sp => 1000.0, AtmosphericDeposition.g_sp => 9.81,
    AtmosphericDeposition.kB_sp => 1.3806488e-23]

# Two collector drop sizes: radius = 0.1 mm (D_p = 0.2 mm) and 1 mm (D_p = 2 mm)
drop_configs = [
    (2e-4, 0.8, "r = 0.1 mm"),    # D_p = 0.2 mm, U_t ≈ 0.8 m/s
    (2e-3, 6.5, "r = 1.0 mm"),    # D_p = 2 mm, U_t ≈ 6.5 m/s
]

# Particle radii from 0.001 to 10 μm
particle_radii = 10 .^ range(-3, 1, length=50) .* 1e-6  # in meters (radii)
particle_diameters = 2 .* particle_radii

p = plot(xscale=:log10, yscale=:log10,
    xlabel="Particle radius (μm)", ylabel="Collection efficiency E",
    title="Slinn Collision Efficiency (Fig. 20.6)",
    ylim=(1e-5, 1.0))

for (Dp, Ut, lbl) in drop_configs
    E_vals = Float64[]
    for dp in particle_diameters
        E_val = Symbolics.value(substitute(E_expr,
            Dict(D_p_rain => Dp, U_t_rain => Ut, d_p_aer => dp,
                ρ_p_aer => 1000.0, T_val => 298.0, defaults_E...)))
        push!(E_vals, max(Float64(E_val), 1e-10))
    end
    plot!(p, particle_radii .* 1e6, E_vals, label=lbl)
end
p
```
