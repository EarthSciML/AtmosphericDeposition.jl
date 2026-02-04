export MassTransferCoeff, GasScavengingCoeff, ReversibleGasScavenging,
    ParticleCollectionEfficiency, ParticleScavengingCoeff,
    WetDepositionFlux, BelowCloudGasScavenging

# ============================================================================
# Chapter 20: Wet Deposition
# Reference: Seinfeld, J.H. and Pandis, S.N. (2006). Atmospheric Chemistry
#   and Physics: From Air Pollution to Climate Change, 2nd ed., Chapter 20:
#   Wet Deposition, pp. 932–979. John Wiley & Sons, Inc.
# ============================================================================

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
@constants μ_air_sp = 1.72e-5 [unit = u"kg/m/s", description = "Dynamic viscosity of air (Table A.7, Seinfeld & Pandis 2006)"]
@constants ρ_air_sp = 1.225 [unit = u"kg/m^3", description = "Air density at standard conditions (Seinfeld & Pandis 2006)"]
@constants μ_w_sp = 1.0e-3 [unit = u"kg/m/s", description = "Dynamic viscosity of water"]
@constants ρ_w_sp = 1000.0 [unit = u"kg/m^3", description = "Water density"]
@constants g_sp = 9.81 [unit = u"m/s^2", description = "Gravitational acceleration"]
@constants kB_sp = 1.3806488e-23 [unit = u"kg*m^2/s^2/K", description = "Boltzmann constant"]

# ---------------------------------------------------------------------------
# Eq. 20.12 — Mass Transfer Coefficient (Sherwood number correlation)
# K_c = (D_g / D_p) * [2 + 0.6 * Re^(1/2) * Sc^(1/3)]
# where Re = ρ_air * U_t * D_p / μ_air, Sc = μ_air / (ρ_air * D_g)
# ---------------------------------------------------------------------------
"""
    mass_transfer_coeff(D_g, D_p, U_t)

Compute the gas-phase mass transfer coefficient to a spherical drop using
the Sherwood number correlation (Eq. 20.12 in Seinfeld & Pandis, 2006).

# Arguments
- `D_g`: gas-phase diffusivity [m²/s]
- `D_p`: raindrop diameter [m]
- `U_t`: terminal velocity of the raindrop [m/s]
"""
function mass_transfer_coeff(D_g, D_p, U_t)
    Re = ρ_air_sp * U_t * D_p / μ_air_sp         # Reynolds number (dimensionless)
    Sc = μ_air_sp / (ρ_air_sp * D_g)              # Schmidt number (dimensionless)
    # Eq. 20.12: Sh = 2 + 0.6 * Re^(1/2) * Sc^(1/3)
    Sh = 2 + 0.6 * sqrt(Re) * cbrt(Sc)
    K_c = D_g / D_p * Sh
    return K_c
end

"""
    MassTransferCoeff(; name = :MassTransferCoeff)

ModelingToolkit component for the gas-phase mass transfer coefficient to a
raindrop (Eq. 20.12, Seinfeld & Pandis, 2006).
"""
@component function MassTransferCoeff(; name = :MassTransferCoeff)
    params = @parameters begin
        D_g = 1.5e-5, [unit = u"m^2/s", description = "Gas-phase diffusivity"]
        D_p = 1.0e-3, [unit = u"m", description = "Raindrop diameter"]
        U_t = 4.0, [unit = u"m/s", description = "Terminal velocity of raindrop"]
    end

    vars = @variables begin
        K_c(t), [unit = u"m/s", description = "Mass transfer coefficient (Eq. 20.12)"]
    end

    eqs = [
        K_c ~ mass_transfer_coeff(D_g, D_p, U_t), # Eq. 20.12
    ]

    return System(eqs, t, vars, params; name)
end

# ---------------------------------------------------------------------------
# Eq. 20.25 — Gas Scavenging Coefficient (Irreversible, monodisperse drops)
# Λ = 6e-3 * p₀ * K_c / (U_t * D_p)
# where p₀ is in mm/h; Λ in s⁻¹
#
# Derivation of the 6e-3 factor:
# From Eq. 20.22, F_bc = N_D * π * D_p² * K_c * h * C_g
# where N_D = 6 * p₀ / (π * D_p³ * U_t) from Eq. 20.56.
# Then Λ = F_bc / (h * C_g) = 6 * p₀ * K_c / (D_p * U_t)
# When p₀ is in mm/h, multiply by (1 mm/h → m/s) = 1/(3600*1000) = 1/3.6e6
# So Λ [s⁻¹] = 6 / 3.6e6 * p₀[mm/h] * K_c[m/s] / (D_p[m] * U_t[m/s])
#            ≈ 1.667e-6 * p₀ * K_c / (D_p * U_t)
#
# However, the text uses the convention Λ [h⁻¹] with a different prefactor.
# We implement the version that directly yields s⁻¹.
# ---------------------------------------------------------------------------
"""
    gas_scavenging_coeff(K_c, U_t, D_p, p₀_SI)

Compute the below-cloud gas scavenging coefficient for irreversible uptake
by monodisperse raindrops (Eq. 20.25 in Seinfeld & Pandis, 2006).

All arguments in SI units. `p₀_SI` is precipitation rate in m/s.
"""
function gas_scavenging_coeff(K_c, U_t, D_p, p₀_SI)
    # From Eq. 20.56: N_D = 6 * p₀ / (π * D_p³ * U_t)
    # Λ = (π/4) * D_p² * U_t * N_D * (4 * K_c / U_t)  -- simplifies to:
    # Λ = 6 * p₀ * K_c / (D_p * U_t)   [all SI]
    Λ = 6 * p₀_SI * K_c / (D_p * U_t)
    return Λ
end

"""
    GasScavengingCoeff(; name = :GasScavengingCoeff)

ModelingToolkit component for the irreversible gas scavenging coefficient
(Eq. 20.25, Seinfeld & Pandis, 2006).
"""
@component function GasScavengingCoeff(; name = :GasScavengingCoeff)
    params = @parameters begin
        K_c = 0.13, [unit = u"m/s", description = "Mass transfer coefficient"]
        U_t = 3.0, [unit = u"m/s", description = "Terminal velocity of raindrop"]
        D_p = 1.0e-3, [unit = u"m", description = "Raindrop diameter"]
        p₀ = 2.778e-7, [unit = u"m/s", description = "Precipitation rate (1 mm/h = 2.778e-7 m/s)"]
    end

    vars = @variables begin
        Λ_gas(t), [unit = u"s^-1", description = "Gas scavenging coefficient (Eq. 20.25)"]
    end

    eqs = [
        Λ_gas ~ gas_scavenging_coeff(K_c, U_t, D_p, p₀), # Eq. 20.25
    ]

    return System(eqs, t, vars, params; name)
end

# ---------------------------------------------------------------------------
# Eq. 20.28, 20.33, 20.35 — Reversible Gas Scavenging
#
# The key quantity is H* R T where:
#   H* = effective Henry's law coefficient [M/atm]
#   R = gas constant = 8.20578e-5 m³·atm/(mol·K)
#   T = temperature [K]
#
# The product H*RT has units of (mol/L/atm)(m³·atm/mol/K)(K) = m³/L = 1e-3
# (for H*=1). For dimensional consistency when C_g is in mol/m³ and C_aq
# in mol/m³, we need HRT_dimless = H* × R × T × 1000 (converting L→m³).
#
# We take HRT as a dimensionless parameter that the user computes externally.
# ---------------------------------------------------------------------------
"""
    reversible_drop_conc(C_g, HRT, K_c, D_p, U_t, z, C_aq0)

Compute the aqueous-phase concentration in a falling raindrop at fall
distance `z` from cloud base for a reversibly-absorbed gas species
(Eq. 20.28, Seinfeld & Pandis, 2006).

`HRT` = H* × R × T × 1000 (dimensionless), representing the ratio
C_aq_eq / C_g when both are in mol/m³.
"""
function reversible_drop_conc(C_g, HRT, K_c, D_p, U_t, z, C_aq0)
    C_eq = C_g * HRT                                   # equilibrium aqueous concentration [mol/m³]
    α = 6 * K_c * z / (D_p * U_t * HRT)                # dimensionless exponent (Eq. 20.28)
    C_aq = C_eq - (C_eq - C_aq0) * exp(-α)             # Eq. 20.28
    return C_aq
end

"""
    reversible_scavenging_flux(C_g, HRT, K_c, D_p, U_t, h, C_aq0, p₀_SI)

Compute the total below-cloud scavenging rate for a reversibly-absorbed gas
(Eq. 20.35, Seinfeld & Pandis, 2006). All SI units.

`HRT` = H* × R × T × 1000 (dimensionless).
`p₀_SI` = precipitation rate [m/s].
"""
function reversible_scavenging_flux(C_g, HRT, K_c, D_p, U_t, h, C_aq0, p₀_SI)
    α = 6 * K_c * h / (D_p * U_t * HRT)
    # Eq. 20.35: F_bc = p₀ * [C_g * HRT - C_aq0] * [1 - exp(-α)]
    # The factor converts from droplet volume flux to mol/m²/s
    F_bc = p₀_SI * (C_g * HRT - C_aq0) * (1 - exp(-α))
    return F_bc
end

"""
    ReversibleGasScavenging(; name = :ReversibleGasScavenging)

ModelingToolkit component for reversible below-cloud gas scavenging
(Eqs. 20.28, 20.33, 20.35, Seinfeld & Pandis, 2006).

The parameter `HRT` is the dimensionless product H* × R × T × 1000
(effective Henry's law coefficient × gas constant × temperature × unit conversion).
For HNO₃ at 298 K: HRT = 2.1e5 × 8.206e-5 × 298 × 1000 ≈ 5.14e6.
"""
@component function ReversibleGasScavenging(; name = :ReversibleGasScavenging)
    params = @parameters begin
        C_g = 1.0e-6, [unit = u"mol/m^3", description = "Gas-phase concentration"]
        HRT = 5.14e6, [unit = u"1", description = "H* × R × T × 1000 (dimensionless)"]
        K_c = 0.13, [unit = u"m/s", description = "Mass transfer coefficient"]
        D_p = 1.0e-3, [unit = u"m", description = "Raindrop diameter"]
        U_t = 4.0, [unit = u"m/s", description = "Terminal velocity of raindrop"]
        h = 1000.0, [unit = u"m", description = "Cloud base height / fall distance"]
        C_aq0 = 0.0, [unit = u"mol/m^3", description = "Initial aqueous-phase concentration"]
        p₀ = 2.778e-7, [unit = u"m/s", description = "Precipitation rate"]
    end

    vars = @variables begin
        C_aq_ground(t), [unit = u"mol/m^3", description = "Ground-level drop concentration (Eq. 20.33)"]
        F_bc(t), [unit = u"mol/m^2/s", description = "Below-cloud scavenging flux (Eq. 20.35)"]
    end

    eqs = [
        C_aq_ground ~ reversible_drop_conc(C_g, HRT, K_c, D_p, U_t, h, C_aq0), # Eq. 20.33
        F_bc ~ reversible_scavenging_flux(C_g, HRT, K_c, D_p, U_t, h, C_aq0, p₀), # Eq. 20.35
    ]

    return System(eqs, t, vars, params; name)
end

# ---------------------------------------------------------------------------
# Eq. 20.22, 20.24 — Below-Cloud Gas Scavenging (Irreversible)
# F_bc = Λ * h * C_g   (Eq. 20.22)
# dC_g/dt = -Λ * C_g   (Eq. 20.24)
# ---------------------------------------------------------------------------
"""
    BelowCloudGasScavenging(; name = :BelowCloudGasScavenging)

ModelingToolkit component for irreversible below-cloud gas scavenging.
Computes the scavenging coefficient (Eq. 20.25), the below-cloud flux
(Eq. 20.22), and the exponential gas concentration decay (Eq. 20.24).
"""
@component function BelowCloudGasScavenging(; name = :BelowCloudGasScavenging)
    params = @parameters begin
        K_c = 0.13, [unit = u"m/s", description = "Mass transfer coefficient"]
        U_t = 3.0, [unit = u"m/s", description = "Terminal velocity of raindrop"]
        D_p = 1.0e-3, [unit = u"m", description = "Raindrop diameter"]
        h = 1000.0, [unit = u"m", description = "Cloud base height"]
        p₀ = 2.778e-7, [unit = u"m/s", description = "Precipitation rate"]
    end

    vars = @variables begin
        C_g(t), [unit = u"mol/m^3", description = "Gas-phase concentration below cloud"]
        Λ_gas(t), [unit = u"s^-1", description = "Gas scavenging coefficient (Eq. 20.25)"]
        F_bc(t), [unit = u"mol/m^2/s", description = "Below-cloud scavenging flux (Eq. 20.22)"]
    end

    eqs = [
        Λ_gas ~ gas_scavenging_coeff(K_c, U_t, D_p, p₀), # Eq. 20.25
        D(C_g) ~ -Λ_gas * C_g,                             # Eq. 20.24
        F_bc ~ Λ_gas * h * C_g,                             # Eq. 20.22
    ]

    return System(eqs, t, vars, params; name)
end

# ---------------------------------------------------------------------------
# Eq. 20.53, 20.54 — Slinn (1983) Collision Efficiency for Particle Scavenging
# ---------------------------------------------------------------------------
"""
    particle_relaxation_time(d_p, ρ_p)

Compute the particle relaxation time τ = ρ_p * d_p² / (18 * μ_air)
for Stokes drag regime.
"""
function particle_relaxation_time(d_p, ρ_p)
    τ = ρ_p * d_p^2 / (18 * μ_air_sp)
    return τ
end

"""
    particle_terminal_velocity(d_p, ρ_p)

Compute the Stokes settling velocity u_t = τ * g for small particles.
"""
function particle_terminal_velocity(d_p, ρ_p)
    τ = particle_relaxation_time(d_p, ρ_p)
    u_t = τ * g_sp
    return u_t
end

"""
    particle_diffusivity(d_p, T)

Compute Brownian diffusivity of an aerosol particle via Stokes-Einstein:
D_diff = k_B * T / (3π * μ_air * d_p).
"""
function particle_diffusivity(d_p, T)
    D_diff = kB_sp * T / (3 * pi * μ_air_sp * d_p)
    return D_diff
end

"""
    slinn_collection_efficiency(D_p, U_t, d_p, ρ_p, T)

Compute the semi-empirical collision (collection) efficiency E between a
raindrop of diameter `D_p` and an aerosol particle of diameter `d_p`
(Eq. 20.53–20.54, Seinfeld & Pandis, 2006).

The three terms represent:
1. Brownian diffusion
2. Interception
3. Inertial impaction (Stokes number term)
"""
function slinn_collection_efficiency(D_p, U_t, d_p, ρ_p, T)
    # Reynolds number based on drop radius: Re = D_p * U_t * ρ_air / (2 * μ_air)
    Re = D_p * U_t * ρ_air_sp / (2 * μ_air_sp)

    # Schmidt number for particle: Sc = μ_air / (ρ_air * D_diff)
    D_diff = particle_diffusivity(d_p, T)
    Sc = μ_air_sp / (ρ_air_sp * D_diff)

    # Stokes number: St = 2 * τ * (U_t - u_t) / D_p
    τ = particle_relaxation_time(d_p, ρ_p)
    u_t_particle = particle_terminal_velocity(d_p, ρ_p)
    St = 2 * τ * (U_t - u_t_particle) / D_p

    # Diameter ratio: φ = d_p / D_p
    φ = d_p / D_p

    # Viscosity ratio: ω = μ_w / μ_air
    ω = μ_w_sp / μ_air_sp

    # Critical Stokes number (Eq. 20.54)
    S_star = (1.2 + (1 // 12) * log(1 + Re)) / (1 + log(1 + Re))

    # Term 1: Brownian diffusion (Eq. 20.53, first line)
    E_diff = (4 / (Re * Sc)) * (1 + 0.4 * sqrt(Re) * cbrt(Sc) + 0.16 * sqrt(Re) * sqrt(Sc))

    # Term 2: Interception (Eq. 20.53, second line)
    E_int = 4 * φ * (1 / ω + (1 + 2 * sqrt(Re)) * φ)

    # Term 3: Inertial impaction (Eq. 20.53, third line)
    # Only active if St - S_star >= 2/3, else zero
    # Use max(0, ...) to avoid negative base in power (both branches are evaluated symbolically)
    # Note: For particles of density different from 1 g/cm³ (= 1000 kg/m³), the impaction term
    # should be scaled by (ρ_w/ρ_p)^(1/2) per note on p. 950
    stokes_diff = St - S_star
    safe_ratio = max(0, stokes_diff) / (max(0, stokes_diff) + 2 // 3)
    density_scaling = sqrt(ρ_w_sp / ρ_p)
    E_imp_unscaled = ifelse(stokes_diff < 2 // 3, zero(safe_ratio),
        safe_ratio * sqrt(safe_ratio))
    E_imp = E_imp_unscaled * density_scaling

    E_total = E_diff + E_int + E_imp
    # Clamp to maximum of 1
    E = ifelse(E_total > 1, one(E_total), E_total)
    return E
end

"""
    ParticleCollectionEfficiency(; name = :ParticleCollectionEfficiency)

ModelingToolkit component for the Slinn (1983) semi-empirical particle-drop
collision (collection) efficiency (Eq. 20.53–20.54, Seinfeld & Pandis, 2006).
"""
@component function ParticleCollectionEfficiency(; name = :ParticleCollectionEfficiency)
    params = @parameters begin
        D_p = 1.0e-3, [unit = u"m", description = "Raindrop diameter"]
        U_t = 4.0, [unit = u"m/s", description = "Terminal velocity of raindrop"]
        d_p = 1.0e-6, [unit = u"m", description = "Aerosol particle diameter"]
        ρ_p = 1000.0, [unit = u"kg/m^3", description = "Aerosol particle density"]
        T = 298.0, [unit = u"K", description = "Temperature"]
    end

    vars = @variables begin
        E(t), [unit = u"1", description = "Collection efficiency (Eq. 20.53) (dimensionless)"]
    end

    eqs = [
        E ~ slinn_collection_efficiency(D_p, U_t, d_p, ρ_p, T), # Eq. 20.53
    ]

    return System(eqs, t, vars, params; name)
end

# ---------------------------------------------------------------------------
# Eq. 20.57 — Particle Scavenging Coefficient (monodisperse drops)
# Λ(d_p) = (3/2) * E(D_p, d_p) * p₀ / D_p
# ---------------------------------------------------------------------------
"""
    particle_scavenging_coeff(E, p₀_SI, D_p)

Compute the particle scavenging coefficient for monodisperse raindrops
(Eq. 20.57, Seinfeld & Pandis, 2006).

# Arguments
- `E`: collection efficiency (dimensionless)
- `p₀_SI`: precipitation rate [m/s]
- `D_p`: raindrop diameter [m]
"""
function particle_scavenging_coeff(E, p₀_SI, D_p)
    # Eq. 20.57: Λ = (3/2) * E * p₀ / D_p
    Λ = (3 // 2) * E * p₀_SI / D_p
    return Λ
end

"""
    ParticleScavengingCoeff(; name = :ParticleScavengingCoeff)

ModelingToolkit component for the particle scavenging coefficient
with monodisperse raindrops (Eq. 20.57, Seinfeld & Pandis, 2006).
"""
@component function ParticleScavengingCoeff(; name = :ParticleScavengingCoeff)
    params = @parameters begin
        D_p = 1.0e-3, [unit = u"m", description = "Raindrop diameter"]
        U_t = 4.0, [unit = u"m/s", description = "Terminal velocity of raindrop"]
        d_p = 1.0e-6, [unit = u"m", description = "Aerosol particle diameter"]
        ρ_p = 1000.0, [unit = u"kg/m^3", description = "Aerosol particle density"]
        T = 298.0, [unit = u"K", description = "Temperature"]
        p₀ = 2.778e-7, [unit = u"m/s", description = "Precipitation rate (1 mm/h = 2.778e-7 m/s)"]
    end

    vars = @variables begin
        E(t), [unit = u"1", description = "Collection efficiency (dimensionless)"]
        Λ_particle(t), [unit = u"s^-1", description = "Particle scavenging coefficient (Eq. 20.57)"]
    end

    eqs = [
        E ~ slinn_collection_efficiency(D_p, U_t, d_p, ρ_p, T),  # Eq. 20.53
        Λ_particle ~ particle_scavenging_coeff(E, p₀, D_p),       # Eq. 20.57
    ]

    return System(eqs, t, vars, params; name)
end

# ---------------------------------------------------------------------------
# Eq. 20.7–20.9 — Wet Deposition Flux and Velocity
# F_w = C_rain * p₀     (Eq. 20.7)
# w_r = C_rain / C_air   (Eq. 20.6)
# u_w = F_w / C_air      (Eq. 20.8, = w_r * p₀, Eq. 20.9)
# ---------------------------------------------------------------------------
"""
    WetDepositionFlux(; name = :WetDepositionFlux)

ModelingToolkit component for the net wet deposition flux and velocity
(Eqs. 20.7–20.9, Seinfeld & Pandis, 2006).
"""
@component function WetDepositionFlux(; name = :WetDepositionFlux)
    params = @parameters begin
        C_rain = 1.0e-3, [unit = u"mol/m^3", description = "Concentration in rainwater at ground"]
        C_air = 1.0e-6, [unit = u"mol/m^3", description = "Gas-phase concentration at ground"]
        p₀ = 2.778e-7, [unit = u"m/s", description = "Precipitation rate"]
    end

    vars = @variables begin
        F_w(t), [unit = u"mol/m^2/s", description = "Wet deposition flux (Eq. 20.7)"]
        w_r(t), [unit = u"1", description = "Washout ratio (Eq. 20.6) (dimensionless)"]
        u_w(t), [unit = u"m/s", description = "Wet deposition velocity (Eq. 20.8)"]
    end

    eqs = [
        F_w ~ C_rain * p₀,          # Eq. 20.7
        w_r ~ C_rain / C_air,        # Eq. 20.6
        u_w ~ F_w / C_air,           # Eq. 20.8 (= w_r * p₀, Eq. 20.9)
    ]

    return System(eqs, t, vars, params; name)
end
