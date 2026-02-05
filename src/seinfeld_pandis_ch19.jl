"""
    DryDeposition

ModelingToolkit.jl implementation of dry deposition equations from:
Seinfeld, J.H. and Pandis, S.N. (2006). Atmospheric Chemistry and Physics,
2nd Edition, Chapter 19: Dry Deposition.

This module implements:
- Gas deposition velocity using the resistance model (Eq. 19.2)
- Particle deposition velocity (Eq. 19.7)
- Aerodynamic resistance (Eq. 19.14)
- Quasi-laminar resistance for gases (Eq. 19.17) and particles (Eq. 19.27)
- Surface resistance using Wesely (1989) model (Eq. 19.50)
"""
module DryDeposition

using ModelingToolkit
using ModelingToolkit: t, D
using DynamicQuantities

export DryDepositionGas, DryDepositionParticle
export AerodynamicResistance, QuasiLaminarResistanceGas, QuasiLaminarResistanceParticle
export SurfaceResistance, ParticleSettling
export GasProperties, LandUseParameters

#=============================================================================
# Physical Constants
=============================================================================#

"""Physical constants used throughout the dry deposition calculations."""
const KAPPA = 0.4        # von Karman constant (dimensionless)
const G_ACCEL = 9.81     # Gravitational acceleration (m/s^2)
const K_BOLTZ = 1.38e-23 # Boltzmann constant (J/K)

#=============================================================================
# Table 19.2: Land-Use Parameters for Zhang et al. (2001) Particle Model
=============================================================================#

"""
Land-use dependent parameters from Table 19.2 (Zhang et al. 2001).

Fields:
- A: characteristic radius of collectors (m)
- alpha: parameter for impaction efficiency
- gamma: parameter for Brownian diffusion efficiency
"""
struct LandUseParameters
    A::Float64      # Characteristic radius (m)
    alpha::Float64  # Impaction parameter
    gamma::Float64  # Brownian diffusion parameter
end

# Table 19.2 values for selected land use categories
const LANDUSE_GRASS = LandUseParameters(2.0e-3, 1.2, 0.54)
const LANDUSE_DECIDUOUS_FOREST = LandUseParameters(5.0e-3, 0.6, 0.56)
const LANDUSE_DESERT = LandUseParameters(10.0e-3, 50.0, 0.54)

#=============================================================================
# Table 19.4: Gas Properties
=============================================================================#

"""
Gas-specific properties from Table 19.4 for surface resistance calculations.

Fields:
- D_ratio: D_H2O/D_gas ratio for diffusivity
- H_star: effective Henry's law constant (M/atm)
- f_0: reactivity factor (0-1)
"""
struct GasProperties
    D_ratio::Float64   # D_H2O/D_gas
    H_star::Float64    # Effective Henry's law constant (M/atm)
    f_0::Float64       # Reactivity factor
end

# Table 19.4 values for selected gases
const GAS_SO2 = GasProperties(1.9, 1.0e5, 0.0)
const GAS_O3 = GasProperties(1.6, 0.01, 1.0)
const GAS_NO2 = GasProperties(1.6, 0.01, 0.1)
const GAS_HNO3 = GasProperties(1.9, 1.0e14, 0.0)

#=============================================================================
# Aerodynamic Resistance (Eq. 19.14)
=============================================================================#

"""
    AerodynamicResistance(; name=:AerodynamicResistance)

Calculates aerodynamic resistance for neutral stability conditions.

Implements Equation 19.14:
    r_a = (1/(kappa * u_star)) * ln(z/z_0)

where:
- kappa = 0.4 (von Karman constant)
- u_star: friction velocity (m/s)
- z: reference height (m)
- z_0: roughness length (m)
"""
@component function AerodynamicResistance(; name=:AerodynamicResistance)
    @constants begin
        kappa = KAPPA, [description = "von Karman constant (dimensionless)", unit = u"1"]
    end

    @parameters begin
        z, [description = "Reference height above surface", unit = u"m"]
        z_0, [description = "Surface roughness length", unit = u"m"]
    end

    @variables begin
        u_star(t), [description = "Friction velocity", unit = u"m/s"]
        r_a(t), [description = "Aerodynamic resistance", unit = u"s/m"]
    end

    eqs = [
        # Eq. 19.14 - Aerodynamic resistance for neutral stability
        r_a ~ (1 / (kappa * u_star)) * log(z / z_0),
    ]

    return System(eqs, t; name)
end

#=============================================================================
# Quasi-Laminar Resistance for Gases (Eq. 19.17)
=============================================================================#

"""
    QuasiLaminarResistanceGas(; name=:QuasiLaminarResistanceGas)

Calculates quasi-laminar (boundary layer) resistance for gases.

Implements Equation 19.17:
    r_b = 5 * Sc^(2/3) / u_star

where Sc = nu/D_g is the Schmidt number.
"""
@component function QuasiLaminarResistanceGas(; name=:QuasiLaminarResistanceGas)
    @parameters begin
        nu, [description = "Kinematic viscosity of air", unit = u"m^2/s"]
        D_g, [description = "Molecular diffusivity of gas in air", unit = u"m^2/s"]
    end

    @variables begin
        u_star(t), [description = "Friction velocity", unit = u"m/s"]
        Sc(t), [description = "Schmidt number (dimensionless)", unit = u"1"]
        r_b(t), [description = "Quasi-laminar resistance for gas", unit = u"s/m"]
    end

    eqs = [
        # Schmidt number definition
        Sc ~ nu / D_g,
        # Eq. 19.17 - Quasi-laminar resistance for gases
        r_b ~ 5 * Sc^(2/3) / u_star,
    ]

    return System(eqs, t; name)
end

#=============================================================================
# Particle Settling Velocity (Eq. 19.18) and Brownian Diffusivity (Eq. 19.20)
=============================================================================#

"""
    ParticleSettling(; name=:ParticleSettling)

Calculates particle settling velocity and Brownian diffusivity.

Implements:
- Equation 19.18 (Settling velocity):
    v_s = rho_p * D_p^2 * g * C_c / (18 * mu)

- Equation 19.20 (Brownian diffusivity):
    D = k * T * C_c / (3 * pi * mu * D_p)
"""
@component function ParticleSettling(; name=:ParticleSettling)
    @constants begin
        g = G_ACCEL, [description = "Gravitational acceleration", unit = u"m/s^2"]
        k_B = K_BOLTZ, [description = "Boltzmann constant", unit = u"J/K"]
        pi_val = Float64(pi), [description = "Mathematical constant pi (dimensionless)", unit = u"1"]
    end

    @parameters begin
        rho_p, [description = "Particle density", unit = u"kg/m^3"]
        D_p, [description = "Particle diameter", unit = u"m"]
        mu, [description = "Dynamic viscosity of air", unit = u"Pa*s"]
        T, [description = "Temperature", unit = u"K"]
        C_c, [description = "Cunningham slip correction factor (dimensionless)", unit = u"1"]
    end

    @variables begin
        v_s(t), [description = "Particle settling velocity", unit = u"m/s"]
        D_diff(t), [description = "Brownian diffusivity of particle", unit = u"m^2/s"]
    end

    eqs = [
        # Eq. 19.18 - Particle settling velocity (Stokes settling)
        v_s ~ rho_p * D_p^2 * g * C_c / (18 * mu),
        # Eq. 19.20 - Brownian diffusivity
        D_diff ~ k_B * T * C_c / (3 * pi_val * mu * D_p),
    ]

    return System(eqs, t; name)
end

#=============================================================================
# Quasi-Laminar Resistance for Particles (Eq. 19.27 - Zhang et al. 2001)
=============================================================================#

"""
    QuasiLaminarResistanceParticle(; name=:QuasiLaminarResistanceParticle)

Calculates quasi-laminar resistance for particles using Zhang et al. (2001) model.

Implements Equation 19.27:
    r_b = 1 / (3 * u_star * (E_B + E_IM + E_IN) * R_1)

where:
- E_B = Sc^(-gamma) (Brownian diffusion, Eq. 19.21)
- E_IM = (St/(alpha + St))^2 (impaction, Eq. 19.22 with beta=2)
- E_IN = 0.5 * (D_p/A)^2 (interception, Eq. 19.25)
- R_1 = exp(-sqrt(St)) (sticking fraction, Eq. 19.26)
- St = v_s * u_star / (g * A) (Stokes number, Eq. 19.24)
"""
@component function QuasiLaminarResistanceParticle(; name=:QuasiLaminarResistanceParticle)
    @constants begin
        g = G_ACCEL, [description = "Gravitational acceleration", unit = u"m/s^2"]
    end

    @parameters begin
        nu, [description = "Kinematic viscosity of air", unit = u"m^2/s"]
        A, [description = "Characteristic radius of collectors", unit = u"m"]
        alpha, [description = "Impaction efficiency parameter (dimensionless)", unit = u"1"]
        gamma, [description = "Brownian diffusion efficiency parameter (dimensionless)", unit = u"1"]
        D_p, [description = "Particle diameter", unit = u"m"]
    end

    @variables begin
        u_star(t), [description = "Friction velocity", unit = u"m/s"]
        v_s(t), [description = "Particle settling velocity", unit = u"m/s"]
        D_diff(t), [description = "Brownian diffusivity of particle", unit = u"m^2/s"]
        Sc(t), [description = "Schmidt number (dimensionless)", unit = u"1"]
        St(t), [description = "Stokes number (dimensionless)", unit = u"1"]
        E_B(t), [description = "Brownian diffusion collection efficiency (dimensionless)", unit = u"1"]
        E_IM(t), [description = "Impaction collection efficiency (dimensionless)", unit = u"1"]
        E_IN(t), [description = "Interception collection efficiency (dimensionless)", unit = u"1"]
        R_1(t), [description = "Sticking fraction (dimensionless)", unit = u"1"]
        r_b(t), [description = "Quasi-laminar resistance for particle", unit = u"s/m"]
    end

    eqs = [
        # Schmidt number for particles
        Sc ~ nu / D_diff,
        # Eq. 19.24 - Stokes number for vegetated surfaces
        St ~ v_s * u_star / (g * A),
        # Eq. 19.21 - Brownian diffusion collection efficiency
        E_B ~ Sc^(-gamma),
        # Eq. 19.22 - Impaction collection efficiency (beta = 2)
        E_IM ~ (St / (alpha + St))^2,
        # Eq. 19.25 - Interception collection efficiency
        E_IN ~ 0.5 * (D_p / A)^2,
        # Eq. 19.26 - Sticking fraction (rebound correction)
        R_1 ~ exp(-sqrt(St)),
        # Eq. 19.27 - Quasi-laminar resistance for particles (Zhang et al. 2001)
        r_b ~ 1 / (3 * u_star * (E_B + E_IM + E_IN) * R_1),
    ]

    return System(eqs, t; name)
end

#=============================================================================
# Surface Resistance - Wesely (1989) Model (Eq. 19.50)
=============================================================================#

"""
    SurfaceResistance(; name=:SurfaceResistance)

Calculates surface (canopy) resistance for gases using Wesely (1989) model.

Implements Equation 19.50:
    r_c = (1/(r_st + r_m) + 1/r_lu + 1/(r_dc + r_cl) + 1/(r_ac + r_gs))^(-1)

and component resistances:
- Eq. 19.51: Stomatal resistance
    r_st = r_i * [1 + (200/(G+0.1))^2] * [400/(T_s*(40-T_s))]
- Eq. 19.52: Mesophyll resistance
    r_m = 1 / (3.3e-4 * H_star + 100 * f_0)

Note: The stomatal resistance formula uses empirical coefficients with implicit units.
Temperature T_s is expected in Celsius (valid range: 0-40C for the formula).
Solar irradiance G is expected in W/m^2.
"""
@component function SurfaceResistance(; name=:SurfaceResistance)
    @constants begin
        # Conversion constant for temperature
        T_ref = 273.15, [description = "Celsius to Kelvin offset", unit = u"K"]
        # Empirical constants for stomatal resistance (Eq. 19.51)
        G_ref = 200.0, [description = "Reference irradiance for stomatal response", unit = u"W/m^2"]
        G_offset = 0.1, [description = "Small offset to prevent division by zero", unit = u"W/m^2"]
        T_coeff = 400.0, [description = "Temperature coefficient for stomatal response", unit = u"K^2"]
        T_max = 40.0, [description = "Maximum temperature for stomatal function", unit = u"K"]
        # Empirical constants for mesophyll resistance (Eq. 19.52)
        # The coefficient 3.3e-4 has units of m/s per (M/atm) = m*atm/(s*M)
        H_coeff = 3.3e-4, [description = "Henry coefficient for mesophyll", unit = u"m/s"]
        f_coeff = 100.0, [description = "Reactivity coefficient for mesophyll", unit = u"m/s"]
    end

    @parameters begin
        # Gas properties (from Table 19.4)
        # H_star is dimensionless in this formulation (effective Henry constant normalized)
        H_star, [description = "Effective Henry's law constant (dimensionless, normalized)", unit = u"1"]
        f_0, [description = "Reactivity factor (dimensionless)", unit = u"1"]
        D_ratio, [description = "D_H2O/D_gas diffusivity ratio (dimensionless)", unit = u"1"]
        # Land-use dependent resistances (from Wesely 1989 tables)
        r_i, [description = "Minimum stomatal resistance", unit = u"s/m"]
        r_lu, [description = "Cuticular/outer surface resistance", unit = u"s/m"]
        r_dc, [description = "Lower canopy buoyant convection resistance", unit = u"s/m"]
        r_cl, [description = "Lower canopy exposed surfaces resistance", unit = u"s/m"]
        r_ac, [description = "In-canopy aerodynamic resistance", unit = u"s/m"]
        r_gs, [description = "Ground surface resistance", unit = u"s/m"]
    end

    @variables begin
        G(t), [description = "Solar irradiance at surface", unit = u"W/m^2"]
        T_s(t), [description = "Surface temperature", unit = u"K"]
        T_s_C(t), [description = "Surface temperature offset from reference", unit = u"K"]
        G_factor(t), [description = "Irradiance factor for stomatal resistance (dimensionless)", unit = u"1"]
        T_factor(t), [description = "Temperature factor for stomatal resistance (dimensionless)", unit = u"1"]
        r_st(t), [description = "Stomatal resistance", unit = u"s/m"]
        r_m(t), [description = "Mesophyll resistance", unit = u"s/m"]
        r_c(t), [description = "Surface (canopy) resistance", unit = u"s/m"]
    end

    eqs = [
        # Temperature offset from 273.15 K (gives value in "Celsius-like" units but with K unit)
        T_s_C ~ T_s - T_ref,
        # Eq. 19.51 - Stomatal resistance factors (dimensionless)
        # G_factor = 1 + (200/(G+0.1))^2
        G_factor ~ 1 + (G_ref / (G + G_offset))^2,
        # T_factor = 400/(T_s_C * (40 - T_s_C)) where T_s_C is in Celsius-equivalent
        # T_coeff has units K^2 so that T_factor becomes dimensionless
        T_factor ~ T_coeff / (T_s_C * (T_max - T_s_C)),
        # Eq. 19.51 - Stomatal resistance
        r_st ~ D_ratio * r_i * G_factor * T_factor,
        # Eq. 19.52 - Mesophyll resistance
        # r_m = 1 / (3.3e-4 * H_star + 100 * f_0) with coefficients having s/m units
        r_m ~ 1 / (H_coeff * H_star + f_coeff * f_0),
        # Eq. 19.50 - Total surface resistance (parallel pathways)
        r_c ~ 1 / (1 / (r_st + r_m) + 1 / r_lu + 1 / (r_dc + r_cl) + 1 / (r_ac + r_gs)),
    ]

    return System(eqs, t; name)
end

#=============================================================================
# Complete Gas Deposition Velocity System (Eq. 19.2)
=============================================================================#

"""
    DryDepositionGas(; name=:DryDepositionGas)

Complete gas dry deposition velocity model using resistance approach.

Implements Equation 19.2:
    v_d = 1 / (r_a + r_b + r_c)

and Equation 19.1 for flux:
    F = -v_d * C

Composes:
- AerodynamicResistance (Eq. 19.14)
- QuasiLaminarResistanceGas (Eq. 19.17)
- SurfaceResistance (Eq. 19.50-19.52)

Inputs (must be specified externally):
- u_star: Friction velocity (m/s)
- C: Gas concentration at reference height (mol/m^3)
- G: Solar irradiance at surface (W/m^2) - via surf subsystem
- T_s: Surface temperature (K) - via surf subsystem
"""
@component function DryDepositionGas(; name=:DryDepositionGas)
    # Create subsystems
    @named aero = AerodynamicResistance()
    @named qlam = QuasiLaminarResistanceGas()
    @named surf = SurfaceResistance()

    @variables begin
        u_star(t), [description = "Friction velocity", unit = u"m/s", input = true]
        C(t), [description = "Gas concentration at reference height", unit = u"mol/m^3", input = true]
        G(t), [description = "Solar irradiance at surface", unit = u"W/m^2", input = true]
        T_s(t), [description = "Surface temperature", unit = u"K", input = true]
        r_t(t), [description = "Total resistance", unit = u"s/m"]
        v_d(t), [description = "Deposition velocity", unit = u"m/s"]
        F(t), [description = "Deposition flux (downward positive)", unit = u"mol/(m^2*s)"]
    end

    eqs = [
        # Connect friction velocity to subsystems
        aero.u_star ~ u_star,
        qlam.u_star ~ u_star,
        # Connect surface conditions to surface resistance subsystem
        surf.G ~ G,
        surf.T_s ~ T_s,
        # Eq. 19.2 - Total resistance (resistance in series)
        r_t ~ aero.r_a + qlam.r_b + surf.r_c,
        # Deposition velocity
        v_d ~ 1 / r_t,
        # Eq. 19.1 - Deposition flux (negative sign indicates downward flux)
        F ~ -v_d * C,
    ]

    return System(eqs, t; systems = [aero, qlam, surf], name)
end

#=============================================================================
# Complete Particle Deposition Velocity System (Eq. 19.7)
=============================================================================#

"""
    DryDepositionParticle(; name=:DryDepositionParticle)

Complete particle dry deposition velocity model.

Implements Equation 19.7:
    v_d = 1/(r_a + r_b + r_a*r_b*v_s) + v_s

and Equation 19.1 for flux:
    F = -v_d * C

Composes:
- AerodynamicResistance (Eq. 19.14)
- ParticleSettling (Eq. 19.18, 19.20)
- QuasiLaminarResistanceParticle (Eq. 19.27)

Inputs (must be specified externally):
- u_star: Friction velocity (m/s)
- C: Particle concentration at reference height (kg/m^3)
"""
@component function DryDepositionParticle(; name=:DryDepositionParticle)
    # Create subsystems
    @named aero = AerodynamicResistance()
    @named settling = ParticleSettling()
    @named qlam = QuasiLaminarResistanceParticle()

    @variables begin
        u_star(t), [description = "Friction velocity", unit = u"m/s", input = true]
        C(t), [description = "Particle concentration at reference height", unit = u"kg/m^3", input = true]
        v_d(t), [description = "Deposition velocity", unit = u"m/s"]
        F(t), [description = "Deposition flux (downward positive)", unit = u"kg/(m^2*s)"]
    end

    eqs = [
        # Connect friction velocity to subsystems
        aero.u_star ~ u_star,
        qlam.u_star ~ u_star,
        # Connect settling velocity and diffusivity to quasi-laminar subsystem
        qlam.v_s ~ settling.v_s,
        qlam.D_diff ~ settling.D_diff,
        # Eq. 19.7 - Particle deposition velocity
        # Accounts for gravitational settling in addition to turbulent deposition
        v_d ~ 1 / (aero.r_a + qlam.r_b + aero.r_a * qlam.r_b * settling.v_s) + settling.v_s,
        # Eq. 19.1 - Deposition flux (negative sign indicates downward flux)
        F ~ -v_d * C,
    ]

    return System(eqs, t; systems = [aero, settling, qlam], name)
end

end # module DryDeposition
