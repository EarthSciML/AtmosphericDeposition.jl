export AirRefreshingLimitation, CloudIceUptakeLimitation, WetScavengingLimitations

# ============================================================================
# Air Refreshing and Cloud Ice Uptake Limitations on Wet Scavenging
# Reference: Luo, G., & Yu, F. (2023). Impact of air refreshing and cloud ice
#   uptake limitations on vertical profiles and wet depositions of nitrate,
#   ammonium, and sulfate. Geophysical Research Letters, 50, e2023GL104258.
#   https://doi.org/10.1029/2023GL104258
# ============================================================================

# ---------------------------------------------------------------------------
# Constants for unit handling
# ---------------------------------------------------------------------------
@constants begin
    T_upper_luo = 220.0,
    [unit = u"K",
        description = "Upper temperature threshold for γ (Eq. 15)"]
    T_lower_luo = 209.0,
    [unit = u"K",
        description = "Lower temperature threshold for γ (Eq. 15)"]
    γ_base = 3e-3,
    [unit = u"1", description = "Base uptake efficiency (Eq. 15) (dimensionless)"]
    γ_delta = 4e-3,
    [unit = u"1",
        description = "Temperature-dependent uptake coefficient (Eq. 15) (dimensionless)"]
    zero_dimless = 0.0, [unit = u"1", description = "Zero (dimensionless)"]
    one_dimless = 1.0, [unit = u"1", description = "One (dimensionless)"]
    kinetic_prefactor_si = 2.749064e-2,
    [unit = u"s/m",
        description = "Kinetic theory prefactor for ice uptake, converted from CGS to SI (Eq. 14)"]
    M_ref = 1.0,
    [unit = u"g/mol",
        description = "Reference molar mass for non-dimensionalization of square root"]
    T_ref = 1.0,
    [unit = u"K",
        description = "Reference temperature for non-dimensionalization of square root"]
end

# ---------------------------------------------------------------------------
# Eq. 10 — Isotropic turbulence velocity from TKE
# u' = v' = w' = √(2/3 · TKE)
# Based on Stull (1988) and Pinto (1998).
# ---------------------------------------------------------------------------
"""
    turbulence_velocity(TKE)

Compute isotropic turbulence velocity from turbulence kinetic energy
(Eq. 10, Luo & Yu, 2023). Assumes u' = v' = w' per Pinto (1998).
"""
function turbulence_velocity(TKE)
    return sqrt(2 // 3 * TKE)
end

# ---------------------------------------------------------------------------
# Eq. 11 — Cloudy/rainy air refreshing rate Kᵢ
# Kᵢ = √(2/3·TKE) · [1/(f^(1/2)·Δx) + 1/(f^(1/2)·Δy) + 1/Δz]
# Derived from Eqs. 6–10 by summing mass exchange across all faces
# of the cloud volume and dividing by the cloud volume.
# ---------------------------------------------------------------------------
"""
    cloudy_air_refreshing_rate(f, TKE, Δx, Δy, Δz)

Compute the cloudy (rainy) air refreshing rate Kᵢ (Eq. 11, Luo & Yu, 2023).
Cloud fractions in x and y are f^(1/2); z-direction fraction is 1.
"""
function cloudy_air_refreshing_rate(f, TKE, Δx, Δy, Δz)
    v_turb = turbulence_velocity(TKE) # Eq. 10
    # Eq. 11: sum of face-area/volume contributions times turbulence velocity
    Kᵢ = v_turb * (1 / (sqrt(f) * Δx) + 1 / (sqrt(f) * Δy) + 1 / Δz)
    return Kᵢ
end

# ---------------------------------------------------------------------------
# Eq. 5 — Grid refreshing time τ_A
# τ_A = (1 - f) / (f · Kᵢ)
# ---------------------------------------------------------------------------
"""
    grid_refreshing_time(f, Kᵢ)

Compute the grid refreshing time τ_A (Eq. 5, Luo & Yu, 2023).
"""
function grid_refreshing_time(f, Kᵢ)
    return (1 - f) / (f * Kᵢ)
end

# ---------------------------------------------------------------------------
# Eq. 2 — Air refreshing limited removal rate R_A
# R_A = 1 / (1/(f·Rᵢ) + τ_A)
# When τ_A → 0 (strong mixing), R_A → f·Rᵢ (standard model).
# When τ_A → ∞ (no mixing), R_A → 0.
# ---------------------------------------------------------------------------
"""
    air_refreshing_limited_rate(f, Rᵢ, τ_A)

Compute the air refreshing limited grid mean mass loss rate R_A
(Eq. 2, Luo & Yu, 2023).
"""
function air_refreshing_limited_rate(f, Rᵢ, τ_A)
    return 1 / (1 / (f * Rᵢ) + τ_A)
end

# ---------------------------------------------------------------------------
# Eq. 15 — Uptake efficiency γ for HNO₃ (Hudson et al., 2002)
# γ = 3×10⁻³ + 4×10⁻³ · max(0, min(1, (220 - T)/(220 - 209)))
# ---------------------------------------------------------------------------
"""
    hno3_uptake_efficiency(T)

Compute the uptake efficiency γ of nitric acid on ice crystals
(Eq. 15, Luo & Yu, 2023; based on Hudson et al., 2002).
Returns γ = 0.003 for T ≥ 220 K and γ = 0.007 for T ≤ 209 K.
"""
function hno3_uptake_efficiency(T)
    # Use constants with units to avoid unit mismatch
    return γ_base +
           γ_delta * max(zero_dimless,
        min(one_dimless, (T_upper_luo - T) / (T_upper_luo - T_lower_luo)))
end

# ---------------------------------------------------------------------------
# Eq. 14 — Cloud ice uptake rate for HNO₃ (Jacob, 2000)
# R_{U,HNO3} = N_I · S_I / (r/D_g + 2.749064e-4 · √M / (γ · √T))
#
# Note: The original equation uses CGS units (cm⁻³, cm², cm, cm²/s).
# We implement in SI with explicit unit conversion.
#
# In CGS: R_U = N_I[cm⁻³] · S_I[cm²] / (r[cm]/D_g[cm²/s] + 2.749064e-4·√(M[g/mol])/(γ·√(T[K])))
# The denominator has units s/cm.
# In SI: r[m]/D_g[m²/s] has units s/m. So we convert the prefactor from
# s/cm to s/m by multiplying by 100 → 2.749064e-2 s/m.
# Then N_I[m⁻³]·S_I[m²] / (s/m) → m⁻¹/(s/m) = s⁻¹ ✓
#
# To handle √M and √T with units, we non-dimensionalize by reference values.
# ---------------------------------------------------------------------------
"""
    cloud_ice_uptake_rate(N_I, S_I, r, D_g, M, T, γ)

Compute the cloud ice uptake rate R_U for HNO₃ (Eq. 14, Luo & Yu, 2023;
based on Jacob, 2000). All inputs in SI units.
"""
function cloud_ice_uptake_rate(N_I, S_I, r, D_g, M, T, γ)
    # Non-dimensionalize M and T before taking square roots to avoid
    # fractional unit powers, then multiply by the SI prefactor (s/m).
    # M/M_ref and T/T_ref are dimensionless; kinetic_prefactor_si is s/m.
    kinetic_term = kinetic_prefactor_si * sqrt(M / M_ref) / (γ * sqrt(T / T_ref))
    R_U = N_I * S_I / (r / D_g + kinetic_term)
    return R_U
end

# ---------------------------------------------------------------------------
# Eq. 13 — Air refreshing limited cloud ice uptake rate R_{A,U}
# R_{A,U} = f · R_U · Kᵢ / (Kᵢ + (1 - f) · R_U)
# This is an application of Eq. 2 structure to R_U.
# ---------------------------------------------------------------------------
"""
    air_refreshing_limited_ice_uptake_rate(f, R_U, Kᵢ)

Compute the air refreshing limited cloud ice uptake rate R_{A,U}
(Eq. 13, Luo & Yu, 2023).
"""
function air_refreshing_limited_ice_uptake_rate(f, R_U, Kᵢ)
    return f * R_U * Kᵢ / (Kᵢ + (1 - f) * R_U)
end

# ---------------------------------------------------------------------------
# Eq. 12 — Cold cloud rainout efficiency F_I
# F_I = 1 - exp(-R_{A,U} · Δt)
# ---------------------------------------------------------------------------
"""
    cold_cloud_rainout_efficiency(R_AU, Δt)

Compute the cold cloud rainout efficiency F_I (Eq. 12, Luo & Yu, 2023).
F_I equals the fraction of water-soluble aerosols in cloud ice due to uptake.
"""
function cold_cloud_rainout_efficiency(R_AU, Δt)
    return 1 - exp(-R_AU * Δt)
end

# ============================================================================
# ModelingToolkit Components
# ============================================================================

"""
    AirRefreshingLimitation(; name = :AirRefreshingLimitation)

ModelingToolkit component implementing the air refreshing limitation on
wet scavenging from Luo & Yu (2023), Eqs. 2, 5, 10, 11.

This parameterization accounts for the fact that species in subgrid
cloud-free (rain-free) air need time to be mixed with those in cloudy
(rainy) air before being influenced by wet scavenging. The standard
well-mixed assumption can overestimate wet scavenging when the subgrid
air mixing time scale is comparable to the model time step.

The air refreshing limited removal rate `R_A` replaces the standard
`f·Rᵢ` in wet scavenging calculations. When turbulent mixing is strong
(large TKE), `R_A ≈ f·Rᵢ`; when mixing is weak, `R_A < f·Rᵢ`.
"""
@component function AirRefreshingLimitation(; name = :AirRefreshingLimitation)
    @parameters begin
        f, [unit = u"1", description = "Cloud or rainfall fraction (dimensionless)"]
        Rᵢ, [unit = u"s^-1", description = "In-cloud or under-rain removing rate"]
        TKE, [unit = u"m^2/s^2", description = "Turbulence kinetic energy"]
        Δx, [unit = u"m", description = "Grid spacing in x direction"]
        Δy, [unit = u"m", description = "Grid spacing in y direction"]
        Δz, [unit = u"m", description = "Grid spacing in z direction"]
    end

    @variables begin
        Kᵢ(t), [unit = u"s^-1",
            description = "Cloudy/rainy air refreshing rate (Eq. 11)"]
        τ_A(t), [unit = u"s",
            description = "Grid refreshing time (Eq. 5)"]
        R_A(t),
        [unit = u"s^-1",
            description = "Air refreshing limited removal rate (Eq. 2)"]
    end

    eqs = [
        Kᵢ ~ cloudy_air_refreshing_rate(f, TKE, Δx, Δy, Δz), # Eq. 11
        τ_A ~ grid_refreshing_time(f, Kᵢ),                      # Eq. 5
        R_A ~ air_refreshing_limited_rate(f, Rᵢ, τ_A)          # Eq. 2
    ]

    return System(eqs, t; name)
end

"""
    CloudIceUptakeLimitation(; name = :CloudIceUptakeLimitation)

ModelingToolkit component implementing the cloud ice uptake limitation
on cold cloud wet scavenging from Luo & Yu (2023), Eqs. 12–15.

In cold clouds, water-soluble aerosols are captured by ice crystals
via coagulation and then removed by precipitation. The cold cloud rainout
efficiency `F_I` is limited by the cloud ice uptake rate, which is
approximated using the HNO₃ uptake rate from Jacob (2000).

The uptake efficiency `γ` depends on temperature following
Hudson et al. (2002): γ = 0.003 at T ≥ 220 K, increasing linearly
to 0.007 at T ≤ 209 K.
"""
@component function CloudIceUptakeLimitation(; name = :CloudIceUptakeLimitation)
    @parameters begin
        f, [unit = u"1", description = "Cloud fraction (dimensionless)"]
        TKE, [unit = u"m^2/s^2", description = "Turbulence kinetic energy"]
        Δx, [unit = u"m", description = "Grid spacing in x direction"]
        Δy, [unit = u"m", description = "Grid spacing in y direction"]
        Δz, [unit = u"m", description = "Grid spacing in z direction"]
        Δt, [unit = u"s", description = "Time step"]
        N_I, [unit = u"m^-3", description = "Ice crystal number concentration"]
        S_I, [unit = u"m^2", description = "Ice crystal surface area"]
        r_ice, [unit = u"m", description = "Ice crystal radius"]
        D_g, [unit = u"m^2/s", description = "Gas-phase diffusion coefficient"]
        M, [unit = u"g/mol", description = "Molar mass"]
        T, [unit = u"K", description = "Temperature"]
    end

    @variables begin
        γ(t),
        [unit = u"1",
            description = "Uptake efficiency of HNO₃ on ice (Eq. 15) (dimensionless)"]
        R_U(t), [unit = u"s^-1",
            description = "Cloud ice uptake rate (Eq. 14)"]
        Kᵢ(t), [unit = u"s^-1",
            description = "Cloudy air refreshing rate (Eq. 11)"]
        R_AU(t),
        [unit = u"s^-1",
            description = "Air refreshing limited cloud ice uptake rate (Eq. 13)"]
        F_I(t),
        [unit = u"1",
            description = "Cold cloud rainout efficiency (Eq. 12) (dimensionless)"]
    end

    eqs = [
        γ ~ hno3_uptake_efficiency(T),                                    # Eq. 15
        R_U ~ cloud_ice_uptake_rate(N_I, S_I, r_ice, D_g, M, T, γ),     # Eq. 14
        Kᵢ ~ cloudy_air_refreshing_rate(f, TKE, Δx, Δy, Δz),           # Eq. 11
        R_AU ~ air_refreshing_limited_ice_uptake_rate(f, R_U, Kᵢ),      # Eq. 13
        F_I ~ cold_cloud_rainout_efficiency(R_AU, Δt)                    # Eq. 12
    ]

    return System(eqs, t; name)
end

"""
    WetScavengingLimitations(; name = :WetScavengingLimitations)

ModelingToolkit component combining both the air refreshing limitation and
cloud ice uptake limitation on wet scavenging from Luo & Yu (2023).

This provides a complete parameterization of the two novel approaches:

 1. Air refreshing limitation (Section 2.1): Reduces wet scavenging rate
    when subgrid air mixing is slow relative to the removal rate.
 2. Cloud ice uptake limitation (Section 2.2): Limits cold cloud rainout
    efficiency by the rate at which aerosols are captured by ice crystals.
"""
@component function WetScavengingLimitations(; name = :WetScavengingLimitations)
    @parameters begin
        f, [unit = u"1", description = "Cloud or rainfall fraction (dimensionless)"]
        Rᵢ, [unit = u"s^-1", description = "In-cloud or under-rain removing rate"]
        TKE, [unit = u"m^2/s^2", description = "Turbulence kinetic energy"]
        Δx, [unit = u"m", description = "Grid spacing in x direction"]
        Δy, [unit = u"m", description = "Grid spacing in y direction"]
        Δz, [unit = u"m", description = "Grid spacing in z direction"]
        Δt, [unit = u"s", description = "Time step"]
        N_I, [unit = u"m^-3", description = "Ice crystal number concentration"]
        S_I, [unit = u"m^2", description = "Ice crystal surface area"]
        r_ice, [unit = u"m", description = "Ice crystal radius"]
        D_g, [unit = u"m^2/s", description = "Gas-phase diffusion coefficient"]
        M, [unit = u"g/mol", description = "Molar mass"]
        T, [unit = u"K", description = "Temperature"]
    end

    @variables begin
        Kᵢ(t), [unit = u"s^-1",
            description = "Cloudy/rainy air refreshing rate (Eq. 11)"]
        τ_A(t), [unit = u"s",
            description = "Grid refreshing time (Eq. 5)"]
        R_A(t),
        [unit = u"s^-1",
            description = "Air refreshing limited removal rate (Eq. 2)"]
        γ(t),
        [unit = u"1",
            description = "Uptake efficiency of HNO₃ on ice (Eq. 15) (dimensionless)"]
        R_U(t), [unit = u"s^-1",
            description = "Cloud ice uptake rate (Eq. 14)"]
        R_AU(t),
        [unit = u"s^-1",
            description = "Air refreshing limited cloud ice uptake rate (Eq. 13)"]
        F_I(t),
        [unit = u"1",
            description = "Cold cloud rainout efficiency (Eq. 12) (dimensionless)"]
    end

    eqs = [
        # Air refreshing limitation (Section 2.1)
        Kᵢ ~ cloudy_air_refreshing_rate(f, TKE, Δx, Δy, Δz),    # Eq. 11
        τ_A ~ grid_refreshing_time(f, Kᵢ),                         # Eq. 5
        R_A ~ air_refreshing_limited_rate(f, Rᵢ, τ_A),             # Eq. 2
        # Cloud ice uptake limitation (Section 2.2)
        γ ~ hno3_uptake_efficiency(T),                               # Eq. 15
        R_U ~ cloud_ice_uptake_rate(N_I, S_I, r_ice, D_g, M, T, γ), # Eq. 14
        R_AU ~ air_refreshing_limited_ice_uptake_rate(f, R_U, Kᵢ),  # Eq. 13
        F_I ~ cold_cloud_rainout_efficiency(R_AU, Δt)               # Eq. 12
    ]

    return System(eqs, t; name)
end
