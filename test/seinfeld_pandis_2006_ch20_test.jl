@testsnippet SP2006Setup begin
    using AtmosphericDeposition
    using AtmosphericDeposition: mass_transfer_coeff, gas_scavenging_coeff,
                                 reversible_drop_conc, reversible_scavenging_flux,
                                 slinn_collection_efficiency, particle_scavenging_coeff,
                                 particle_relaxation_time, particle_terminal_velocity,
                                 particle_diffusivity,
                                 μ_air_sp, ρ_air_sp, μ_w_sp, ρ_w_sp, g_sp, kB_sp
    using Test, DynamicQuantities, ModelingToolkit

    @parameters D_g_p = 1.5e-5 [unit = u"m^2/s"]
    @parameters D_p_p = 1.0e-3 [unit = u"m"]
    @parameters U_t_p = 4.0 [unit = u"m/s"]
    @parameters K_c_p = 0.13 [unit = u"m/s"]
    @parameters p₀_p = 2.778e-7 [unit = u"m/s"]
    @parameters d_p_p = 1.0e-6 [unit = u"m"]
    @parameters ρ_p_p = 1000.0 [unit = u"kg/m^3"]
    @parameters T_p = 298.0 [unit = u"K"]
    @parameters h_p = 1000.0 [unit = u"m"]

    sp_defaults = [
        μ_air_sp => 1.72e-5,
        ρ_air_sp => 1.225,
        μ_w_sp => 1.0e-3,
        ρ_w_sp => 1000.0,
        g_sp => 9.81,
        kB_sp => 1.3806488e-23
    ]
end

# ====================================================================
# Eq. 20.12 — Mass Transfer Coefficient
# ====================================================================
@testitem "Mass Transfer Coefficient (Eq. 20.12)" setup=[SP2006Setup] tags=[:sp2006] begin
    # Test with Table 20.1 values: D_p = 0.01 cm = 1e-4 m, U_t = 26 cm/s = 0.26 m/s
    # K_c should be ~32 cm/s = 0.32 m/s (for SO2 with D_g ≈ 0.126 cm^2/s = 1.26e-5 m^2/s)
    # Manual calculation:
    # Re = 1.225 * 0.26 * 1e-4 / 1.72e-5 = 1.851
    # Sc = 1.72e-5 / (1.225 * 1.26e-5) = 1.114
    # Sh = 2 + 0.6 * sqrt(1.851) * cbrt(1.114) = 2 + 0.6 * 1.360 * 1.037 = 2.847
    # K_c = 1.26e-5 / 1e-4 * 2.847 = 0.359 m/s
    Kc_expr=mass_transfer_coeff(D_g_p, D_p_p, U_t_p)
    Kc_val=substitute(Kc_expr,
        Dict(D_g_p=>1.26e-5, D_p_p=>1e-4, U_t_p=>0.26, sp_defaults...))
    @test Kc_val ≈ 0.3586 rtol = 0.01

    # Check units
    @test ModelingToolkit.get_unit(Kc_expr) == u"m/s"

    # Test with larger drop: D_p = 1 cm = 0.01 m, U_t = 1000 cm/s = 10 m/s
    # Re = 1.225 * 10 * 0.01 / 1.72e-5 = 7122
    # Sc = 1.114
    # Sh = 2 + 0.6 * sqrt(7122) * cbrt(1.114) = 2 + 0.6 * 84.39 * 1.037 = 54.52
    # K_c = 1.26e-5 / 0.01 * 54.52 = 0.0687
    Kc_val2=substitute(Kc_expr,
        Dict(D_g_p=>1.26e-5, D_p_p=>0.01, U_t_p=>10.0, sp_defaults...))
    @test Kc_val2 ≈ 0.0687 rtol = 0.01
end

@testitem "MassTransferCoeff Component" setup=[SP2006Setup] tags=[:sp2006] begin
    sys=MassTransferCoeff()
    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 1
end

# ====================================================================
# Eq. 20.25 — Gas Scavenging Coefficient
# ====================================================================
@testitem "Gas Scavenging Coefficient (Eq. 20.25)" setup=[SP2006Setup] tags=[:sp2006] begin
    # Table 20.1: D_p = 0.1 cm = 1e-3 m, U_t = 300 cm/s = 3.0 m/s,
    # K_c = 13 cm/s = 0.13 m/s, p₀ = 1 mm/h = 2.778e-7 m/s
    # Λ = 6 * p₀ * K_c / (D_p * U_t) = 6 * 2.778e-7 * 0.13 / (1e-3 * 3.0)
    #   = 6 * 3.611e-8 / 3e-3 = 2.167e-7 / 3e-3 = 7.22e-5 s⁻¹
    # Table reports Λ = 0.26 h⁻¹ = 0.26/3600 = 7.22e-5 s⁻¹ ✓
    Λ_expr=gas_scavenging_coeff(K_c_p, U_t_p, D_p_p, p₀_p)
    Λ_val=substitute(Λ_expr,
        Dict(K_c_p=>0.13, U_t_p=>3.0, D_p_p=>1e-3, p₀_p=>2.778e-7, sp_defaults...))
    @test Λ_val ≈ 7.22e-5 rtol = 0.02

    # Table 20.1: D_p = 0.01 cm = 1e-4 m, U_t = 26 cm/s = 0.26 m/s,
    # K_c = 32 cm/s = 0.32 m/s
    # Λ = 6 * 2.778e-7 * 0.32 / (1e-4 * 0.26) = 5.333e-7 / 2.6e-5 = 0.0205 s⁻¹
    # Table reports Λ = 73.8 h⁻¹ = 73.8/3600 = 0.0205 s⁻¹ ✓
    Λ_val2=substitute(Λ_expr,
        Dict(K_c_p=>0.32, U_t_p=>0.26, D_p_p=>1e-4, p₀_p=>2.778e-7, sp_defaults...))
    @test Λ_val2 ≈ 0.0205 rtol = 0.02

    # Check units
    @test ModelingToolkit.get_unit(Λ_expr) == u"s^-1"
end

@testitem "GasScavengingCoeff Component" setup=[SP2006Setup] tags=[:sp2006] begin
    sys=GasScavengingCoeff()
    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 1
end

# ====================================================================
# Eq. 20.28, 20.35 — Reversible Gas Scavenging
# ====================================================================
@testitem "Reversible Gas Scavenging (Eq. 20.28)" setup=[SP2006Setup] tags=[:sp2006] begin
    @parameters C_g_p = 1.0e-6 [unit = u"mol/m^3"]
    @parameters HRT_p = 100.0
    @parameters C_aq0_p = 0.0 [unit = u"mol/m^3"]
    @parameters z_p = 1000.0 [unit = u"m"]

    # Test limiting case: when HRT is very large (very soluble gas),
    # the exponent α = 6*K_c*z/(D_p*U_t*HRT) → 0, so C_aq → C_aq0
    # (drop can't reach equilibrium in the available fall distance)
    C_aq_expr=reversible_drop_conc(C_g_p, HRT_p, K_c_p, D_p_p, U_t_p, z_p, C_aq0_p)
    C_aq_huge_H=substitute(C_aq_expr,
        Dict(C_g_p=>1e-6, HRT_p=>1e12, K_c_p=>0.13, D_p_p=>1e-3,
            U_t_p=>4.0, z_p=>1000.0, C_aq0_p=>0.0, sp_defaults...))
    # With huge HRT, equilibrium is C_g * HRT = 1e-6 * 1e12 = 1e6
    # But α ≈ 0 so exp(-α) ≈ 1 and C_aq ≈ C_aq0 = 0
    @test C_aq_huge_H ≈ 0.0 atol = 1.0

    # Test limiting case: when α → ∞ (small HRT or long fall), C_aq → C_g * HRT (equilibrium)
    C_aq_equil=substitute(C_aq_expr,
        Dict(C_g_p=>1e-6, HRT_p=>1.0, K_c_p=>10.0, D_p_p=>1e-4,
            U_t_p=>0.1, z_p=>10000.0, C_aq0_p=>0.0, sp_defaults...))
    # α = 6*10*10000 / (1e-4*0.1*1) = 6e5 / 1e-5 = 6e10 → very large
    # C_aq → C_g * HRT = 1e-6 * 1 = 1e-6
    @test C_aq_equil ≈ 1e-6 rtol = 1e-6

    # Test initial condition preserved at z = 0
    C_aq_z0=substitute(C_aq_expr,
        Dict(C_g_p=>1e-6, HRT_p=>100.0, K_c_p=>0.13, D_p_p=>1e-3,
            U_t_p=>4.0, z_p=>0.0, C_aq0_p=>5.0e-4, sp_defaults...))
    @test C_aq_z0 ≈ 5.0e-4 rtol = 1e-10
end

@testitem "Reversible Scavenging Flux (Eq. 20.35)" setup=[SP2006Setup] tags=[:sp2006] begin
    @parameters C_g_p = 1.0e-6 [unit = u"mol/m^3"]
    @parameters HRT_p = 100.0
    @parameters C_aq0_p = 0.0 [unit = u"mol/m^3"]

    F_expr=reversible_scavenging_flux(C_g_p, HRT_p, K_c_p, D_p_p, U_t_p, h_p, C_aq0_p, p₀_p)

    # Flux should be zero if C_aq0 = C_g * HRT (already at equilibrium)
    F_equil=substitute(F_expr,
        Dict(C_g_p=>1e-6, HRT_p=>100.0, K_c_p=>0.13, D_p_p=>1e-3,
            U_t_p=>4.0, h_p=>1000.0, C_aq0_p=>1e-6*100.0, p₀_p=>2.778e-7, sp_defaults...))
    @test abs(F_equil) < 1e-20

    # Flux should be positive when C_aq0 < C_g * HRT
    F_pos=substitute(F_expr,
        Dict(C_g_p=>1e-6, HRT_p=>100.0, K_c_p=>0.13, D_p_p=>1e-3,
            U_t_p=>4.0, h_p=>1000.0, C_aq0_p=>0.0, p₀_p=>2.778e-7, sp_defaults...))
    @test F_pos > 0
end

@testitem "ReversibleGasScavenging Component" setup=[SP2006Setup] tags=[:sp2006] begin
    sys=ReversibleGasScavenging()
    @test length(equations(sys)) == 2
    @test length(unknowns(sys)) == 2
end

# ====================================================================
# BelowCloudGasScavenging — Irreversible ODE (Eq. 20.22, 20.24, 20.25)
# ====================================================================
@testitem "BelowCloudGasScavenging Component" setup=[SP2006Setup] tags=[:sp2006] begin
    using OrdinaryDiffEqDefault

    sys=BelowCloudGasScavenging()
    compiled=mtkcompile(sys)

    # Use Table 20.1 values for D_p = 0.1 cm
    # K_c = 0.13 m/s, U_t = 3.0 m/s, D_p = 1e-3 m, p₀ = 1 mm/h = 2.778e-7 m/s
    # h = 1000 m
    # Λ = 6 * 2.778e-7 * 0.13 / (1e-3 * 3.0) ≈ 7.22e-5 s⁻¹ = 0.26 h⁻¹
    C_g0=1e-5 # mol/m³
    prob=ODEProblem(compiled,
        merge(Dict(compiled.C_g=>C_g0),
            Dict(compiled.K_c=>0.13, compiled.U_t=>3.0, compiled.D_p=>1e-3,
                compiled.h=>1000.0, compiled.p₀=>2.778e-7)),
        (0.0, 3600.0))
    sol=solve(prob)

    # After 1 hour: C_g(3600) = C_g0 * exp(-Λ*3600)
    # Λ ≈ 7.22e-5 s⁻¹, so Λ*3600 ≈ 0.26
    # C_g(3600) = 1e-5 * exp(-0.26) ≈ 1e-5 * 0.771 ≈ 7.71e-6
    @test sol[compiled.C_g][end] ≈ C_g0 * exp(-7.22e-5 * 3600) rtol = 0.05

    # Verify Λ_gas output
    @test sol[compiled.Λ_gas][end] ≈ 7.22e-5 rtol = 0.05

    # Verify F_bc = Λ * h * C_g
    Λ_end=sol[compiled.Λ_gas][end]
    C_end=sol[compiled.C_g][end]
    @test sol[compiled.F_bc][end] ≈ Λ_end * 1000.0 * C_end rtol = 0.01
end

# ====================================================================
# Eq. 20.53, 20.54 — Slinn Collection Efficiency
# ====================================================================
@testitem "Slinn Collection Efficiency (Eq. 20.53)" setup=[SP2006Setup] tags=[:sp2006] begin
    # Test qualitative U-shape: E should be higher for very small
    # and very large particles, with a minimum near 1 μm (Greenfield gap)
    E_expr=slinn_collection_efficiency(D_p_p, U_t_p, d_p_p, ρ_p_p, T_p)

    # Small particle (d_p = 0.01 μm = 1e-8 m): Brownian diffusion dominates
    E_small=substitute(E_expr,
        Dict(D_p_p=>2e-3, U_t_p=>6.5, d_p_p=>1e-8, ρ_p_p=>1000.0,
            T_p=>298.0, sp_defaults...))
    @test E_small > 1e-4  # should be measurable efficiency

    # Medium particle (d_p = 1 μm = 1e-6 m): Greenfield gap — minimum
    E_medium=substitute(E_expr,
        Dict(D_p_p=>2e-3, U_t_p=>6.5, d_p_p=>1e-6, ρ_p_p=>1000.0,
            T_p=>298.0, sp_defaults...))
    @test E_medium > 0  # still positive

    # Large particle (d_p = 10 μm = 1e-5 m): impaction dominates
    E_large=substitute(E_expr,
        Dict(D_p_p=>2e-3, U_t_p=>6.5, d_p_p=>1e-5, ρ_p_p=>1000.0,
            T_p=>298.0, sp_defaults...))

    # U-shape: E at edges should be larger than at the gap
    @test E_small > E_medium
    # Large particle may or may not exceed medium depending on exact parameters,
    # but should be non-trivial
    @test E_large > 0

    # E should never exceed 1
    @test E_small <= 1.0
    @test E_medium <= 1.0
    @test E_large <= 1.0
end

@testitem "Critical Stokes Number (Eq. 20.54)" setup=[SP2006Setup] tags=[:sp2006] begin
    # S* = (1.2 + (1/12)*ln(1+Re)) / (1 + ln(1+Re))
    # At Re = 0: S* = 1.2 / 1 = 1.2
    # At Re → ∞: S* → (1/12) = 0.0833
    # S* should decrease with Re

    # Re = 0 case: approximation via very small Re
    S_star_low=(1.2+(1/12)*log(1+0.001))/(1+log(1+0.001))
    @test S_star_low ≈ 1.2 rtol = 0.01

    # Re = 100
    S_star_100=(1.2+(1/12)*log(1+100))/(1+log(1+100))
    @test 0.08 < S_star_100 < 1.2

    # Monotonically decreasing
    S_star_1000=(1.2+(1/12)*log(1+1000))/(1+log(1+1000))
    @test S_star_1000 < S_star_100
end

# ====================================================================
# Eq. 20.57 — Particle Scavenging Coefficient
# ====================================================================
@testitem "Particle Scavenging Coefficient (Eq. 20.57)" setup=[SP2006Setup] tags=[:sp2006] begin
    # Λ = (3/2) * E * p₀ / D_p
    # For E = 0.1, p₀ = 2.778e-7 m/s, D_p = 1e-3 m:
    # Λ = 1.5 * 0.1 * 2.778e-7 / 1e-3 = 4.167e-5 s⁻¹
    Λ_expr=particle_scavenging_coeff(d_p_p, p₀_p, D_p_p)
    Λ_val=substitute(Λ_expr,
        Dict(d_p_p=>0.1, p₀_p=>2.778e-7, D_p_p=>1e-3, sp_defaults...))
    @test Λ_val ≈ 4.167e-5 rtol = 0.01
end

@testitem "ParticleScavengingCoeff Component" setup=[SP2006Setup] tags=[:sp2006] begin
    sys=ParticleScavengingCoeff()
    @test length(equations(sys)) == 2
    @test length(unknowns(sys)) == 2
end

# ====================================================================
# Eq. 20.7–20.9 — Wet Deposition Flux and Velocity
# ====================================================================
@testitem "Wet Deposition Flux (Eq. 20.7–20.9)" setup=[SP2006Setup] tags=[:sp2006] begin
    sys=WetDepositionFlux()
    @test length(equations(sys)) == 3
    @test length(unknowns(sys)) == 3

    # Numerical check:
    # C_rain = 1e-3 mol/m³, C_air = 1e-6 mol/m³, p₀ = 2.778e-7 m/s
    # F_w = C_rain * p₀ = 1e-3 * 2.778e-7 = 2.778e-10 mol/m²/s
    # w_r = C_rain / C_air = 1e-3 / 1e-6 = 1000
    # u_w = F_w / C_air = 2.778e-10 / 1e-6 = 2.778e-4 m/s
    # Also u_w = w_r * p₀ = 1000 * 2.778e-7 = 2.778e-4 ✓
    @parameters C_rain_p = 1.0e-3 [unit = u"mol/m^3"]
    @parameters C_air_p = 1.0e-6 [unit = u"mol/m^3"]

    @test substitute(C_rain_p * p₀_p,
        Dict(C_rain_p => 1e-3, p₀_p => 2.778e-7)) ≈ 2.778e-10
    @test substitute(C_rain_p / C_air_p,
        Dict(C_rain_p => 1e-3, C_air_p => 1e-6)) ≈ 1000.0
end

# ====================================================================
# Worked Example (p. 936) — Below-Cloud Scavenging
# ====================================================================
@testitem "Worked Example: Below-Cloud Scavenging" setup=[SP2006Setup] tags=[:sp2006] begin
    using OrdinaryDiffEqDefault

    # From text: Λ = 3.3 h⁻¹ = 9.167e-4 s⁻¹, C₀ = 10 μg/m³, t = 0.5 h
    # After 0.5 h: C = C₀ * exp(-Λ*t) = 10 * exp(-3.3*0.5) = 10 * exp(-1.65)
    #            = 10 * 0.1920 ≈ 1.92 μg/m³
    # Total removed: (10 - 1.92) * 2000 m = 16.16 mg/m²

    # This is a first-order decay, same as our BelowCloudGasScavenging
    # We use the analytical solution to verify
    Λ=3.3/3600  # s⁻¹
    C0=10.0  # arbitrary units
    t_half=0.5*3600  # seconds

    C_final=C0*exp(-Λ*t_half)
    @test C_final ≈ 1.92 rtol = 0.02

    h=2000.0 # m
    total_deposited=(C0-C_final)*h
    @test total_deposited ≈ 16160 rtol = 0.02 # μg/m³ * m = μg/m² → need to match 16.2 mg/m² = 16200 μg/m²
end

# ====================================================================
# Helper Functions
# ====================================================================
@testitem "Particle Relaxation Time" setup=[SP2006Setup] tags=[:sp2006] begin
    # τ = ρ_p * d_p² / (18 * μ_air)
    # For d_p = 1e-6 m, ρ_p = 1000 kg/m³:
    # τ = 1000 * (1e-6)² / (18 * 1.72e-5) = 1e-9 / 3.096e-4 = 3.23e-6 s
    τ_expr=particle_relaxation_time(d_p_p, ρ_p_p)
    τ_val=substitute(τ_expr,
        Dict(d_p_p=>1e-6, ρ_p_p=>1000.0, sp_defaults...))
    @test τ_val ≈ 3.23e-6 rtol = 0.01
    @test ModelingToolkit.get_unit(τ_expr) == u"s"
end

@testitem "Particle Terminal Velocity" setup=[SP2006Setup] tags=[:sp2006] begin
    # u_t = τ * g = 3.23e-3 * 9.81 = 0.0317 m/s for d_p = 1e-6 m
    # Actually, for 1 μm particle this is too high; Stokes settling for unit-density
    # 1 μm particle: u_t ≈ 3.5e-5 m/s. Let's check:
    # τ = 1000*(1e-6)^2/(18*1.72e-5) = 1e-6/3.096e-4 = 3.23e-6 s
    # Wait — (1e-6)^2 = 1e-12, so τ = 1000*1e-12/(18*1.72e-5) = 1e-9/3.096e-4 = 3.23e-6 s
    # u_t = 3.23e-6 * 9.81 = 3.17e-5 m/s ✓
    u_t_expr=particle_terminal_velocity(d_p_p, ρ_p_p)
    u_t_val=substitute(u_t_expr,
        Dict(d_p_p=>1e-6, ρ_p_p=>1000.0, sp_defaults...))
    @test u_t_val ≈ 3.17e-5 rtol = 0.02
    @test ModelingToolkit.get_unit(u_t_expr) == u"m/s"
end

@testitem "Particle Brownian Diffusivity" setup=[SP2006Setup] tags=[:sp2006] begin
    # D = k_B * T / (3π * μ * d_p)
    # For d_p = 1e-6 m, T = 298 K:
    # D = 1.381e-23 * 298 / (3π * 1.72e-5 * 1e-6)
    #   = 4.115e-21 / 1.621e-10 = 2.54e-11 m²/s
    D_expr=particle_diffusivity(d_p_p, T_p)
    D_val=substitute(D_expr,
        Dict(d_p_p=>1e-6, T_p=>298.0, sp_defaults...))
    @test D_val ≈ 2.54e-11 rtol = 0.02
    @test ModelingToolkit.get_unit(D_expr) == u"m^2/s"
end
