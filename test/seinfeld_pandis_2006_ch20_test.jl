@testsnippet SP2006Setup begin
    using AtmosphericDeposition
    using AtmosphericDeposition: mass_transfer_coeff, gas_scavenging_coeff,
        reversible_drop_conc, reversible_scavenging_flux,
        slinn_collection_efficiency, particle_scavenging_coeff,
        particle_relaxation_time, particle_terminal_velocity,
        particle_diffusivity,
        μ_air_sp, ρ_air_sp, μ_w_sp, ρ_w_sp, g_sp, kB_sp
    using Test, DynamicQuantities, ModelingToolkit
    using Symbolics
    using Symbolics: value

    function to_float(x)
        v = value(x)
        v isa Number && return Float64(v)
        return Float64(eval(Symbolics.toexpr(Symbolics.unwrap(x))))
    end

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
        kB_sp => 1.3806488e-23,
    ]
end

# ====================================================================
# Eq. 20.12 — Mass Transfer Coefficient
# ====================================================================
@testitem "Mass Transfer Coefficient (Eq. 20.12)" setup = [SP2006Setup] tags = [:sp2006] begin
    Kc_expr = mass_transfer_coeff(D_g_p, D_p_p, U_t_p)
    Kc_val = to_float(
        substitute(
            Kc_expr,
            Dict(D_g_p => 1.26e-5, D_p_p => 1.0e-4, U_t_p => 0.26, sp_defaults...)
        )
    )
    @test Kc_val ≈ 0.3586 rtol = 0.01

    # Check units
    @test ModelingToolkit.get_unit(Kc_expr) == u"m/s"

    Kc_val2 = to_float(
        substitute(
            Kc_expr,
            Dict(D_g_p => 1.26e-5, D_p_p => 0.01, U_t_p => 10.0, sp_defaults...)
        )
    )
    @test Kc_val2 ≈ 0.0687 rtol = 0.01
end

@testitem "MassTransferCoeff Component" setup = [SP2006Setup] tags = [:sp2006] begin
    sys = MassTransferCoeff()
    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 1
end

# ====================================================================
# Eq. 20.25 — Gas Scavenging Coefficient
# ====================================================================
@testitem "Gas Scavenging Coefficient (Eq. 20.25)" setup = [SP2006Setup] tags = [:sp2006] begin
    Λ_expr = gas_scavenging_coeff(K_c_p, U_t_p, D_p_p, p₀_p)
    Λ_val = to_float(
        substitute(
            Λ_expr,
            Dict(K_c_p => 0.13, U_t_p => 3.0, D_p_p => 1.0e-3, p₀_p => 2.778e-7, sp_defaults...)
        )
    )
    @test Λ_val ≈ 7.22e-5 rtol = 0.02

    Λ_val2 = to_float(
        substitute(
            Λ_expr,
            Dict(K_c_p => 0.32, U_t_p => 0.26, D_p_p => 1.0e-4, p₀_p => 2.778e-7, sp_defaults...)
        )
    )
    @test Λ_val2 ≈ 0.0205 rtol = 0.02

    # Check units
    @test ModelingToolkit.get_unit(Λ_expr) == u"s^-1"
end

@testitem "GasScavengingCoeff Component" setup = [SP2006Setup] tags = [:sp2006] begin
    sys = GasScavengingCoeff()
    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 1
end

# ====================================================================
# Eq. 20.28, 20.35 — Reversible Gas Scavenging
# ====================================================================
@testitem "Reversible Gas Scavenging (Eq. 20.28)" setup = [SP2006Setup] tags = [:sp2006] begin
    @parameters C_g_p = 1.0e-6 [unit = u"mol/m^3"]
    @parameters HRT_p = 100.0
    @parameters C_aq0_p = 0.0 [unit = u"mol/m^3"]
    @parameters z_p = 1000.0 [unit = u"m"]

    C_aq_expr = reversible_drop_conc(C_g_p, HRT_p, K_c_p, D_p_p, U_t_p, z_p, C_aq0_p)
    C_aq_huge_H = to_float(
        substitute(
            C_aq_expr,
            Dict(
                C_g_p => 1.0e-6, HRT_p => 1.0e12, K_c_p => 0.13, D_p_p => 1.0e-3,
                U_t_p => 4.0, z_p => 1000.0, C_aq0_p => 0.0, sp_defaults...
            )
        )
    )
    @test C_aq_huge_H ≈ 0.0 atol = 1.0

    C_aq_equil = to_float(
        substitute(
            C_aq_expr,
            Dict(
                C_g_p => 1.0e-6, HRT_p => 1.0, K_c_p => 10.0, D_p_p => 1.0e-4,
                U_t_p => 0.1, z_p => 10000.0, C_aq0_p => 0.0, sp_defaults...
            )
        )
    )
    @test C_aq_equil ≈ 1.0e-6 rtol = 1.0e-6

    C_aq_z0 = to_float(
        substitute(
            C_aq_expr,
            Dict(
                C_g_p => 1.0e-6, HRT_p => 100.0, K_c_p => 0.13, D_p_p => 1.0e-3,
                U_t_p => 4.0, z_p => 0.0, C_aq0_p => 5.0e-4, sp_defaults...
            )
        )
    )
    @test C_aq_z0 ≈ 5.0e-4 rtol = 1.0e-10
end

@testitem "Reversible Scavenging Flux (Eq. 20.35)" setup = [SP2006Setup] tags = [:sp2006] begin
    @parameters C_g_p = 1.0e-6 [unit = u"mol/m^3"]
    @parameters HRT_p = 100.0
    @parameters C_aq0_p = 0.0 [unit = u"mol/m^3"]

    F_expr = reversible_scavenging_flux(C_g_p, HRT_p, K_c_p, D_p_p, U_t_p, h_p, C_aq0_p, p₀_p)

    F_equil = to_float(
        substitute(
            F_expr,
            Dict(
                C_g_p => 1.0e-6, HRT_p => 100.0, K_c_p => 0.13, D_p_p => 1.0e-3,
                U_t_p => 4.0, h_p => 1000.0, C_aq0_p => 1.0e-6 * 100.0, p₀_p => 2.778e-7, sp_defaults...
            )
        )
    )
    @test abs(F_equil) < 1.0e-20

    F_pos = to_float(
        substitute(
            F_expr,
            Dict(
                C_g_p => 1.0e-6, HRT_p => 100.0, K_c_p => 0.13, D_p_p => 1.0e-3,
                U_t_p => 4.0, h_p => 1000.0, C_aq0_p => 0.0, p₀_p => 2.778e-7, sp_defaults...
            )
        )
    )
    @test F_pos > 0
end

@testitem "ReversibleGasScavenging Component" setup = [SP2006Setup] tags = [:sp2006] begin
    sys = ReversibleGasScavenging()
    @test length(equations(sys)) == 2
    @test length(unknowns(sys)) == 2
end

# ====================================================================
# BelowCloudGasScavenging — Irreversible ODE (Eq. 20.22, 20.24, 20.25)
# ====================================================================
@testitem "BelowCloudGasScavenging Component" setup = [SP2006Setup] tags = [:sp2006] begin
    using OrdinaryDiffEqDefault

    sys = BelowCloudGasScavenging()
    compiled = mtkcompile(sys)

    C_g0 = 1.0e-5 # mol/m³
    prob = ODEProblem(
        compiled,
        merge(
            Dict(compiled.C_g => C_g0),
            Dict(
                compiled.K_c => 0.13, compiled.U_t => 3.0, compiled.D_p => 1.0e-3,
                compiled.h => 1000.0, compiled.p₀ => 2.778e-7
            )
        ),
        (0.0, 3600.0)
    )
    sol = solve(prob)

    @test sol[compiled.C_g][end] ≈ C_g0 * exp(-7.22e-5 * 3600) rtol = 0.05

    @test sol[compiled.Λ_gas][end] ≈ 7.22e-5 rtol = 0.05

    Λ_end = sol[compiled.Λ_gas][end]
    C_end = sol[compiled.C_g][end]
    @test sol[compiled.F_bc][end] ≈ Λ_end * 1000.0 * C_end rtol = 0.01
end

# ====================================================================
# Eq. 20.53, 20.54 — Slinn Collection Efficiency
# ====================================================================
@testitem "Slinn Collection Efficiency (Eq. 20.53)" setup = [SP2006Setup] tags = [:sp2006] begin
    E_expr = slinn_collection_efficiency(D_p_p, U_t_p, d_p_p, ρ_p_p, T_p)

    E_small = to_float(
        substitute(
            E_expr,
            Dict(
                D_p_p => 2.0e-3, U_t_p => 6.5, d_p_p => 1.0e-8, ρ_p_p => 1000.0,
                T_p => 298.0, sp_defaults...
            )
        )
    )
    @test E_small > 1.0e-4

    E_medium = to_float(
        substitute(
            E_expr,
            Dict(
                D_p_p => 2.0e-3, U_t_p => 6.5, d_p_p => 1.0e-6, ρ_p_p => 1000.0,
                T_p => 298.0, sp_defaults...
            )
        )
    )
    @test E_medium > 0

    E_large = to_float(
        substitute(
            E_expr,
            Dict(
                D_p_p => 2.0e-3, U_t_p => 6.5, d_p_p => 1.0e-5, ρ_p_p => 1000.0,
                T_p => 298.0, sp_defaults...
            )
        )
    )

    @test E_small > E_medium
    @test E_large > 0

    @test E_small <= 1.0
    @test E_medium <= 1.0
    @test E_large <= 1.0
end

@testitem "Critical Stokes Number (Eq. 20.54)" setup = [SP2006Setup] tags = [:sp2006] begin
    S_star_low = (1.2 + (1 / 12) * log(1 + 0.001)) / (1 + log(1 + 0.001))
    @test S_star_low ≈ 1.2 rtol = 0.01

    S_star_100 = (1.2 + (1 / 12) * log(1 + 100)) / (1 + log(1 + 100))
    @test 0.08 < S_star_100 < 1.2

    S_star_1000 = (1.2 + (1 / 12) * log(1 + 1000)) / (1 + log(1 + 1000))
    @test S_star_1000 < S_star_100
end

# ====================================================================
# Eq. 20.57 — Particle Scavenging Coefficient
# ====================================================================
@testitem "Particle Scavenging Coefficient (Eq. 20.57)" setup = [SP2006Setup] tags = [:sp2006] begin
    Λ_expr = particle_scavenging_coeff(d_p_p, p₀_p, D_p_p)
    Λ_val = to_float(
        substitute(
            Λ_expr,
            Dict(d_p_p => 0.1, p₀_p => 2.778e-7, D_p_p => 1.0e-3, sp_defaults...)
        )
    )
    @test Λ_val ≈ 4.167e-5 rtol = 0.01
end

@testitem "ParticleScavengingCoeff Component" setup = [SP2006Setup] tags = [:sp2006] begin
    sys = ParticleScavengingCoeff()
    @test length(equations(sys)) == 2
    @test length(unknowns(sys)) == 2
end

# ====================================================================
# Eq. 20.7–20.9 — Wet Deposition Flux and Velocity
# ====================================================================
@testitem "Wet Deposition Flux (Eq. 20.7–20.9)" setup = [SP2006Setup] tags = [:sp2006] begin
    sys = WetDepositionFlux()
    @test length(equations(sys)) == 3
    @test length(unknowns(sys)) == 3

    @parameters C_rain_p = 1.0e-3 [unit = u"mol/m^3"]
    @parameters C_air_p = 1.0e-6 [unit = u"mol/m^3"]

    @test to_float(
        substitute(
            C_rain_p * p₀_p,
            Dict(C_rain_p => 1.0e-3, p₀_p => 2.778e-7)
        )
    ) ≈ 2.778e-10
    @test to_float(
        substitute(
            C_rain_p / C_air_p,
            Dict(C_rain_p => 1.0e-3, C_air_p => 1.0e-6)
        )
    ) ≈ 1000.0
end

# ====================================================================
# Worked Example (p. 936) — Below-Cloud Scavenging
# ====================================================================
@testitem "Worked Example: Below-Cloud Scavenging" setup = [SP2006Setup] tags = [:sp2006] begin
    using OrdinaryDiffEqDefault

    Λ = 3.3 / 3600  # s⁻¹
    C0 = 10.0  # arbitrary units
    t_half = 0.5 * 3600  # seconds

    C_final = C0 * exp(-Λ * t_half)
    @test C_final ≈ 1.92 rtol = 0.02

    h = 2000.0 # m
    total_deposited = (C0 - C_final) * h
    @test total_deposited ≈ 16160 rtol = 0.02
end

# ====================================================================
# Helper Functions
# ====================================================================
@testitem "Particle Relaxation Time" setup = [SP2006Setup] tags = [:sp2006] begin
    τ_expr = particle_relaxation_time(d_p_p, ρ_p_p)
    τ_val = to_float(
        substitute(
            τ_expr,
            Dict(d_p_p => 1.0e-6, ρ_p_p => 1000.0, sp_defaults...)
        )
    )
    @test τ_val ≈ 3.23e-6 rtol = 0.01
    @test ModelingToolkit.get_unit(τ_expr) == u"s"
end

@testitem "Particle Terminal Velocity" setup = [SP2006Setup] tags = [:sp2006] begin
    u_t_expr = particle_terminal_velocity(d_p_p, ρ_p_p)
    u_t_val = to_float(
        substitute(
            u_t_expr,
            Dict(d_p_p => 1.0e-6, ρ_p_p => 1000.0, sp_defaults...)
        )
    )
    @test u_t_val ≈ 3.17e-5 rtol = 0.02
    @test ModelingToolkit.get_unit(u_t_expr) == u"m/s"
end

@testitem "Particle Brownian Diffusivity" setup = [SP2006Setup] tags = [:sp2006] begin
    D_expr = particle_diffusivity(d_p_p, T_p)
    D_val = to_float(
        substitute(
            D_expr,
            Dict(d_p_p => 1.0e-6, T_p => 298.0, sp_defaults...)
        )
    )
    @test D_val ≈ 2.54e-11 rtol = 0.02
    @test ModelingToolkit.get_unit(D_expr) == u"m^2/s"
end
