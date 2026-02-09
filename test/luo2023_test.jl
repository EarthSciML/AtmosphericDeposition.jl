@testsnippet Luo2023Setup begin
    using AtmosphericDeposition
    using AtmosphericDeposition: turbulence_velocity, cloudy_air_refreshing_rate,
                                 grid_refreshing_time, air_refreshing_limited_rate,
                                 hno3_uptake_efficiency,
                                 cloud_ice_uptake_rate,
                                 air_refreshing_limited_ice_uptake_rate,
                                 cold_cloud_rainout_efficiency,
                                 T_upper_luo, T_lower_luo, γ_base, γ_delta, zero_dimless,
                                 one_dimless,
                                 kinetic_prefactor_si, M_ref, T_ref
    using Test, DynamicQuantities, ModelingToolkit
    using ModelingToolkit: t

    # Helper to get a Float64 from a substituted symbolic expression
    # This resolves @constants that substitute alone doesn't evaluate
    to_float(x) = Float64(ModelingToolkit.Symbolics.value(x))

    # Common parameters for testing
    @parameters begin
        f, [description = "Cloud fraction (dimensionless)"]
        Rᵢ, [unit = u"s^-1", description = "Removing rate"]
        TKE, [unit = u"m^2/s^2", description = "Turbulence kinetic energy"]
        Δx, [unit = u"m", description = "Grid spacing x"]
        Δy, [unit = u"m", description = "Grid spacing y"]
        Δz, [unit = u"m", description = "Grid spacing z"]
        Δt, [unit = u"s", description = "Time step"]
        N_I, [unit = u"m^-3", description = "Ice number concentration"]
        S_I, [unit = u"m^2", description = "Ice surface area"]
        r_ice, [unit = u"m", description = "Ice crystal radius"]
        D_g, [unit = u"m^2/s", description = "Gas diffusion coefficient"]
        M, [unit = u"g/mol", description = "Molar mass"]
        T, [unit = u"K", description = "Temperature"]
    end

    # Constant defaults for full substitution
    const_defaults = Dict(
        T_upper_luo => 220.0, T_lower_luo => 209.0,
        γ_base => 3e-3, γ_delta => 4e-3,
        zero_dimless => 0.0, one_dimless => 1.0,
        kinetic_prefactor_si => 2.749064e-2,
        M_ref => 1.0, T_ref => 1.0
    )
end

# ============================================================================
# Structural Tests
# ============================================================================
@testitem "Structural - AirRefreshingLimitation" setup=[Luo2023Setup] tags=[:luo2023] begin
    sys=AirRefreshingLimitation()
    @test length(equations(sys)) == 3
    @test length(unknowns(sys)) == 3

    vars=unknowns(sys)
    var_names=[string(ModelingToolkit.Symbolics.tosymbol(v, escape = false)) for v in vars]
    @test "Kᵢ" in var_names
    @test "τ_A" in var_names
    @test "R_A" in var_names
end

@testitem "Structural - CloudIceUptakeLimitation" setup=[Luo2023Setup] tags=[:luo2023] begin
    sys=CloudIceUptakeLimitation()
    @test length(equations(sys)) == 5
    @test length(unknowns(sys)) == 5

    vars=unknowns(sys)
    var_names=[string(ModelingToolkit.Symbolics.tosymbol(v, escape = false)) for v in vars]
    @test "γ" in var_names
    @test "R_U" in var_names
    @test "Kᵢ" in var_names
    @test "R_AU" in var_names
    @test "F_I" in var_names
end

@testitem "Structural - WetScavengingLimitations" setup=[Luo2023Setup] tags=[:luo2023] begin
    sys=WetScavengingLimitations()
    @test length(equations(sys)) == 7
    @test length(unknowns(sys)) == 7

    vars=unknowns(sys)
    var_names=[string(ModelingToolkit.Symbolics.tosymbol(v, escape = false)) for v in vars]
    @test "Kᵢ" in var_names
    @test "τ_A" in var_names
    @test "R_A" in var_names
    @test "γ" in var_names
    @test "R_U" in var_names
    @test "R_AU" in var_names
    @test "F_I" in var_names
end

# ============================================================================
# Unit Tests
# ============================================================================
@testitem "Units - AirRefreshingLimitation" setup=[Luo2023Setup] tags=[:luo2023] begin
    sys=AirRefreshingLimitation()
    vars=unknowns(sys)
    for v in vars
        name=string(ModelingToolkit.Symbolics.tosymbol(v, escape = false))
        if name=="Kᵢ"||name=="R_A"
            @test ModelingToolkit.get_unit(v) == u"s^-1"
        elseif name=="τ_A"
            @test ModelingToolkit.get_unit(v) == u"s"
        end
    end
end

@testitem "Units - CloudIceUptakeLimitation" setup=[Luo2023Setup] tags=[:luo2023] begin
    sys=CloudIceUptakeLimitation()
    vars=unknowns(sys)
    for v in vars
        name=string(ModelingToolkit.Symbolics.tosymbol(v, escape = false))
        if name in ["R_U", "Kᵢ", "R_AU"]
            @test ModelingToolkit.get_unit(v) == u"s^-1"
        end
    end
end

@testitem "Units - WetScavengingLimitations" setup=[Luo2023Setup] tags=[:luo2023] begin
    sys=WetScavengingLimitations()
    vars=unknowns(sys)
    for v in vars
        name=string(ModelingToolkit.Symbolics.tosymbol(v, escape = false))
        if name in ["Kᵢ", "R_A", "R_U", "R_AU"]
            @test ModelingToolkit.get_unit(v) == u"s^-1"
        elseif name == "τ_A"
            @test ModelingToolkit.get_unit(v) == u"s"
        end
    end
end

# ============================================================================
# Compilation Tests
# ============================================================================
@testitem "Compilation - AirRefreshingLimitation" setup=[Luo2023Setup] tags=[:luo2023] begin
    sys = AirRefreshingLimitation()
    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.System
end

@testitem "Compilation - CloudIceUptakeLimitation" setup=[Luo2023Setup] tags=[:luo2023] begin
    sys = CloudIceUptakeLimitation()
    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.System
end

@testitem "Compilation - WetScavengingLimitations" setup=[Luo2023Setup] tags=[:luo2023] begin
    sys = WetScavengingLimitations()
    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.System
end

# ============================================================================
# Equation Verification Tests
# ============================================================================
@testitem "Eq. 10 - Turbulence velocity" setup=[Luo2023Setup] tags=[:luo2023] begin
    # u' = √(2/3 · TKE)
    # For TKE = 1.5 m²/s²: u' = √(1.0) = 1.0 m/s
    v=turbulence_velocity(TKE)
    result=to_float(substitute(v, Dict(TKE=>1.5)))
    @test result ≈ 1.0

    # For TKE = 6.0 m²/s²: u' = √(4.0) = 2.0 m/s
    result2=to_float(substitute(v, Dict(TKE=>6.0)))
    @test result2 ≈ 2.0
end

@testitem "Eq. 11 - Cloudy air refreshing rate" setup=[Luo2023Setup] tags=[:luo2023] begin
    # Kᵢ = √(2/3·TKE) · [1/(√f·Δx) + 1/(√f·Δy) + 1/Δz]
    Ki=cloudy_air_refreshing_rate(f, TKE, Δx, Δy, Δz)

    # Test case: f=0.25, TKE=1.5, Δx=Δy=200000, Δz=1000
    # v_turb = √(1.0) = 1.0
    # 1/(0.5*200000) + 1/(0.5*200000) + 1/1000 = 1e-5 + 1e-5 + 1e-3 = 1.02e-3
    # Kᵢ = 1.0 * 1.02e-3 = 1.02e-3
    vals=Dict(f=>0.25, TKE=>1.5, Δx=>200000.0, Δy=>200000.0, Δz=>1000.0)
    result=to_float(substitute(Ki, vals))
    expected=1.0*(1/(0.5*200000)+1/(0.5*200000)+1/1000.0)
    @test result ≈ expected rtol = 1e-10
end

@testitem "Eq. 5 - Grid refreshing time" setup=[Luo2023Setup] tags=[:luo2023] begin
    # τ_A = (1 - f) / (f · Kᵢ)
    @parameters Kᵢ_param [unit = u"s^-1"]
    τ=grid_refreshing_time(f, Kᵢ_param)

    # f=0.3, Kᵢ=0.001: τ = 0.7 / (0.3 * 0.001) = 2333.33
    result=to_float(substitute(τ, Dict(f=>0.3, Kᵢ_param=>0.001)))
    @test result ≈ 0.7 / (0.3 * 0.001) rtol = 1e-10
end

@testitem "Eq. 2 - Air refreshing limited rate" setup=[Luo2023Setup] tags=[:luo2023] begin
    # R_A = 1 / (1/(f·Rᵢ) + τ_A)
    @parameters τ_A_param [unit = u"s"]
    R_A=air_refreshing_limited_rate(f, Rᵢ, τ_A_param)

    # f=0.5, Rᵢ=0.001, τ_A=100
    # R_A = 1 / (1/(0.5*0.001) + 100) = 1 / (2000 + 100) = 1/2100
    result=to_float(substitute(R_A, Dict(f=>0.5, Rᵢ=>0.001, τ_A_param=>100.0)))
    @test result ≈ 1.0 / 2100.0 rtol = 1e-10

    # When τ_A = 0, R_A = f·Rᵢ
    result0=to_float(substitute(R_A, Dict(f=>0.5, Rᵢ=>0.001, τ_A_param=>0.0)))
    @test result0 ≈ 0.5 * 0.001 rtol = 1e-10
end

@testitem "Eq. 15 - HNO3 uptake efficiency" setup=[Luo2023Setup] tags=[:luo2023] begin
    γ_val=hno3_uptake_efficiency(T)

    # T ≥ 220 K: γ = 0.003
    @test to_float(substitute(γ_val, merge(const_defaults, Dict(T => 220.0)))) ≈ 0.003
    @test to_float(substitute(γ_val, merge(const_defaults, Dict(T => 250.0)))) ≈ 0.003

    # T ≤ 209 K: γ = 0.007
    @test to_float(substitute(γ_val, merge(const_defaults, Dict(T => 209.0)))) ≈ 0.007
    @test to_float(substitute(γ_val, merge(const_defaults, Dict(T => 200.0)))) ≈ 0.007

    # T = 214.5 K (midpoint): γ = 0.003 + 0.004 * 0.5 = 0.005
    @test to_float(substitute(γ_val, merge(const_defaults, Dict(T => 214.5)))) ≈ 0.005 rtol = 1e-10
end

@testitem "Eq. 14 - Cloud ice uptake rate" setup=[Luo2023Setup] tags=[:luo2023] begin
    @parameters γ_param
    R_U_expr=cloud_ice_uptake_rate(N_I, S_I, r_ice, D_g, M, T, γ_param)

    # Test with typical values:
    # N_I = 1e6 m⁻³, S_I = 1e-8 m², r = 1e-4 m, D_g = 1e-5 m²/s,
    # M = 63.0 g/mol (HNO₃), T = 210 K, γ = 0.005
    vals=merge(const_defaults,
        Dict(
            N_I=>1e6, S_I=>1e-8, r_ice=>1e-4, D_g=>1e-5,
            M=>63.0, T=>210.0, γ_param=>0.005
        ))
    result=to_float(substitute(R_U_expr, vals))

    # Manual calculation:
    # diffusion_term = r/D_g = 1e-4 / 1e-5 = 10 s/m
    # kinetic_term = 2.749064e-2 * √(63/1) / (0.005 * √(210/1))
    #             = 2.749064e-2 * 7.9373 / (0.005 * 14.4914)
    #             = 0.21822 / 0.072457 = 3.0117 s/m
    # R_U = 1e6 * 1e-8 / (10 + 3.0117) = 0.01 / 13.0117 = 7.688e-4 s⁻¹
    diffusion_term=1e-4/1e-5
    kinetic_term=2.749064e-2*sqrt(63.0)/(0.005*sqrt(210.0))
    expected=1e6*1e-8/(diffusion_term+kinetic_term)
    @test result ≈ expected rtol = 1e-6
end

@testitem "Eq. 13 - Air refreshing limited ice uptake rate" setup=[Luo2023Setup] tags=[:luo2023] begin
    # R_{A,U} = f · R_U · Kᵢ / (Kᵢ + (1-f) · R_U)
    @parameters R_U_param [unit = u"s^-1"]
    @parameters Kᵢ_param [unit = u"s^-1"]
    R_AU=air_refreshing_limited_ice_uptake_rate(f, R_U_param, Kᵢ_param)

    # f=0.3, R_U=0.001, Kᵢ=0.01
    vals=Dict(f=>0.3, R_U_param=>0.001, Kᵢ_param=>0.01)
    result=to_float(substitute(R_AU, vals))
    expected=0.3*0.001*0.01/(0.01+0.7*0.001)
    @test result ≈ expected rtol = 1e-10
end

@testitem "Eq. 12 - Cold cloud rainout efficiency" setup=[Luo2023Setup] tags=[:luo2023] begin
    # F_I = 1 - exp(-R_{A,U} · Δt)
    @parameters R_AU_param [unit = u"s^-1"]
    F_I=cold_cloud_rainout_efficiency(R_AU_param, Δt)

    # R_AU=0.001, Δt=600: F_I = 1 - exp(-0.6) ≈ 0.4512
    result=to_float(substitute(F_I, Dict(R_AU_param=>0.001, Δt=>600.0)))
    @test result ≈ 1 - exp(-0.6) rtol = 1e-10
end

# ============================================================================
# Limiting Behavior Tests
# ============================================================================
@testitem "Limiting behavior - strong mixing" setup=[Luo2023Setup] tags=[:luo2023] begin
    # When TKE is very large, τ_A → 0, so R_A → f·Rᵢ
    R_A_expr=air_refreshing_limited_rate(f, Rᵢ,
        grid_refreshing_time(f, cloudy_air_refreshing_rate(f, TKE, Δx, Δy, Δz)))

    # Very large TKE = 1000 m²/s², small grid for extra strong mixing
    vals=Dict(f=>0.3, Rᵢ=>0.001, TKE=>1000.0,
        Δx=>50000.0, Δy=>50000.0, Δz=>500.0)
    R_A_val=to_float(substitute(R_A_expr, vals))
    f_Ri=0.3*0.001

    # With large TKE and small grid, R_A should be close to f·Rᵢ
    @test R_A_val ≈ f_Ri rtol = 0.05
end

@testitem "Limiting behavior - weak mixing" setup=[Luo2023Setup] tags=[:luo2023] begin
    # When TKE is very small, τ_A → large, so R_A → 0
    R_A_expr=air_refreshing_limited_rate(f, Rᵢ,
        grid_refreshing_time(f, cloudy_air_refreshing_rate(f, TKE, Δx, Δy, Δz)))

    vals=Dict(f=>0.3, Rᵢ=>0.001, TKE=>1e-6,
        Δx=>200000.0, Δy=>200000.0, Δz=>1000.0)
    R_A_val=to_float(substitute(R_A_expr, vals))
    f_Ri=0.3*0.001

    # With tiny TKE, R_A should be much smaller than f·Rᵢ
    @test R_A_val < 0.1 * f_Ri
end

@testitem "Limiting behavior - full cloud coverage" setup=[Luo2023Setup] tags=[:luo2023] begin
    # When f → 1, τ_A → 0 (Eq. 5: (1-f)/(f·Kᵢ) → 0), so R_A → Rᵢ
    R_A_expr=air_refreshing_limited_rate(f, Rᵢ,
        grid_refreshing_time(f, cloudy_air_refreshing_rate(f, TKE, Δx, Δy, Δz)))

    vals=Dict(f=>0.999, Rᵢ=>0.001, TKE=>1.0,
        Δx=>200000.0, Δy=>200000.0, Δz=>1000.0)
    R_A_val=to_float(substitute(R_A_expr, vals))

    @test R_A_val ≈ 0.001 rtol = 0.01
end

@testitem "Limiting behavior - F_I limits" setup=[Luo2023Setup] tags=[:luo2023] begin
    @parameters R_AU_param [unit = u"s^-1"]

    # F_I → 0 when R_AU·Δt → 0
    F_I_expr=cold_cloud_rainout_efficiency(R_AU_param, Δt)
    @test to_float(substitute(F_I_expr, Dict(R_AU_param => 1e-10, Δt => 600.0))) ≈ 0.0 atol = 1e-6

    # F_I → 1 when R_AU·Δt → large
    @test to_float(substitute(F_I_expr, Dict(R_AU_param => 1.0, Δt => 6000.0))) ≈ 1.0 atol = 1e-6
end

# ============================================================================
# Qualitative Behavior Tests
# ============================================================================
@testitem "Qualitative - R_A monotonic in TKE" setup=[Luo2023Setup] tags=[:luo2023] begin
    # R_A should increase with TKE (stronger turbulence → more mixing → more scavenging)
    R_A_expr=air_refreshing_limited_rate(f, Rᵢ,
        grid_refreshing_time(f, cloudy_air_refreshing_rate(f, TKE, Δx, Δy, Δz)))

    base=Dict(f=>0.3, Rᵢ=>0.001, Δx=>200000.0, Δy=>200000.0, Δz=>1000.0)

    tke_vals=[0.01, 0.1, 1.0, 10.0, 100.0]
    R_A_vals=[to_float(substitute(R_A_expr, merge(base, Dict(TKE=>tke))))
              for tke in tke_vals]

    # Check monotonically increasing
    for i in 2:length(R_A_vals)
        @test R_A_vals[i] > R_A_vals[i - 1]
    end
end

@testitem "Qualitative - R_A bounded" setup=[Luo2023Setup] tags=[:luo2023] begin
    # R_A should always be ≤ f·Rᵢ and ≥ 0
    R_A_expr=air_refreshing_limited_rate(f, Rᵢ,
        grid_refreshing_time(f, cloudy_air_refreshing_rate(f, TKE, Δx, Δy, Δz)))

    base=Dict(Δx=>200000.0, Δy=>200000.0, Δz=>1000.0)
    for f_val in [0.1, 0.3, 0.5, 0.8]
        for Ri_val in [0.0001, 0.001, 0.01]
            for tke_val in [0.1, 1.0, 10.0]
                vals=merge(base, Dict(f=>f_val, Rᵢ=>Ri_val, TKE=>tke_val))
                R_A_val=to_float(substitute(R_A_expr, vals))
                @test R_A_val >= 0
                @test R_A_val <= f_val * Ri_val * 1.0001  # small tolerance
            end
        end
    end
end

@testitem "Qualitative - F_I monotonic in Δt" setup=[Luo2023Setup] tags=[:luo2023] begin
    # F_I should increase with Δt
    @parameters R_AU_param [unit = u"s^-1"]
    F_I_expr=cold_cloud_rainout_efficiency(R_AU_param, Δt)

    dt_vals=[60.0, 300.0, 600.0, 1200.0, 3600.0]
    F_I_vals=[to_float(substitute(F_I_expr, Dict(R_AU_param=>0.001, Δt=>dt_val)))
              for dt_val in dt_vals]

    for i in 2:length(F_I_vals)
        @test F_I_vals[i] > F_I_vals[i - 1]
    end
end

@testitem "Qualitative - γ bounded" setup=[Luo2023Setup] tags=[:luo2023] begin
    # γ should be between 0.003 and 0.007 for all temperatures
    γ_expr=hno3_uptake_efficiency(T)

    for T_val in [150.0, 180.0, 200.0, 209.0, 214.5, 220.0, 250.0, 300.0]
        γ_val=to_float(substitute(γ_expr, merge(const_defaults, Dict(T=>T_val))))
        @test γ_val >= 0.003 - 1e-10
        @test γ_val <= 0.007 + 1e-10
    end
end
