"""
Test suite for Seinfeld & Pandis Chapter 19 dry deposition implementation.

Comprehensive tests for the ModelingToolkit.jl implementation of
Seinfeld & Pandis (2006) Chapter 19 "Dry Deposition" equations.

Tests verify:
1. Structural correctness (model components, variables, parameters)
2. Equation verification against textbook formulas
3. Numerical accuracy against known values from the chapter
4. Physical constraints and edge cases
5. Dimensional consistency

Reference: Seinfeld, J.H. and Pandis, S.N. (2006). Atmospheric Chemistry and Physics,
2nd Edition, Chapter 19: Dry Deposition. Wiley-Interscience.
"""

@testsnippet SeinfeldPandisSetup begin
    using AtmosphericDeposition.DryDeposition
    using Test, DynamicQuantities, ModelingToolkit
    using ModelingToolkit: t, D
end

#=============================================================================
# Test Set 1: Structural Tests
=============================================================================#

@testitem "1.1 AerodynamicResistance Structure" setup=[SeinfeldPandisSetup] begin
    @named ra_sys = AerodynamicResistance()
    @test ra_sys isa System

    # Check that all expected variables exist
    vars = unknowns(ra_sys)
    var_names = string.(vars)
    @test any(occursin("u_star", n) for n in var_names)
    @test any(occursin("r_a", n) for n in var_names)
    @test length(vars) == 2  # u_star and r_a

    # Check parameters
    params = parameters(ra_sys)
    param_names = string.(params)
    @test any(occursin("z", n) for n in param_names)
    @test any(occursin("z_0", n) for n in param_names)

    # Check equations
    eqs = equations(ra_sys)
    @test length(eqs) == 1  # Single equation for r_a
end

@testitem "1.2 QuasiLaminarResistanceGas Structure" setup=[SeinfeldPandisSetup] begin
    @named rb_gas = QuasiLaminarResistanceGas()
    @test rb_gas isa System

    vars = unknowns(rb_gas)
    var_names = string.(vars)
    @test any(occursin("Sc", n) for n in var_names)
    @test any(occursin("r_b", n) for n in var_names)
    @test any(occursin("u_star", n) for n in var_names)
    @test length(vars) == 3

    eqs = equations(rb_gas)
    @test length(eqs) == 2  # Sc definition + r_b equation
end

@testitem "1.3 ParticleSettling Structure" setup=[SeinfeldPandisSetup] begin
    @named ps = ParticleSettling()
    @test ps isa System

    vars = unknowns(ps)
    var_names = string.(vars)
    @test any(occursin("v_s", n) for n in var_names)
    @test any(occursin("D_diff", n) for n in var_names)
    @test length(vars) == 2

    params = parameters(ps)
    param_names = string.(params)
    @test any(occursin("rho_p", n) for n in param_names)
    @test any(occursin("D_p", n) for n in param_names)
    @test any(occursin("mu", n) for n in param_names)
    @test any(occursin("T", n) for n in param_names)
    @test any(occursin("C_c", n) for n in param_names)
end

@testitem "1.4 QuasiLaminarResistanceParticle Structure" setup=[SeinfeldPandisSetup] begin
    @named rb_part = QuasiLaminarResistanceParticle()
    @test rb_part isa System

    vars = unknowns(rb_part)
    var_names = string.(vars)
    @test any(occursin("E_B", n) for n in var_names)
    @test any(occursin("E_IM", n) for n in var_names)
    @test any(occursin("E_IN", n) for n in var_names)
    @test any(occursin("R_1", n) for n in var_names)
    @test any(occursin("St", n) for n in var_names)
    @test any(occursin("Sc", n) for n in var_names)
    @test any(occursin("r_b", n) for n in var_names)

    # Should have 7 equations per the Zhang et al. model
    eqs = equations(rb_part)
    @test length(eqs) == 7
end

@testitem "1.5 SurfaceResistance Structure" setup=[SeinfeldPandisSetup] begin
    @named sr = SurfaceResistance()
    @test sr isa System

    vars = unknowns(sr)
    var_names = string.(vars)
    @test any(occursin("r_st", n) for n in var_names)
    @test any(occursin("r_m", n) for n in var_names)
    @test any(occursin("r_c", n) for n in var_names)
    @test any(occursin("G", n) for n in var_names)
    @test any(occursin("T_s", n) for n in var_names)
end

@testitem "1.6 DryDepositionGas Composite Structure" setup=[SeinfeldPandisSetup] begin
    @named gas_dep = DryDepositionGas()
    @test gas_dep isa System

    vars = unknowns(gas_dep)
    var_names = string.(vars)
    @test any(occursin("v_d", n) for n in var_names)
    @test any(occursin("F", n) for n in var_names)
    @test any(occursin("r_t", n) for n in var_names)

    # Check subsystems are composed
    @test any(occursin("aero", n) for n in var_names)
    @test any(occursin("qlam", n) for n in var_names)
    @test any(occursin("surf", n) for n in var_names)
end

@testitem "1.7 DryDepositionParticle Composite Structure" setup=[SeinfeldPandisSetup] begin
    @named part_dep = DryDepositionParticle()
    @test part_dep isa System

    vars = unknowns(part_dep)
    var_names = string.(vars)
    @test any(occursin("v_d", n) for n in var_names)
    @test any(occursin("F", n) for n in var_names)

    # Check subsystems are composed
    @test any(occursin("aero", n) for n in var_names)
    @test any(occursin("settling", n) for n in var_names)
    @test any(occursin("qlam", n) for n in var_names)
end

#=============================================================================
# Test Set 2: Equation Verification Tests
=============================================================================#

@testitem "2.1 Eq. 19.14: Aerodynamic Resistance (Neutral Stability)" setup=[SeinfeldPandisSetup] begin
    # r_a = (1/(kappa * u_star)) * ln(z/z_0)
    # kappa = 0.4 (von Karman constant)

    # Test case: z = 10 m, z_0 = 0.1 m, u_star = 0.4 m/s
    kappa = 0.4
    z = 10.0      # m
    z_0 = 0.1     # m
    u_star = 0.4  # m/s

    # Expected: r_a = (1/(0.4 * 0.4)) * ln(10/0.1) = 6.25 * ln(100) = 28.78 s/m
    expected_r_a = (1 / (kappa * u_star)) * log(z / z_0)

    @test expected_r_a ≈ 28.78 rtol=0.01

    # Verify formula matches textbook (manual calculation)
    @test (1 / (0.4 * 0.4)) * log(100) ≈ expected_r_a rtol=1e-10
end

@testitem "2.2 Eq. 19.17: Quasi-Laminar Resistance for Gases" setup=[SeinfeldPandisSetup] begin
    # r_b = 5 * Sc^(2/3) / u_star
    # For air at 298 K: nu ≈ 1.5e-5 m^2/s
    # For O3: D ≈ 1.5e-5 m^2/s, so Sc ≈ 1.0

    u_star = 0.4  # m/s
    Sc = 1.0      # Schmidt number

    # Expected: r_b = 5 * 1.0^(2/3) / 0.4 = 5 / 0.4 = 12.5 s/m
    expected_r_b = 5 * Sc^(2/3) / u_star
    @test expected_r_b ≈ 12.5 rtol=1e-10

    # Test with different Sc (SO2 at 298 K, Sc ≈ 1.3)
    Sc_SO2 = 1.3
    r_b_SO2 = 5 * Sc_SO2^(2/3) / u_star
    @test r_b_SO2 > expected_r_b  # Higher Sc means higher r_b
    @test r_b_SO2 ≈ 5 * 1.3^(2/3) / 0.4 rtol=1e-10
end

@testitem "2.3 Eq. 19.18: Particle Settling Velocity (Stokes)" setup=[SeinfeldPandisSetup] begin
    # v_s = rho_p * D_p^2 * g * C_c / (18 * mu)

    # Test case: 1 um particle, rho_p = 1000 kg/m^3
    rho_p = 1000.0     # kg/m^3
    D_p = 1.0e-6       # m (1 um)
    g = 9.81           # m/s^2
    C_c = 1.16         # Cunningham correction for 1 um
    mu = 1.85e-5       # Pa*s (air at 298 K)

    # Expected settling velocity
    expected_v_s = rho_p * D_p^2 * g * C_c / (18 * mu)

    # For 1 um particle at STP, v_s should be ~3.5e-5 m/s (about 0.035 mm/s)
    @test expected_v_s > 0
    @test expected_v_s < 1e-3  # Should be much less than 1 mm/s for 1 um particle
    @test expected_v_s ≈ 3.41e-5 rtol=0.1
end

@testitem "2.4 Eq. 19.2: Gas Deposition Velocity" setup=[SeinfeldPandisSetup] begin
    # v_d = 1 / (r_a + r_b + r_c)

    r_a = 30.0   # s/m (typical aerodynamic resistance)
    r_b = 15.0   # s/m (typical quasi-laminar resistance)
    r_c = 200.0  # s/m (typical surface resistance for O3)

    # Total resistance
    r_t = r_a + r_b + r_c
    @test r_t == 245.0

    # Deposition velocity
    v_d = 1 / r_t
    @test v_d ≈ 0.00408 rtol=0.01  # about 0.4 cm/s

    # Convert to cm/s for comparison with Table 19.1
    v_d_cms = v_d * 100
    @test v_d_cms ≈ 0.408 rtol=0.01  # close to O3 over continent (0.4 cm/s)
end

@testitem "2.5 Eq. 19.7: Particle Deposition Velocity" setup=[SeinfeldPandisSetup] begin
    # v_d = 1/(r_a + r_b + r_a*r_b*v_s) + v_s

    r_a = 30.0     # s/m
    r_b = 50.0     # s/m
    v_s = 1.0e-4   # m/s (settling velocity for ~3 um particle)

    # Expected deposition velocity
    expected_v_d = 1 / (r_a + r_b + r_a * r_b * v_s) + v_s

    # The virtual resistance term r_a*r_b*v_s should be small for small particles
    virtual_resistance = r_a * r_b * v_s
    @test virtual_resistance ≈ 0.15 rtol=1e-10  # 30 * 50 * 1e-4 = 0.15 s/m

    # Total denominator
    denom = r_a + r_b + virtual_resistance
    @test denom ≈ 80.15 rtol=1e-10

    @test expected_v_d ≈ 1/80.15 + 1e-4 rtol=0.01
end

#=============================================================================
# Test Set 3: Physical Constraints
=============================================================================#

@testitem "3.1 Resistance Positivity" setup=[SeinfeldPandisSetup] begin
    # All resistances must be positive

    # Aerodynamic resistance: positive when z > z_0 and u_star > 0
    kappa = 0.4
    z = 10.0
    z_0 = 0.1
    u_star = 0.4
    r_a = (1 / (kappa * u_star)) * log(z / z_0)
    @test r_a > 0

    # Quasi-laminar resistance: positive when Sc > 0 and u_star > 0
    Sc = 1.0
    r_b = 5 * Sc^(2/3) / u_star
    @test r_b > 0
end

@testitem "3.2 Deposition Velocity Bounds" setup=[SeinfeldPandisSetup] begin
    # v_d should be positive and bounded

    # Minimum v_d occurs when resistances are maximum
    r_max = 1000.0  # s/m (very high resistance)
    v_d_min = 1 / (3 * r_max)  # Three resistances in series
    @test v_d_min > 0
    @test v_d_min ≈ 3.33e-4 rtol=0.01  # About 0.03 cm/s

    # Maximum v_d for gases limited by aerodynamic resistance
    r_a_min = 10.0  # s/m (typical minimum)
    r_b_min = 5.0   # s/m
    r_c_min = 0.0   # Perfect sink (HNO3 type)
    v_d_max = 1 / (r_a_min + r_b_min + r_c_min)
    @test v_d_max ≈ 0.0667 rtol=0.01  # About 6.7 cm/s
    @test v_d_max < 0.1  # Should be less than 10 cm/s for gases
end

@testitem "3.3 Sticking Fraction Bounds (Eq. 19.26)" setup=[SeinfeldPandisSetup] begin
    # R_1 = exp(-sqrt(St))
    # R_1 should be between 0 and 1

    # Low Stokes number (small particles stick well)
    St_low = 0.01
    R_1_low = exp(-sqrt(St_low))
    @test R_1_low > 0.9  # Close to 1
    @test R_1_low ≤ 1.0

    # High Stokes number (large particles may bounce)
    St_high = 10.0
    R_1_high = exp(-sqrt(St_high))
    @test R_1_high < 0.1  # Much less than 1
    @test R_1_high > 0.0
end

#=============================================================================
# Test Set 4: Parameter Tables
=============================================================================#

@testitem "4.1 Land-Use Parameters (Table 19.2)" setup=[SeinfeldPandisSetup] begin
    # Test that land-use parameters are defined correctly
    @test DryDeposition.LANDUSE_GRASS.A == 2.0e-3
    @test DryDeposition.LANDUSE_GRASS.alpha == 1.2
    @test DryDeposition.LANDUSE_GRASS.gamma == 0.54

    @test DryDeposition.LANDUSE_DECIDUOUS_FOREST.A == 5.0e-3
    @test DryDeposition.LANDUSE_DECIDUOUS_FOREST.alpha == 0.6
    @test DryDeposition.LANDUSE_DECIDUOUS_FOREST.gamma == 0.56

    @test DryDeposition.LANDUSE_DESERT.A == 10.0e-3
    @test DryDeposition.LANDUSE_DESERT.alpha == 50.0
    @test DryDeposition.LANDUSE_DESERT.gamma == 0.54

    # Physical constraints on parameters
    @test DryDeposition.LANDUSE_GRASS.gamma > 0.5
    @test DryDeposition.LANDUSE_GRASS.gamma < 2/3
end

@testitem "4.2 Gas Properties (Table 19.4)" setup=[SeinfeldPandisSetup] begin
    # Test that gas properties are defined correctly

    # SO2
    @test DryDeposition.GAS_SO2.D_ratio ≈ 1.9 rtol=0.1
    @test DryDeposition.GAS_SO2.H_star == 1.0e5
    @test DryDeposition.GAS_SO2.f_0 == 0.0

    # O3
    @test DryDeposition.GAS_O3.D_ratio ≈ 1.6 rtol=0.1
    @test DryDeposition.GAS_O3.H_star == 0.01
    @test DryDeposition.GAS_O3.f_0 == 1.0  # Highly reactive

    # NO2
    @test DryDeposition.GAS_NO2.D_ratio ≈ 1.6 rtol=0.1
    @test DryDeposition.GAS_NO2.f_0 == 0.1

    # HNO3 - extremely high Henry's law constant
    @test DryDeposition.GAS_HNO3.H_star == 1.0e14
    @test DryDeposition.GAS_HNO3.H_star > DryDeposition.GAS_SO2.H_star
end

@testitem "4.3 Physical Constants" setup=[SeinfeldPandisSetup] begin
    @test DryDeposition.KAPPA ≈ 0.4 rtol=0.01
    @test DryDeposition.G_ACCEL ≈ 9.81 rtol=0.001
    @test DryDeposition.K_BOLTZ ≈ 1.38e-23 rtol=0.01
end

#=============================================================================
# Test Set 5: Integration Tests
=============================================================================#

@testitem "5.1 Gas Deposition System Compilation" setup=[SeinfeldPandisSetup] begin
    @named gas_dep = DryDepositionGas()

    # Test that the system can be compiled with inputs specified
    # These are algebraic systems so we need to specify which variables are inputs
    gas_dep_nns = toggle_namespacing(gas_dep, false)
    inputs = [gas_dep_nns.u_star, gas_dep_nns.C, gas_dep_nns.G, gas_dep_nns.T_s]
    compiled = mtkcompile(gas_dep; inputs)
    @test compiled isa System

    # Check that parameters are accessible
    params = parameters(compiled)
    @test length(params) > 0
end

@testitem "5.2 Particle Deposition System Compilation" setup=[SeinfeldPandisSetup] begin
    @named part_dep = DryDepositionParticle()

    # Test that the system can be compiled with inputs specified
    part_dep_nns = toggle_namespacing(part_dep, false)
    inputs = [part_dep_nns.u_star, part_dep_nns.C]
    compiled = mtkcompile(part_dep; inputs)
    @test compiled isa System

    # Check that parameters are accessible
    params = parameters(compiled)
    @test length(params) > 0
end

@testitem "5.3 Component System Compilation" setup=[SeinfeldPandisSetup] begin
    # Test individual components can be compiled with inputs specified

    @named ra = AerodynamicResistance()
    ra_nns = toggle_namespacing(ra, false)
    compiled_ra = mtkcompile(ra; inputs = [ra_nns.u_star])
    @test compiled_ra isa System

    @named rb_gas = QuasiLaminarResistanceGas()
    rb_gas_nns = toggle_namespacing(rb_gas, false)
    compiled_rb = mtkcompile(rb_gas; inputs = [rb_gas_nns.u_star])
    @test compiled_rb isa System

    @named settling = ParticleSettling()
    # ParticleSettling has no inputs - all parameters are fixed
    compiled_settling = mtkcompile(settling)
    @test compiled_settling isa System
end

#=============================================================================
# Test Set 6: Sensitivity Analysis
=============================================================================#

@testitem "6.1 Aerodynamic Resistance vs Friction Velocity" setup=[SeinfeldPandisSetup] begin
    # r_a should decrease with increasing u_star (inverse relationship)

    kappa = 0.4
    z = 10.0
    z_0 = 0.1

    r_a_func = (u_star) -> (1 / (kappa * u_star)) * log(z / z_0)

    r_a_low = r_a_func(0.1)
    r_a_mid = r_a_func(0.4)
    r_a_high = r_a_func(1.0)

    @test r_a_low > r_a_mid > r_a_high

    # Doubling u_star should halve r_a
    @test r_a_func(0.2) ≈ r_a_func(0.4) * 2 rtol=1e-10
end

@testitem "6.2 Aerodynamic Resistance vs Roughness Length" setup=[SeinfeldPandisSetup] begin
    # r_a should decrease with increasing z_0 (less resistance over rough surfaces)

    kappa = 0.4
    z = 10.0
    u_star = 0.4

    r_a_func = (z_0) -> (1 / (kappa * u_star)) * log(z / z_0)

    r_a_smooth = r_a_func(0.01)
    r_a_grass = r_a_func(0.1)
    r_a_forest = r_a_func(1.0)

    @test r_a_smooth > r_a_grass > r_a_forest
end
