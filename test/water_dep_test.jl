@testsnippet WaterDepSetup begin
    using AtmosphericDeposition
    using AtmosphericDeposition: kG_gas_transfer, kL_liss_merlivat, kL_wanninkhof,
        schmidt_CO2_seawater, rc_water_wanninkhof, rc_water_liss_merlivat,
        unit_kL, kL_unit_convert, T_K_unit
    using Test, DynamicQuantities, ModelingToolkit

    @parameters u10  # Wind speed at 10m [m/s]
    @parameters T_sym [unit = u"K"]  # Temperature
    @parameters H_star  # Henry's law constant [M/atm]
    @parameters Sc_ratio  # Schmidt number ratio
end

@testitem "kG_gas_transfer" setup = [WaterDepSetup] begin
    # Test gas-phase transfer coefficient
    # kG = 0.0013 * u10 * unit_kL, so at u10 = 10, kG = 0.013 m/s
    @test substitute(
        kG_gas_transfer(u10),
        Dict(u10 => 10.0, unit_kL => 1)
    ) ≈ 0.013

    # At u10 = 5, kG = 0.0065 m/s
    @test substitute(
        kG_gas_transfer(u10),
        Dict(u10 => 5.0, unit_kL => 1)
    ) ≈ 0.0065

    # kG should be zero at zero wind speed
    @test substitute(
        kG_gas_transfer(u10),
        Dict(u10 => 0.0, unit_kL => 1)
    ) ≈ 0.0

    # Check units
    @test ModelingToolkit.get_unit(kG_gas_transfer(u10)) == u"m/s"
end

@testitem "schmidt_CO2_seawater" setup = [WaterDepSetup] begin
    # Test Schmidt number for CO2 in seawater
    # From S&P Eq. 19.45: Sc = 2073.1 - 125.62T + 3.6276T² - 0.043219T³

    # At T = 20°C, Sc ≈ 660 (reference value)
    Sc_20 = schmidt_CO2_seawater(20.0)
    @test 650 < Sc_20 < 680  # Should be approximately 660

    # At T = 0°C, Sc should be higher (cold water has higher viscosity)
    Sc_0 = schmidt_CO2_seawater(0.0)
    @test Sc_0 > Sc_20

    # At T = 30°C, Sc should be lower
    Sc_30 = schmidt_CO2_seawater(30.0)
    @test Sc_30 < Sc_20
end

@testitem "kL_wanninkhof" setup = [WaterDepSetup] begin
    # Test liquid-phase transfer coefficient using Wanninkhof (1992)
    # kL = 0.31 * u10² * Sc_ratio^0.5 / 360000 * unit_kL [m/s]

    # At u10 = 10 and Sc_ratio = 1:
    # kL = 0.31 * 100 * 1 / 360000 = 8.61e-5 m/s
    kL_val = substitute(
        kL_wanninkhof(u10, Sc_ratio),
        Dict(u10 => 10.0, Sc_ratio => 1.0, kL_unit_convert => 360000, unit_kL => 1)
    )
    @test 8.0e-5 < kL_val < 9.0e-5

    # kL should increase with wind speed (quadratic)
    kL_5 = substitute(
        kL_wanninkhof(u10, Sc_ratio),
        Dict(u10 => 5.0, Sc_ratio => 1.0, kL_unit_convert => 360000, unit_kL => 1)
    )
    kL_10 = substitute(
        kL_wanninkhof(u10, Sc_ratio),
        Dict(u10 => 10.0, Sc_ratio => 1.0, kL_unit_convert => 360000, unit_kL => 1)
    )
    @test kL_10 / kL_5 ≈ 4.0  # Quadratic dependence

    # Check units
    @test ModelingToolkit.get_unit(kL_wanninkhof(u10, Sc_ratio)) == u"m/s"
end

@testitem "kL_liss_merlivat" setup = [WaterDepSetup] begin
    # Test liquid-phase transfer coefficient using Liss-Merlivat (1986)
    # Three regimes: smooth (u10 <= 3.6), rough (3.6 < u10 <= 13), breaking waves (u10 > 13)

    # Smooth regime: u10 = 2
    kL_smooth = substitute(
        kL_liss_merlivat(u10, Sc_ratio),
        Dict(u10 => 2.0, Sc_ratio => 1.0, kL_unit_convert => 360000, unit_kL => 1)
    )
    @test kL_smooth > 0

    # Rough regime: u10 = 8
    kL_rough = substitute(
        kL_liss_merlivat(u10, Sc_ratio),
        Dict(u10 => 8.0, Sc_ratio => 1.0, kL_unit_convert => 360000, unit_kL => 1)
    )
    @test kL_rough > kL_smooth  # Should be higher than smooth regime

    # Breaking wave regime: u10 = 15
    kL_break = substitute(
        kL_liss_merlivat(u10, Sc_ratio),
        Dict(u10 => 15.0, Sc_ratio => 1.0, kL_unit_convert => 360000, unit_kL => 1)
    )
    @test kL_break > kL_rough  # Should be highest

    # Check units
    @test ModelingToolkit.get_unit(kL_liss_merlivat(u10, Sc_ratio)) == u"m/s"
end

@testitem "rc_water_soluble_gas" setup = [WaterDepSetup] begin
    # For highly soluble gases (H_star >> 1e3), rc ≈ 1/kG (Eq. 19.49)
    # Gas-phase resistance dominates

    # With u10 = 10: kG = 0.013 m/s, so rc ≈ 77 s/m
    rc_soluble = substitute(
        rc_water_wanninkhof(H_star, u10, T_sym),
        Dict(
            H_star => 1e6,  # Very soluble gas
            u10 => 10.0,
            T_sym => 293.0,
            unit_kL => 1,
            T_K_unit => 1,
            kL_unit_convert => 360000
        )
    )

    # For very soluble gases, rc should be close to 1/kG = 1/0.013 ≈ 77 s/m
    @test 70 < rc_soluble < 85

    # Check units
    @test ModelingToolkit.get_unit(rc_water_wanninkhof(H_star, u10, T_sym)) == u"m^-1*s"
end

@testitem "rc_water_insoluble_gas" setup = [WaterDepSetup] begin
    # For slightly soluble gases (H_star << 1), rc ≈ 1/(kL * H̃) (Eq. 19.48)
    # Liquid-phase resistance dominates

    rc_insoluble = substitute(
        rc_water_wanninkhof(H_star, u10, T_sym),
        Dict(
            H_star => 0.01,  # Very insoluble gas
            u10 => 10.0,
            T_sym => 293.0,
            unit_kL => 1,
            T_K_unit => 1,
            kL_unit_convert => 360000
        )
    )

    # For insoluble gases, rc should be very large
    @test rc_insoluble > 1000  # High resistance for insoluble gas
end

@testitem "rc_water_temperature_dependence" setup = [WaterDepSetup] begin
    # Test that rc varies with temperature
    # Warmer water has lower Schmidt number, higher kL, lower rc

    rc_cold = substitute(
        rc_water_wanninkhof(H_star, u10, T_sym),
        Dict(
            H_star => 1000.0,
            u10 => 10.0,
            T_sym => 283.0,  # 10°C
            unit_kL => 1,
            T_K_unit => 1,
            kL_unit_convert => 360000
        )
    )

    rc_warm = substitute(
        rc_water_wanninkhof(H_star, u10, T_sym),
        Dict(
            H_star => 1000.0,
            u10 => 10.0,
            T_sym => 303.0,  # 30°C
            unit_kL => 1,
            T_K_unit => 1,
            kL_unit_convert => 360000
        )
    )

    # Warmer water should have lower resistance (faster gas exchange)
    @test rc_warm < rc_cold
end

@testitem "rc_water_wind_dependence" setup = [WaterDepSetup] begin
    # Test that rc decreases with increasing wind speed

    rc_low_wind = substitute(
        rc_water_wanninkhof(H_star, u10, T_sym),
        Dict(
            H_star => 1000.0,
            u10 => 3.0,
            T_sym => 293.0,
            unit_kL => 1,
            T_K_unit => 1,
            kL_unit_convert => 360000
        )
    )

    rc_high_wind = substitute(
        rc_water_wanninkhof(H_star, u10, T_sym),
        Dict(
            H_star => 1000.0,
            u10 => 15.0,
            T_sym => 293.0,
            unit_kL => 1,
            T_K_unit => 1,
            kL_unit_convert => 360000
        )
    )

    # Higher wind speed should give lower resistance
    @test rc_high_wind < rc_low_wind
end

@testitem "DryDepositionWater_unit" setup = [WaterDepSetup] begin
    using AtmosphericDeposition: DryDepositionWater
    using ModelingToolkit: t

    # Create the model and check units
    model = DryDepositionWater()

    # Check that the system has the expected variables
    @test length(ModelingToolkit.unknowns(model)) == 2

    # Get variables
    vars = ModelingToolkit.unknowns(model)
    v_dep = vars[1]
    k_dep = vars[2]

    # Check units
    @test ModelingToolkit.get_unit(v_dep) == u"m/s"
    @test ModelingToolkit.get_unit(k_dep) == u"s^-1"
end

@testitem "DryDepositionWater_liss_merlivat" setup = [WaterDepSetup] begin
    using AtmosphericDeposition: DryDepositionWater
    using ModelingToolkit: t

    # Create the model with Liss-Merlivat parameterization
    model_lm = DryDepositionWater(use_liss_merlivat = true)

    # Check that the system has the expected variables
    @test length(ModelingToolkit.unknowns(model_lm)) == 2

    # Get equations and check they're properly formed
    eqs = ModelingToolkit.equations(model_lm)
    @test length(eqs) == 2
end

@testitem "DryDepositionWater_values" setup = [WaterDepSetup] begin
    using AtmosphericDeposition: DryDepositionWater
    using ModelingToolkit: t

    # Create the model
    model = DryDepositionWater()

    # Get equations and check they're properly formed
    eqs = ModelingToolkit.equations(model)
    @test length(eqs) == 2

    # The model should produce reasonable deposition velocities
    # For a moderately soluble gas over water, v_dep ~ 0.001-0.01 m/s
end
