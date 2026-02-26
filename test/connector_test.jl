@testsnippet ConnectorSetup begin
    using AtmosphericDeposition
    using Test, ModelingToolkit, Dates, EarthSciMLBase
    using OrdinaryDiffEqRosenbrock
    using EarthSciData, Aerosol

    domain = DomainInfo(
        DateTime(2016, 2, 1),
        DateTime(2016, 2, 2);
        latrange = deg2rad(-85.0f0):deg2rad(2):deg2rad(85.0f0),
        lonrange = deg2rad(-180.0f0):deg2rad(2.5):deg2rad(175.0f0),
        levrange = 1:10
    )
end

# GasChem-dependent tests are temporarily disabled because Catalyst.jl
# (a GasChem dependency) does not yet support ModelingToolkit v11.
# These tests should be re-enabled once GasChem v0.12 is released.

@testitem "AerosolExt" setup = [ConnectorSetup] begin
    model = couple(
        GEOSFP("4x5", domain),
        WetDeposition(),
        ElementalCarbon(),
        DryDepositionAerosol()
    )
    sys = convert(System, model)

    eqs = equations(sys)
    @test contains(string(eqs), "ElementalCarbon₊DryDepositionAerosol_k")
    @test contains(string(eqs), "ElementalCarbon₊WetDeposition_k_particle")
end

@testitem "EarthSciDataExt" setup = [ConnectorSetup] begin
    model = couple(
        GEOSFP("4x5", domain),
        WetDeposition(),
        DryDepositionAerosol(),
        ElementalCarbon(),
    )

    sys = convert(System, model)

    eqs = string(observed(sys))
    wanteq = "GEOSFP₊A1₊USTAR(t)"
    @test contains(eqs, wanteq)
    wanteq = "WetDeposition₊cloudFrac(t) ~ GEOSFP₊A3cld₊CLOUD(t)"
    @test contains(eqs, wanteq)
end
