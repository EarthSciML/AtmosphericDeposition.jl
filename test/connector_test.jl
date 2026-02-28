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

@testitem "GasChemExt SuperFast DryDeposition" begin
    using AtmosphericDeposition, GasChem, EarthSciMLBase, ModelingToolkit
    using Test

    model = couple(SuperFast(), DryDepositionGas())
    sys = convert(System, model)
    eqs = string(equations(sys))

    # Verify that GasChem species are coupled to dry deposition rate constants
    @test contains(eqs, "SuperFast₊DryDepositionGas_k_HNO3")
    @test contains(eqs, "SuperFast₊DryDepositionGas_k_NO2")
    @test contains(eqs, "SuperFast₊DryDepositionGas_k_O3")
    @test contains(eqs, "SuperFast₊DryDepositionGas_k_H2O2")
    @test contains(eqs, "SuperFast₊DryDepositionGas_k_HCHO")
end

@testitem "GasChemExt SuperFast WetDeposition" begin
    using AtmosphericDeposition, GasChem, EarthSciMLBase, ModelingToolkit
    using Test

    model = couple(SuperFast(), WetDeposition())
    sys = convert(System, model)
    eqs = string(equations(sys))

    # Verify that GasChem species are coupled to wet deposition rate constants
    @test contains(eqs, "SuperFast₊WetDeposition_k_othergas")
end

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
