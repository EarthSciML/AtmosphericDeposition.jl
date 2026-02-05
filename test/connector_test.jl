@testsnippet ConnectorSetup begin
    using AtmosphericDeposition
    using Test, ModelingToolkit, Dates, EarthSciMLBase
    using OrdinaryDiffEqRosenbrock
    using EarthSciData, GasChem, Aerosol

    domain = DomainInfo(
        DateTime(2016, 2, 1),
        DateTime(2016, 2, 2);
        latrange = deg2rad(-85.0f0):deg2rad(2):deg2rad(85.0f0),
        lonrange = deg2rad(-180.0f0):deg2rad(2.5):deg2rad(175.0f0),
        levrange = 1:10
    )
end

@testitem "GasChemExt" begin
    using Dates, ModelingToolkit
    using GasChem, EarthSciMLBase, AtmosphericDeposition
    start = Dates.datetime2unix(Dates.DateTime(2016, 5, 1))
    composed_ode = couple(SuperFast(), FastJX(0.0), DryDepositionGas(), WetDeposition())
    sys = convert(System, composed_ode)

    eqs = string(equations(sys))
    @test contains(string(eqs), "SuperFast₊DryDepositionGas_k_O3(t)")
    @test contains(string(eqs), "SuperFast₊WetDeposition_k_othergas(t)")
    @test contains(
        string(observed(sys)),
        "SuperFast₊DryDepositionGas_k_O3(t) ~ -DryDepositionGas₊k_O3(t)*SuperFast₊O3(t)"
    )
end

@testitem "AerosolExt" setup=[ConnectorSetup] begin
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

@testitem "EarthSciDataExt" setup=[ConnectorSetup] begin
    model = couple(
        SuperFast(),
        FastJX(0.0),
        GEOSFP("4x5", domain),
        NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain),
        WetDeposition(),
        DryDepositionGas(),
        ElementalCarbon(),
        DryDepositionAerosol()
    )

    sys = convert(System, model)

    @test length(unknowns(sys)) ≈ 13

    eqs = string(observed(sys))
    wanteq = "DryDepositionGas₊G(t) ~ GEOSFP₊A1₊SWGDN(t)"
    @test contains(eqs, wanteq)
    wanteq = "GEOSFP₊A1₊USTAR(t)"
    @test contains(eqs, wanteq)
    wanteq = "WetDeposition₊cloudFrac(t) ~ GEOSFP₊A3cld₊CLOUD(t)"
    @test contains(eqs, wanteq)
    wanted = "WetDeposition₊ρA(t) ~ GEOSFP₊P/(GEOSFP₊I3₊T*R)*kgperg*MW_air"
    @test contains(eqs, wanteq)
    @test contains(eqs, "NEI2016MonthlyEmis")
end

@testitem "GEOSChemGasPhase_coupling" setup=[ConnectorSetup] begin
    model = couple(
        GEOSChemGasPhase(),
        WetDeposition(),
        DryDepositionGas()
    )
    model_sys = convert(System, model)
    prob = ODEProblem(model_sys, [], get_tspan(domain))

    sol = solve(prob, Rosenbrock23(), tspan = (0, 3600), saveat = 60.0)
    @test sol.retcode == ReturnCode.Success
end
