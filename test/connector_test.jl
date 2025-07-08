using AtmosphericDeposition
using Test, DynamicQuantities, ModelingToolkit, Dates, EarthSciMLBase
using OrdinaryDiffEqRosenbrock
using EarthSciData, GasChem, Aerosol
using ModelingToolkit: t

domain = DomainInfo(
    DateTime(2016, 2, 1),
    DateTime(2016, 2, 2);
    latrange = deg2rad(-85.0f0):deg2rad(2):deg2rad(85.0f0),
    lonrange = deg2rad(-180.0f0):deg2rad(2.5):deg2rad(175.0f0),
    levrange = 1:10,
)

@testset "GasChemExt" begin
    start = Dates.datetime2unix(Dates.DateTime(2016, 5, 1))
    composed_ode = couple(SuperFast(), FastJX(0.0), DryDepositionGas(), WetDeposition())
    sys = convert(ODESystem, composed_ode)

    eqs = string(equations(sys))
    @test contains(string(eqs), "SuperFast₊DryDepositionGas_k_O3(t)")
    @test contains(string(eqs), "SuperFast₊WetDeposition_k_othergas(t)")
    @test contains(
        string(observed(sys)),
        "SuperFast₊DryDepositionGas_k_O3(t) ~ -DryDepositionGas₊k_O3(t)*SuperFast₊O3(t)"
    )
end

@testset "AerosolExt" begin
    model = couple(
        GEOSFP("4x5", domain),
        WetDeposition(),
        ElementalCarbon(),
        DryDepositionAerosol()
    )
    sys = convert(ODESystem, model)

    eqs = equations(sys)
    @test contains(string(eqs), "ElementalCarbon₊DryDepositionAerosol_k")
    @test contains(string(eqs), "ElementalCarbon₊WetDeposition_k_particle")
end

@testset "EarthSciDataExt" begin
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

    sys = convert(ODESystem, model)

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

@testset "GEOSChemGasPhase_coupling" begin
    model = couple(
        GEOSChemGasPhase(),
        WetDeposition(),
        DryDepositionGas(),
        )
        model_sys = convert(ODESystem, model)
        vals = ModelingToolkit.get_defaults(model_sys)
        for key in setdiff(unknowns(model_sys),keys(vals))
            vals[key] = 0 # Set variables with no default to zero.
        end
        prob = ODEProblem(model_sys, u0 = vals,  p = vals)

        sol = solve(prob, Rosenbrock23(), tspan=(0,3600), saveat = 60.)
        @test sol.retcode == ReturnCode.Success

end