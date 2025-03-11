using AtmosphericDeposition
using Test, DynamicQuantities, ModelingToolkit, GasChem, Dates, EarthSciMLBase, EarthSciData
using ModelingToolkit:t

domain = DomainInfo(DateTime(2022, 1, 1), DateTime(2022, 1, 3);
    latrange=deg2rad(-85.0f0):deg2rad(2):deg2rad(85.0f0),
    lonrange=deg2rad(-180.0f0):deg2rad(2.5):deg2rad(175.0f0),
    levrange=1:10, dtype=Float64)

@testset "GasChemExt" begin
    start = Dates.datetime2unix(Dates.DateTime(2016, 5, 1))
    composed_ode = couple(SuperFast(), FastJX(), DrydepositionG(), Wetdeposition())
    sys = convert(ODESystem, composed_ode)
    print(unknowns(sys))
    @test length(unknowns(sys)) ≈ 12

    eqs = string(equations(sys))
    wanteqs = ["Differential(t)(SuperFast₊O3(t)) ~ SuperFast₊DrydepositionG_ddt_O3ˍt(t) + SuperFast₊Wetdeposition_ddt_O3ˍt(t)"]
    @test contains(string(eqs), wanteqs[1])
end

@testset "EarthSciDataExt" begin
    @parameters lat = deg2rad(40.0f0) [unit=u"rad"]
    @parameters lon = deg2rad(-97.0f0) [unit=u"rad"]
    @parameters lev = 1

    geosfp = GEOSFP("4x5", domain)

    model = couple(SuperFast(), FastJX(), geosfp, Wetdeposition(), DrydepositionG())

    sys = convert(ODESystem, model)
    @test length(unknowns(sys)) ≈ 12

    eqs = string(observed(sys))
    wanteq = "DrydepositionG₊G(t) ~ GEOSFP₊A1₊SWGDN(t)"
    @test contains(eqs, wanteq)
    wanteq = "Wetdeposition₊cloudFrac(t) ~ GEOSFP₊A3cld₊CLOUD(t)"
    @test contains(eqs, wanteq)
    wanted = "Wetdeposition₊ρA(t) ~ GEOSFP₊P/(GEOSFP₊I3₊T*R)*kgperg*MW_air"
    @test contains(eqs, wanteq)
end
