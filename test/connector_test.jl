using AtmosphericDeposition
using Test, Unitful, ModelingToolkit, GasChem, Dates, EarthSciMLBase, EarthSciData

@testset "GasChemExt" begin
    start = Dates.datetime2unix(Dates.DateTime(2016, 5, 1))
    @parameters t [unit = u"s"]
    composed_ode = couple(SuperFast(t), FastJX(t), DrydepositionG(t), Wetdeposition(t))
    tspan = (start, start+3600*24*3)
    sys = structural_simplify(get_mtk(composed_ode))
    @test length(states(sys)) ≈ 18

    eqs = string(equations(sys))
    wanteqs = ["Differential(t)(SuperFast₊O3(t)) ~ SuperFast₊DrydepositionG_ddt_O3ˍt(t) + SuperFast₊Wetdeposition_ddt_O3ˍt(t)"]
    @test contains(string(eqs), wanteqs[1])
end

@testset "EarthSciDataExt" begin
    @parameters t [unit = u"s"]

    @parameters lat = 40
    @parameters lon = -97
    @parameters lev = 1
    geosfp = GEOSFP("4x5", t)

    model = couple(SuperFast(t), FastJX(t), geosfp, Wetdeposition(t), DrydepositionG(t))

    sys = structural_simplify(get_mtk(model))
    @test length(states(sys)) ≈ 18

    eqs = string(equations(get_mtk(model)))
    wanteq = "DrydepositionG₊G(t) ~ GEOSFP₊A1₊SWGDN(t)"
    @test contains(eqs, wanteq)
    wanteq = "Wetdeposition₊cloudFrac(t) ~ GEOSFP₊A3cld₊CLOUD(t)"
    @test contains(eqs, wanteq)
    wanted = "Wetdeposition₊ρA(t) ~ GEOSFP₊P * PaPerhPa/(GEOSFP₊I3₊T*R)*kgperg*MW_air"
    @test contains(eqs, wanteq)
end