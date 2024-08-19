using AtmosphericDeposition
using Test, Unitful, ModelingToolkit, GasChem, Dates, EarthSciMLBase

@testset "Connector" begin
    ModelingToolkit.check_units(eqs...) = nothing
    start = Dates.datetime2unix(Dates.DateTime(2016, 5, 1))
    @parameters t
    composed_ode = couple(SuperFast(t), FastJX(t), DrydepositionG(t), Wetdeposition(t))
    tspan = (start, start+3600*24*3)
    sys = structural_simplify(get_mtk(composed_ode))
    @test length(states(sys)) ≈ 18

    eqs = string(equations(sys))
    wanteqs = ["Differential(t)(superfast₊O3(t)) ~ superfast₊DrydepositionG_ddt_O3ˍt(t) + superfast₊WetDeposition_ddt_O3ˍt(t)"]
    @test contains(string(eqs), wanteqs[1])
end
