using AtmosphericDeposition
using Test, Unitful, ModelingToolkit, GasChem, Dates, EarthSciMLBase

@testset "Connector" begin
    ModelingToolkit.check_units(eqs...) = nothing
    start = Dates.datetime2unix(Dates.DateTime(2016, 5, 1))
    @parameters t
    composed_ode = SuperFast(t) + FastJX(t) + DrydepositionG(t) + Wetdeposition(t)
    tspan = (start, start+3600*24*3)
    sys = structural_simplify(get_mtk(composed_ode))
    @test length(states(sys)) ≈ 18
end
