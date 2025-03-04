using AtmosphericDeposition
using Test, DynamicQuantities, ModelingToolkit

@parameters cloudFrac = 0.5
@parameters qrain = 0.5
@parameters ρ_air = 1.204 [unit = u"kg*m^-3"]
@parameters Δz = 200 [unit = u"m"]

@testset "unit" begin
    @test substitute(WetDeposition(cloudFrac, qrain, ρ_air, Δz)[1], Dict(cloudFrac => 0.5, qrain => 0.5, ρ_air => 1.204, Δz => 200, wd_defaults...)) ≈ 0.313047525
    @test ModelingToolkit.get_unit(WetDeposition(cloudFrac, qrain, ρ_air, Δz)[1]) == u"s^-1"
    @test ModelingToolkit.get_unit(WetDeposition(cloudFrac, qrain, ρ_air, Δz)[2]) == u"s^-1"
    @test ModelingToolkit.get_unit(WetDeposition(cloudFrac, qrain, ρ_air, Δz)[3]) == u"s^-1"
end

@testset "WetDeposition" begin
    @test substitute(WetDeposition(cloudFrac, qrain, ρ_air, Δz)[1], Dict(cloudFrac => 0.5, qrain => -1e5, ρ_air => 1.204, Δz => 200, wd_defaults...)) ≈ 0.0
end