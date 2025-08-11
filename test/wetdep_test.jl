@testsnippet WetDepSetup begin
using AtmosphericDeposition
using AtmosphericDeposition: _WetDeposition, get_lev_depth, wd_defaults
using Test, DynamicQuantities, ModelingToolkit

@parameters cloudFrac = 0.5
@parameters qrain = 0.5
@parameters ρ_air = 1.204 [unit = u"kg*m^-3"]
@parameters Δz = 200 [unit = u"m"]
@parameters lev

@constants Δz_unit = 1 [unit = u"m", description = "unit of depth"]
end

@testitem "unit" setup=[WetDepSetup] begin
    @test substitute(
        _WetDeposition(cloudFrac, qrain, ρ_air, Δz)[1],
        Dict(cloudFrac => 0.5, qrain => 0.5, ρ_air => 1.204, Δz => 200, wd_defaults...)
    ) ≈ 0.313047525
    @test ModelingToolkit.get_unit(_WetDeposition(cloudFrac, qrain, ρ_air, Δz)[1]) ==
          u"s^-1"
    @test ModelingToolkit.get_unit(_WetDeposition(cloudFrac, qrain, ρ_air, Δz)[2]) ==
          u"s^-1"
    @test ModelingToolkit.get_unit(_WetDeposition(cloudFrac, qrain, ρ_air, Δz)[3]) ==
          u"s^-1"
    @test ModelingToolkit.get_unit(
        _WetDeposition(cloudFrac, qrain, ρ_air, get_lev_depth(lev) * Δz_unit)[3],
    ) == u"s^-1"
end

@testitem "WetDeposition" setup=[WetDepSetup] begin
    @test substitute(
        _WetDeposition(cloudFrac, qrain, ρ_air, Δz)[1],
        Dict(cloudFrac => 0.5, qrain => -1e5, ρ_air => 1.204, Δz => 200, wd_defaults...)
    ) ≈ 0.0
    @test substitute(get_lev_depth(lev), Dict(lev => 3)) ≈ 127.81793001768432
end
