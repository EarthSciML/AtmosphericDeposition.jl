@testsnippet WetDepSetup begin
    using AtmosphericDeposition
    using AtmosphericDeposition: _WetDeposition, get_lev_depth, wd_defaults
    using Test, DynamicQuantities, ModelingToolkit
    using Symbolics
    using Symbolics: value

    function to_float(x)
        v = value(x)
        v isa Number && return Float64(v)
        return Float64(eval(Symbolics.toexpr(Symbolics.unwrap(x))))
    end

    @parameters cloudFrac = 0.5
    @parameters qrain = 0.5
    @parameters ρ_air = 1.204 [unit = u"kg*m^-3"]
    @parameters Δz = 200 [unit = u"m"]
    @parameters lev

    @constants Δz_unit = 1 [unit = u"m", description = "unit of depth"]
end

@testitem "unit" setup = [WetDepSetup] begin
    @test to_float(
        substitute(
            _WetDeposition(cloudFrac, qrain, ρ_air, Δz)[1],
            Dict(cloudFrac => 0.5, qrain => 0.5, ρ_air => 1.204, Δz => 200, wd_defaults...)
        )
    ) ≈ 1.81804  # = qrain*ρ_air*(A_wd*E + cf*W_in,particle*Vdr/(ρwater*h_s)) with the
    # EMEP 1e6 scavenging-ratio scale + fixed h_s=1000 m (was 0.313047525 when the
    # in-cloud term used the bare Table-9.1 entry and the local Δz)
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

@testitem "WetDeposition" setup = [WetDepSetup] begin
    @test to_float(
        substitute(
            _WetDeposition(cloudFrac, qrain, ρ_air, Δz)[1],
            Dict(cloudFrac => 0.5, qrain => -1.0e5, ρ_air => 1.204, Δz => 200, wd_defaults...)
        )
    ) ≈ 0.0
    @test to_float(substitute(get_lev_depth(lev), Dict(lev => 3))) ≈ 127.81793001768432
end

@testitem "physical magnitude at 1 mm/hr" setup = [WetDepSetup] begin
    # 1 mm/hr => P = 2.78e-4 kg m^-2 s^-1, and qrain = P / (Vdr * ρ_air).
    P = 2.78e-4
    ρ = 1.2
    q = P / (5.0 * ρ)
    subs = Dict(cloudFrac => 0.5, qrain => q, ρ_air => ρ, Δz => 124.0, wd_defaults...)
    wd = _WetDeposition(cloudFrac, qrain, ρ_air, Δz)

    k_particle = to_float(substitute(wd[1], subs))
    k_SO2 = to_float(substitute(wd[2], subs))
    k_othergas = to_float(substitute(wd[3], subs))

    # Moderate rain scavenges on timescales of minutes to hours.
    @test 1.0e-5 < k_particle < 1.0e-1
    @test 1.0e-5 < k_SO2 < 1.0e-1
    @test 1.0e-5 < k_othergas < 1.0e-1
    @test k_SO2 < k_othergas
end

@testitem "in-cloud term uses the in-cloud scavenging ratio" setup = [WetDepSetup] begin
    # Clear sky -> fully overcast adds exactly one in-cloud term, so
    # k(cloudFrac=1) / k(cloudFrac=0) = (W_sub + W_in) / W_sub.
    base = Dict(qrain => 1.0e-4, ρ_air => 1.2, Δz => 124.0, wd_defaults...)
    wd = _WetDeposition(cloudFrac, qrain, ρ_air, Δz)
    clear(k) = to_float(substitute(wd[k], merge(base, Dict(cloudFrac => 0.0))))
    overcast(k) = to_float(substitute(wd[k], merge(base, Dict(cloudFrac => 1.0))))

    @test overcast(2) / clear(2) ≈ (0.15 + 0.3) / 0.15   # SO2:       3.0
    @test overcast(3) / clear(3) ≈ (0.5 + 1.4) / 0.5     # other gas: 3.8
end
