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
    wanteq = "WetDeposition₊cloudFrac(t) ~ GEOSFP₊A3cld₊CLOUD_itp(GEOSFP₊t_ref + t, GEOSFP₊lon, GEOSFP₊lat, GEOSFP₊lev)"
    @test contains(eqs, wanteq)
end


@testitem "dry-dep reference height is the constant floor (no Z_agl in RHS)" begin
    using EarthSciMLBase, EarthSciData, Dates
    using EarthSciMLBase: DomainInfo
    t0 = DateTime(2016, 3, 10)
    dom = DomainInfo(t0, t0 + Hour(3) + Hour(36);
        lonrange = deg2rad(-88):deg2rad(0.625):deg2rad(-86),
        latrange = deg2rad(32):deg2rad(0.5):deg2rad(33.5),
        levrange = 1:2, u_proto = zeros(Float64, 1, 1, 1, 1))
    gfp = EarthSciData.GEOSFP("0.25x0.3125", dom)
    d = DryDepositionGas()
    cs = EarthSciMLBase.couple2(EarthSciMLBase.get_coupletype(d)(d),
                                EarthSciMLBase.get_coupletype(gfp)(gfp))
    zeq = [e for e in cs.eqs if endswith(string(e.lhs), "₊z(t)")]
    @test length(zeq) == 1
    rhs = string(only(zeq).rhs)
    # The reference height must stay a plain constant: referencing the
    # height-above-ground field would be 0 at lev=1 (log(0) = NaN) and would inline
    # the multi-layer interpolation cascade into the deposition RHS.
    @test !occursin("Z_agl", rhs)
    @test length(rhs) < 60
end
