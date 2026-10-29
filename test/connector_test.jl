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

    model = couple(SuperFast(), DryDepositionGasFractional())
    sys = convert(System, model)
    eqs = string(equations(sys))

    # Verify that GasChem species are coupled to dry deposition rate constants
    @test contains(eqs, "SuperFast₊DryDepositionGasFractional_k_SO2")
    @test contains(eqs, "SuperFast₊DryDepositionGasFractional_k_NH3")
    @test contains(eqs, "SuperFast₊DryDepositionGasFractional_k_HNO3")
    @test contains(eqs, "SuperFast₊DryDepositionGasFractional_k_NO2")
    @test contains(eqs, "SuperFast₊DryDepositionGasFractional_k_O3")
    @test contains(eqs, "SuperFast₊DryDepositionGasFractional_k_H2O2")
    # CH2O uses the GEOS-Chem-aligned k_CH2O channel (CH2OData: Hstar=3e3, f0=1.0),
    # not the legacy Wesely-table k_HCHO (f0=0) - see the coupling in ext/GasChemExt.jl.
    @test contains(eqs, "SuperFast₊DryDepositionGasFractional_k_CH2O")
    @test !contains(eqs, "SuperFast₊DryDepositionGasFractional_k_HCHO")
    # NO is not dry-deposited, matching GEOS-Chem.
    @test !occursin(r"SuperFast₊DryDepositionGasFractional_k_NO(?!\w)", eqs)
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
    # The interpolation syntax on the RHS varies across EarthSciData versions,
    # so only check that cloudFrac is coupled to the GEOS-FP CLOUD variable.
    @test occursin(r"WetDeposition₊cloudFrac\(t\) ~ [^,]*GEOSFP₊A3cld₊CLOUD", eqs)
end

@testitem "EarthSciDataExt fractional gas dry deposition" setup = [ConnectorSetup] begin
    # Verify the new couple2 method binds the land-use fractions to a lookup at
    # the GEOS-FP lon/lat, binds season to the GEOS-FP date, and that the
    # surface-meteorology bindings carry over from the scalar-landuse coupler.
    model = couple(
        GEOSFP("4x5", domain),
        DryDepositionGasFractional(),
    )
    sys = convert(System, model)
    eqs = string(equations(sys)) * "\n" * string(observed(sys))

    # New (fractional-specific) bindings. The coupler passes the GEOS-FP
    # coordinates straight to the lookup; the deposition system has no lon/lat.
    for (i, f) in enumerate(
            (
                :f_urban, :f_agricultural, :f_range, :f_deciduous, :f_coniferous,
                :f_mixedforest, :f_water, :f_barren, :f_wetland, :f_rangeag, :f_rockyshrubs,
            )
        )
        @test contains(
            eqs,
            "DryDepositionGasFractional₊$(f)(t) ~ " *
                "AtmosphericDeposition.landuse_frac_at(GEOSFP₊lon, GEOSFP₊lat, $i)"
        )
    end
    @test contains(
        eqs,
        "DryDepositionGasFractional₊season(t) ~ " *
            "AtmosphericDeposition.season_at(GEOSFP₊t_ref + t)"
    )

    # Carried over from the scalar-landuse coupler.
    @test contains(eqs, "GEOSFP₊A1₊USTAR")
    @test contains(eqs, "GEOSFP₊A1₊TS")
    @test contains(eqs, "GEOSFP₊A1₊SWGDN")

    # No scalar landuse parameter survives (the fractional system never
    # declared one — sanity guard against silent regression).
    @test !contains(eqs, "DryDepositionGasFractional₊landuse")
end

@testitem "dry-dep reference height is the level-1 midpoint (no Z_agl in RHS)" begin
    using EarthSciMLBase, EarthSciData, Dates, ModelingToolkit, Symbolics
    using EarthSciMLBase: DomainInfo
    t0 = DateTime(2016, 3, 10)
    dom = DomainInfo(
        t0, t0 + Hour(3) + Hour(36);
        lonrange = deg2rad(-88):deg2rad(0.625):deg2rad(-86),
        latrange = deg2rad(32):deg2rad(0.5):deg2rad(33.5),
        levrange = 1:2, u_proto = zeros(Float64, 1, 1, 1, 1)
    )
    gfp = EarthSciData.GEOSFP("0.25x0.3125", dom)
    # `DryDepositionGas` here would resolve to the unrelated Seinfeld-Pandis
    # ch. 19 component; the fractional system carries `DryDepositionGasCoupler`.
    d = DryDepositionGasFractional()
    cs = EarthSciMLBase.couple2(EarthSciMLBase.get_coupletype(d)(d),
                                EarthSciMLBase.get_coupletype(gfp)(gfp))
    zeq = [e for e in cs.eqs if endswith(string(e.lhs), "₊z(t)")]
    @test length(zeq) == 1
    rhs = only(zeq).rhs
    # The reference height is computed from level-1 and 2 m fields: referencing the
    # height-above-ground field would be 0 at lev=1 on EarthSciData 0.16 and would
    # inline the multi-layer interpolation cascade into the deposition RHS.
    @test !occursin("Z_agl", string(rhs))

    # Standard surface conditions put the level-1 midpoint at 64.07 m, matching the
    # hydrostatic half-layer depth (ΔP/2)·Rd·Tv/(g·P̄).
    met = Dict("I3₊PS" => 101325.0, "I3₊T" => 288.0, "I3₊QV" => 0.01,
        "A1₊T2M" => 290.0, "A1₊QV2M" => 0.01)
    function value_for(v)
        name = string(Symbolics.getname(v))
        for (k, val) in met
            endswith(name, k) && return val
        end
        return ModelingToolkit.getdefault(v)
    end
    subs = Dict(v => value_for(v) for v in Symbolics.get_variables(rhs))
    z = Float64(eval(Symbolics.toexpr(Symbolics.unwrap(substitute(rhs, subs)))))
    @test z ≈ 64.07 rtol = 1e-3
end
