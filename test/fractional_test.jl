@testsnippet FractionalSetup begin
    using AtmosphericDeposition
    using AtmosphericDeposition: FractionalWesleyVd, landuse_frac_at,
        DryDepGasFractional, reload_landuse_fractions!, set_landuse_data_path!,
        set_landuse_target_grid!, landuse_coverage
    using Test, DynamicQuantities, ModelingToolkit
    using StaticArrays
    using NCDatasets

    const DEFAULT_LANDUSE_NC = joinpath(
        @__DIR__, "..", "data", "landuse_wesely_conus_0p05.nc")

    # Restore the bundled source on a known grid. A target grid is mandatory —
    # without one there is no cell to area-average over — so tests that only
    # care about lookups declare the 0.625° × 0.5° CONUS grid.
    const TEST_LON = collect(-125.625:0.625:-66.25)
    const TEST_LAT = collect(24.5:0.5:49.5)

    function reset_landuse!()
        set_landuse_target_grid!(TEST_LON, TEST_LAT)
        set_landuse_data_path!(DEFAULT_LANDUSE_NC)
    end
end

@testitem "landuse_frac_at bundled CONUS fractions" setup=[FractionalSetup] begin
    # Bundled NetCDF is MODIS MCD12C1 v6.1 re-classed to Wesely classes at the
    # native 0.05°. Verify the physical sanity of fractions at well-known
    # locations: every cell must sum to 1.0, and the dominant class must match
    # the geography. Tolerances are loose — we want catastrophic-bug regression
    # protection, not a pixel-exact comparison to MODIS.
    reset_landuse!()

    function frac(lon_deg, lat_deg)
        [landuse_frac_at(deg2rad(lon_deg), deg2rad(lat_deg), i) for i in 1:11]
    end

    # All cells in the CONUS interior sum to exactly 1.0 (the area weighting in
    # FractionalWesleyVd relies on this; any leftover is silently dropped).
    for (lon, lat) in [(-100.0, 40.0), (-90.0, 35.0), (-118.0, 34.0), (-75.0, 40.0)]
        @test sum(frac(lon, lat)) ≈ 1.0 atol = 1.0e-5
    end

    # Class-1 (urban) dominant over LA basin.
    @test frac(-118.25, 34.05)[1] > 0.5
    # Class-2 (agricultural) dominant over Iowa.
    @test frac(-93.5, 42.0)[2] > 0.5
    # Class-7 (water) ≈ 1.0 over Gulf of Mexico.
    @test frac(-88.0, 27.0)[7] > 0.95
end

@testitem "landuse conservative regrid onto a coarse grid" setup=[FractionalSetup] begin
    # Remapping the 0.05° source onto a 5° × 4° grid must area-average, not
    # point-sample: a coastal cell has to come out genuinely mixed. Point
    # sampling would return a single source cell, i.e. 0.0/1.0 fractions.
    try
        set_landuse_target_grid!(collect(-127.5:5.0:-62.5), collect(22.0:4.0:50.0))
        reload_landuse_fractions!()

        frac(lo, la) = [landuse_frac_at(deg2rad(lo), deg2rad(la), i) for i in 1:11]

        # Mid-Atlantic coast: land and water in the same 5° × 4° box.
        v = frac(-77.5, 38.0)
        @test sum(v) ≈ 1.0 atol = 1.0e-5
        @test count(>(0.001), v) >= 4        # genuinely mixed
        @test 0.02 < v[7] < 0.9              # some water, not all water

        # Conservation holds on every cell of the coarse grid.
        for lo in -127.5:5.0:-62.5, la in 22.0:4.0:50.0
            @test sum(frac(lo, la)) ≈ 1.0 atol = 1.0e-5
        end
    finally
        reset_landuse!()
    end
end

@testitem "landuse remap conserves the area integral" setup=[FractionalSetup] begin
    # The defining property of a conservative remap: the area-weighted integral
    # of each class is unchanged. Area on a uniform lon/lat grid is
    # proportional to dlon * (sin(lat_hi) - sin(lat_lo)).
    using AtmosphericDeposition: _read_landuse_nc, _remap_landuse
    src, slon, slat = _read_landuse_nc(
        joinpath(@__DIR__, "..", "data", "landuse_wesely_conus_0p05.nc"))

    function integral(f, lon, lat, dlo, dla)
        a = [dlo * (sind(la + dla / 2) - sind(la - dla / 2)) for _ in lon, la in lat]
        [sum(a .* @view f[:, :, c]) for c in 1:size(f, 3)]
    end
    Is = integral(src, slon, slat, 0.05, 0.05)

    # Targets that exactly tile the source window, coarse and very coarse.
    for (dlo, dla) in ((2.0, 2.5), (3.5, 5.0))
        tlon = collect((-130.0 + dlo / 2):dlo:(-60.0 - dlo / 2))
        tlat = collect((20.0 + dla / 2):dla:(55.0 - dla / 2))
        out = _remap_landuse(src, slon, slat, tlon, tlat, dlo, dla)
        It = integral(out, tlon, tlat, dlo, dla)
        @test maximum(abs.(It .- Is) ./ Is) < 1.0e-6
    end

    # A uniform field must survive exactly, at any target resolution.
    flat = fill(0.25f0, size(src, 1), size(src, 2), 4)
    o = _remap_landuse(flat, slon, slat, collect(-129.0:2.0:-61.0),
        collect(21.25:2.5:53.75), 2.0, 2.5)
    @test maximum(abs.(o .- 0.25f0)) == 0.0
end

@testitem "landuse single-cell axis needs an explicit cell size" setup=[FractionalSetup] begin
    # Cell size cannot be recovered from the centres of a one-cell axis. Guessing
    # it silently produced a point sample where a whole-window average was meant.
    using AtmosphericDeposition: _read_landuse_nc, _remap_landuse
    src, slon, slat = _read_landuse_nc(
        joinpath(@__DIR__, "..", "data", "landuse_wesely_conus_0p05.nc"))

    @test_throws ErrorException _remap_landuse(src, slon, slat, [-95.0], [37.5])

    # Given the size, one cell spanning the window equals the window mean.
    o = _remap_landuse(src, slon, slat, [-95.0], [37.5], 70.0, 35.0)
    a = [0.05 * (sind(la + 0.025) - sind(la - 0.025)) for _ in slon, la in slat]
    tot = sum(a)
    expected = [sum(a .* @view src[:, :, c]) / tot for c in 1:11]
    @test maximum(abs.(vec(o) .- expected)) < 1.0e-6
end

@testitem "landuse without a target grid raises" setup=[FractionalSetup] begin
    # No grid means no cell to average over, so a lookup would hand a coarse
    # cell whichever single 0.05° source pixel sits at its centroid. That is
    # less representative than the pre-baked file this source replaced, and
    # nothing downstream could tell, so it must fail loudly instead.
    try
        set_landuse_target_grid!(nothing)
        set_landuse_data_path!(DEFAULT_LANDUSE_NC)   # populates -> must raise
        @test false  # unreachable
    catch e
        msg = sprint(showerror, e)
        @test occursin("simulation grid", msg)
        @test occursin("DryDepositionGasFractional(domain)", msg)
    finally
        reset_landuse!()
    end

    # And it recovers as soon as a grid is declared.
    @test 0.0 <= landuse_frac_at(deg2rad(-100.0), deg2rad(40.0), 2) <= 1.0
end

@testitem "landuse out-of-coverage raises instead of clamping" setup=[FractionalSetup] begin
    # Silently clamping an out-of-range query to the nearest boundary cell is
    # what turns an unsupported domain into plausible-looking wrong numbers.
    reset_landuse!()
    cov = landuse_coverage()
    @test cov.lon[1] ≈ -130.0 atol = 1.0e-3
    @test cov.lat[2] ≈ 55.0 atol = 1.0e-3

    # Inside the declared grid still resolves. Note the lookup is bounded by
    # the target grid, not by the wider source window: -129° is within the
    # source but outside the CONUS grid declared here, so it raises too.
    @test landuse_frac_at(deg2rad(-120.0), deg2rad(40.0), 7) >= 0.0
    @test_throws ErrorException landuse_frac_at(deg2rad(-129.0), deg2rad(21.0), 7)

    # Outside on any side throws, and names the offending coordinate.
    @test_throws ErrorException landuse_frac_at(deg2rad(-145.0), deg2rad(40.0), 1)
    @test_throws ErrorException landuse_frac_at(deg2rad(-100.0), deg2rad(5.0), 1)
    @test_throws ErrorException landuse_frac_at(deg2rad(-40.0), deg2rad(40.0), 1)
    @test_throws ErrorException landuse_frac_at(deg2rad(-100.0), deg2rad(70.0), 1)

    err = try
        landuse_frac_at(deg2rad(-145.0), deg2rad(40.0), 1)
        nothing
    catch e
        sprint(showerror, e)
    end
    @test occursin("-145.0", err)
end

@testitem "landuse_frac_at lookup varies with lat/lon" setup=[FractionalSetup] begin
    # Synthesise a NetCDF with class 1 = 1.0 in the western half, class 7 =
    # 1.0 in the eastern half, so the lookup must actually use lon. Bin to
    # the same grid as the bundled file (96 × 51 cells).
    mktempdir() do dir
        path = joinpath(dir, "split.nc")
        n_lon, n_lat, n_class = 96, 51, 11
        lon = collect(Float32, range(-125.625f0; step = 0.625f0, length = n_lon))
        lat = collect(Float32, range(24.5f0; step = 0.5f0, length = n_lat))
        frac = zeros(Float32, n_lon, n_lat, n_class)
        for j in 1:n_lat, i in 1:n_lon
            if lon[i] < -96.0
                frac[i, j, 1] = 1.0   # urban (western half)
            else
                frac[i, j, 7] = 1.0   # water (eastern half)
            end
        end
        NCDataset(path, "c") do ds
            defDim(ds, "lon", n_lon); defDim(ds, "lat", n_lat); defDim(ds, "class", n_class)
            defVar(ds, "lon", lon, ("lon",))
            defVar(ds, "lat", lat, ("lat",))
            defVar(ds, "fraction", frac, ("lon", "lat", "class"))
        end

        try
            set_landuse_target_grid!(collect(lon), collect(lat))
            set_landuse_data_path!(path)
            # Western point — urban
            @test landuse_frac_at(deg2rad(-120.0), deg2rad(40.0), 1) == 1.0
            @test landuse_frac_at(deg2rad(-120.0), deg2rad(40.0), 7) == 0.0
            # Eastern point — water
            @test landuse_frac_at(deg2rad(-70.0), deg2rad(40.0), 1) == 0.0
            @test landuse_frac_at(deg2rad(-70.0), deg2rad(40.0), 7) == 1.0
            # Out-of-range query (south of grid) raises rather than clamping.
            @test_throws ErrorException landuse_frac_at(
                deg2rad(-120.0), deg2rad(10.0), 1)
        finally
            reset_landuse!()
        end
    end
end

@testitem "FractionalWesleyVd single-class equality" setup=[FractionalSetup] begin
    # For each Wesely class i, putting all weight on i must give exactly the
    # single-class deposition velocity 1/(RaRb + Rc_i). This is the degenerate
    # case in which the area-weighting cancels entirely.
    G, Ts, θ, iSeason = 800.0, 25.0, 0.0, 1
    RaRb = 50.0

    for i in 1:11
        fractions = zeros(11)
        fractions[i] = 1.0

        rc_single = WesleySurfaceResistance(
            AtmosphericDeposition.O3Data, G, Ts, θ, iSeason, i,
            false, false, false, true
        )
        vd_frac = FractionalWesleyVd(
            RaRb, AtmosphericDeposition.O3Data, G, Ts, θ, iSeason,
            fractions...,
            false, false, false, true
        )
        @test vd_frac ≈ 1 / (RaRb + rc_single)
    end
end

@testitem "FractionalWesleyVd two-class arithmetic" setup=[FractionalSetup] begin
    # 50% agricultural, 50% mixed forest — area-weighted Vd check for O3.
    G, Ts, θ, iSeason = 800.0, 25.0, 0.0, 1
    RaRb = 50.0
    rc_ag = WesleySurfaceResistance(
        AtmosphericDeposition.O3Data, G, Ts, θ, iSeason, 2,
        false, false, false, true
    )
    rc_mf = WesleySurfaceResistance(
        AtmosphericDeposition.O3Data, G, Ts, θ, iSeason, 6,
        false, false, false, true
    )
    vd_expected = 0.5 / (RaRb + rc_ag) + 0.5 / (RaRb + rc_mf)

    vd_frac = FractionalWesleyVd(
        RaRb, AtmosphericDeposition.O3Data, G, Ts, θ, iSeason,
        0.0, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0,
        false, false, false, true
    )
    @test vd_frac ≈ vd_expected
end

@testitem "FractionalWesleyVd three-class arithmetic" setup=[FractionalSetup] begin
    # Realistic mixed cell for HNO3: 50% deciduous, 30% coniferous, 20% urban.
    G, Ts, θ, iSeason = 800.0, 25.0, 0.0, 1
    RaRb = 50.0
    rc_dec = WesleySurfaceResistance(
        AtmosphericDeposition.HNO3Data, G, Ts, θ, iSeason, 4,
        false, false, false, false
    )
    rc_con = WesleySurfaceResistance(
        AtmosphericDeposition.HNO3Data, G, Ts, θ, iSeason, 5,
        false, false, false, false
    )
    rc_urb = WesleySurfaceResistance(
        AtmosphericDeposition.HNO3Data, G, Ts, θ, iSeason, 1,
        false, false, false, false
    )
    vd_expected = 0.5 / (RaRb + rc_dec) + 0.3 / (RaRb + rc_con) +
        0.2 / (RaRb + rc_urb)

    vd_frac = FractionalWesleyVd(
        RaRb, AtmosphericDeposition.HNO3Data, G, Ts, θ, iSeason,
        0.2, 0.0, 0.0, 0.5, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        false, false, false, false
    )
    @test vd_frac ≈ vd_expected
end

@testitem "FractionalWesleyVd is below the Rc-blending form" setup=[FractionalSetup] begin
    # Direction guard. Blending 1/Rc first and adding a single Ra+Rb afterwards
    # is the Ra+Rb -> 0 limit of the correct expression, and by concavity of
    # c -> c/(1+ac) it always gives Vd >= this one, with equality only when
    # every class shares one Rc. A regression to that form must fail here.
    G, Ts, θ, iSeason = 800.0, 25.0, 0.0, 1
    # 96% water + 4% wetland: the widest Rc spread in the bundled CONUS data.
    fr = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.96, 0.0, 0.04, 0.0, 0.0)
    rc = [WesleySurfaceResistance(AtmosphericDeposition.O3Data, G, Ts, θ, iSeason, i,
        false, false, false, true) for i in 1:11]

    for RaRb in (10.0, 50.0, 100.0, 300.0, 500.0)
        vd_new = FractionalWesleyVd(
            RaRb, AtmosphericDeposition.O3Data, G, Ts, θ, iSeason,
            fr..., false, false, false, true
        )
        vd_old = 1 / (RaRb + 1 / sum(fr[i] / rc[i] for i in 1:11))
        @test vd_new < vd_old
    end

    # Equal Rc across the populated classes => the two forms coincide.
    fr_eq = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0)  # water, barren
    if rc[7] ≈ rc[8]
        vd_new = FractionalWesleyVd(
            100.0, AtmosphericDeposition.O3Data, G, Ts, θ, iSeason,
            fr_eq..., false, false, false, true
        )
        @test vd_new ≈ 1 / (100.0 + rc[7])
    end
end

@testitem "DryDepGasFractional unit" setup=[FractionalSetup] begin
    @parameters T [unit = u"K"]
    @parameters z [unit = u"m"]
    @parameters z₀ [unit = u"m"]
    @parameters u_star [unit = u"m/s"]
    @parameters L [unit = u"m"]
    @parameters ρA [unit = u"kg*m^-3"]
    @parameters G [unit = u"W*m^-2"]
    @parameters θ
    @parameters iSeason, lev

    # Concrete fractions: f_mixedforest = 1.0, others = 0.0.
    # GasData overload — this is the active path used by the
    # `DryDepositionGasFractional` broadcast (via `datas = [...]`).
    @test ModelingToolkit.get_unit(
        DryDepGasFractional(
            lev, z, z₀, u_star, L, ρA,
            AtmosphericDeposition.So2Data,
            G, T, θ, iSeason,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            false, false, true, false
        )
    ) == u"m/s"
end

@testitem "landuse lazy population is thread-safe" setup=[FractionalSetup] begin
    # SolverStrangThreads means several threads can reach the first lookup at
    # once. The population must happen once and all threads must agree.
    try
        set_landuse_target_grid!(collect(-120.0:1.0:-80.0), collect(30.0:1.0:45.0))
        # Deliberately do NOT call reload_landuse_fractions!: leave the cache
        # cold so the threads race on the lazy path.
        n = max(4, Threads.nthreads())
        out = Vector{Float64}(undef, n)
        Threads.@threads for k in 1:n
            out[k] = landuse_frac_at(deg2rad(-100.0), deg2rad(40.0), 2)
        end
        @test all(==(out[1]), out)
        @test 0.0 <= out[1] <= 1.0
    finally
        reset_landuse!()
    end
end

@testitem "DryDepositionGasFractional(domain) validates coverage" setup=[FractionalSetup] begin
    using EarthSciMLBase, Dates

    dom(lonrange, latrange) = DomainInfo(
        DateTime(2016, 2, 1), DateTime(2016, 2, 2);
        lonrange = lonrange, latrange = latrange, levrange = 1:3
    )

    try
        # Inside coverage: constructs, and records the grid so lookups land on
        # the simulation cells.
        @test DryDepositionGasFractional(
            dom(deg2rad(-120.0):deg2rad(1.0):deg2rad(-80.0),
                deg2rad(30.0):deg2rad(1.0):deg2rad(45.0))
        ) isa ModelingToolkit.AbstractSystem

        # West of -130: must fail at construction, not mid-solve.
        @test_throws ErrorException DryDepositionGasFractional(
            dom(deg2rad(-145.0):deg2rad(1.0):deg2rad(-80.0),
                deg2rad(30.0):deg2rad(1.0):deg2rad(45.0))
        )
        # North of 55.
        @test_throws ErrorException DryDepositionGasFractional(
            dom(deg2rad(-120.0):deg2rad(1.0):deg2rad(-80.0),
                deg2rad(30.0):deg2rad(1.0):deg2rad(60.0))
        )
    finally
        reset_landuse!()
    end
end

@testitem "DryDepositionGasFractional constructor" setup=[FractionalSetup] begin
    reset_landuse!()
    sys = DryDepositionGasFractional()

    @test sys isa ModelingToolkit.AbstractSystem

    # Fractions are now ordinary parameters (set per cell by the coupler),
    # not system unknowns. This matches the scalar-landuse baseline's
    # symbolic shape: only v_* / k_* are unknowns; everything else is a
    # parameter. lon/lat live on the GEOS-FP side and aren't dep-system
    # parameters anymore.
    param_str = string.(parameters(sys))
    @test any(contains.(param_str, "f_urban"))
    @test any(contains.(param_str, "f_agricultural"))
    @test any(contains.(param_str, "f_mixedforest"))
    @test any(contains.(param_str, "f_rockyshrubs"))
    @test any(contains.(param_str, "season"))
    @test !any(contains.(param_str, "landuse"))
    @test !any(contains.(param_str, "lon"))
    @test !any(contains.(param_str, "lat"))

    var_str = string.(unknowns(sys))
    @test !any(contains.(var_str, "f_urban"))
end
