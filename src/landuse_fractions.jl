using NCDatasets
using Dates: unix2datetime, month

export reload_landuse_fractions!, set_landuse_data_path!, set_landuse_target_grid!

# NOTE on `data/landuse_wesely_conus_0p05.nc`
#
# Wesely (1989) land-use percentages at the native MCD12C1 resolution, held at
# source resolution rather than pre-binned to any simulation grid: the binning
# happens here at run time, once the target grid is known.
#
#   source     MODIS MCD12C1 v6.1, year 2016, `Land_Cover_Type_1_Percent`
#              (granule MCD12C1.A2016001.061.2022168010533, NASA Earthdata)
#   window     -130…-60 °E, 20…55 °N -- bounds are cell *edges* on exact
#              0.05° multiples, so the window aligns with MCD12C1 cells
#   grid       1400 × 700 × 11, cell size 0.05° (~5.6 km)
#   variable   `percent`, UInt8 0..100, dims (lon, lat, class)
#   classes    1 urban, 2 agricultural, 3 range, 4 deciduous, 5 coniferous,
#              6 mixed forest, 7 water, 8 barren, 9 wetland, 10 rangeag,
#              11 rocky shrubs -- matching the `wesleyLandUse` enum
#   attribs    `cell_size`, `window_{lon,lat}_{min,max}` as Float64; the grid
#              is rebuilt from these because the Float32 `lon`/`lat` variables
#              put cell edges ~1e-6° off and cost the remap four digits of
#              conservation
#
# To regenerate: read the 17 IGBP percentage bands, sum the bands mapping to
# each Wesely class (the mapping is many-to-one and covers all 17, so the
# result is integer and sums to exactly 100 per cell), slice the window, and
# write with lat increasing. Reading MCD12C1 needs GDAL with HDF4 support and
# NASA Earthdata credentials, so this is an offline step and its 1.2 GB
# granule is not part of the repo.
const _DEFAULT_LANDUSE_NC_PATH = joinpath(@__DIR__, "..", "data", "landuse_wesely_conus_0p05.nc")
const _LANDUSE_NC_PATH = Ref{String}(_DEFAULT_LANDUSE_NC_PATH)

# Target grid recorded by `DryDepositionGasFractional(domain)`. When set, the
# source is conservatively remapped onto it on first use; when unset, the
# lookup raises.
const _LANDUSE_TARGET = Ref{Union{Nothing, NamedTuple{
    (:lon, :lat, :dlon, :dlat),
    Tuple{Vector{Float64}, Vector{Float64}, Float64, Float64}}}}(nothing)

# Lazy-population guard. `_LANDUSE_READY[]` is the only thing the hot path
# reads; the lock is taken solely on a miss.
const _LANDUSE_READY = Ref{Bool}(false)
const _LANDUSE_LOCK = ReentrantLock()

# Concretely-typed cache fields, populated on first use. This mirrors
# the GasChem FastJX pattern (e.g. `const top_flux::SVector{18, Float32}`) so
# `landuse_frac_at` is type-stable and allocation-free in the hot RHS path.
# A `Ref{Union{Nothing, NamedTuple{...}}}` would force every cache access to
# narrow the union at runtime, causing per-call boxing and GC pressure
# (~2.8× slowdown observed in 9-day grid runs prior to this change).
const _LANDUSE_FRAC     = Ref{Array{Float32, 3}}()  # (lon, lat, class)
const _LANDUSE_LON_MIN  = Ref{Float64}()
const _LANDUSE_LON_STEP = Ref{Float64}()
const _LANDUSE_LAT_MIN  = Ref{Float64}()
const _LANDUSE_LAT_STEP = Ref{Float64}()
const _LANDUSE_N_LON    = Ref{Int}()
const _LANDUSE_N_LAT    = Ref{Int}()

# Read the source NetCDF. Accepts either `percent` (UInt8, 0..100) as written
# by the offline producer, or `fraction` (0..1) as used by test fixtures.
function _read_landuse_nc(path)
    ds = NCDataset(path)
    try
        lon = Float64.(ds["lon"][:])
        lat = Float64.(ds["lat"][:])
        # `lon`/`lat` are stored as Float32, which puts cell edges ~1e-6° off
        # the true window bounds and costs the remap several digits of
        # conservation. Rebuild them exactly from the Float64 window attributes
        # when the file carries them.
        if haskey(ds.attrib, "cell_size")
            c = Float64(ds.attrib["cell_size"])
            lon = [Float64(ds.attrib["window_lon_min"]) + c * (k - 0.5)
                   for k in 1:length(lon)]
            lat = [Float64(ds.attrib["window_lat_min"]) + c * (k - 0.5)
                   for k in 1:length(lat)]
        end
        frac = if haskey(ds, "percent")
            Float32.(ds["percent"][:, :, :]) ./ 100.0f0
        else
            Float32.(ds["fraction"][:, :, :])
        end
        return frac, lon, lat
    finally
        close(ds)
    end
end

# Cells of a uniform axis (first cell starting at `edge0`, width `step`) that
# the target interval [t_lo, t_hi] touches.
function _overlap_range(t_lo, t_hi, edge0, step, n)
    i_first = max(1, floor(Int, (t_lo - edge0) / step) + 1)
    i_last = min(n, ceil(Int, (t_hi - edge0) / step))
    return i_first:i_last
end

# Overlap of [t_lo, t_hi] with each touched cell. `lat = true` weights by true
# spherical area, whose meridional extent between two latitudes is proportional
# to sin(φ_hi) - sin(φ_lo); longitude overlap is simply linear.
function _overlap_1d(t_lo, t_hi, edge0, step, n; lat::Bool = false)
    rng = _overlap_range(t_lo, t_hi, edge0, step, n)
    isempty(rng) && return rng, Float64[]
    w = Vector{Float64}(undef, length(rng))
    for (k, i) in enumerate(rng)
        c_lo = edge0 + (i - 1) * step
        lo = max(t_lo, c_lo)
        hi = min(t_hi, c_lo + step)
        w[k] = lat ? sind(hi) - sind(lo) : hi - lo
    end
    return rng, w
end

# Conservatively remap class fractions from a uniform source grid onto uniform
# target cell centres: each source cell contributes in proportion to the area
# it shares with the target cell.
function _remap_landuse(src, src_lon, src_lat, tgt_lon, tgt_lat,
        tgt_dlon = _spacing(tgt_lon), tgt_dlat = _spacing(tgt_lat))
    n_class = size(src, 3)
    src_dlon = src_lon[2] - src_lon[1]
    src_dlat = src_lat[2] - src_lat[1]
    lon_edge0 = src_lon[1] - src_dlon / 2
    lat_edge0 = src_lat[1] - src_dlat / 2

    isnan(tgt_dlon) && error("Target lon cell size is unknown for a \
        single-cell axis; pass `dlon` to `set_landuse_target_grid!`.")
    isnan(tgt_dlat) && error("Target lat cell size is unknown for a \
        single-cell axis; pass `dlat` to `set_landuse_target_grid!`.")

    out = zeros(Float32, length(tgt_lon), length(tgt_lat), n_class)
    acc = Vector{Float64}(undef, n_class)
    for (jt, latc) in enumerate(tgt_lat)
        jr, jw = _overlap_1d(latc - tgt_dlat / 2, latc + tgt_dlat / 2,
            lat_edge0, src_dlat, length(src_lat); lat = true)
        isempty(jw) && _landuse_oob(tgt_lon[1], latc)
        for (it, lonc) in enumerate(tgt_lon)
            ir, iw = _overlap_1d(lonc - tgt_dlon / 2, lonc + tgt_dlon / 2,
                lon_edge0, src_dlon, length(src_lon))
            isempty(iw) && _landuse_oob(lonc, latc)
            fill!(acc, 0.0)
            total = 0.0
            for (kj, j) in enumerate(jr), (ki, i) in enumerate(ir)
                w = iw[ki] * jw[kj]
                total += w
                @inbounds for c in 1:n_class
                    acc[c] += w * src[i, j, c]
                end
            end
            @inbounds for c in 1:n_class
                out[it, jt, c] = Float32(acc[c] / total)
            end
        end
    end
    return out
end

# Without a target grid there is no cell to average over, and a bare lookup
# would hand a ~50 km cell whichever single ~5 km source pixel sits at its
# centroid. Refuse rather than return that: it is not a representative
# fraction for the cell, and nothing downstream can tell the difference.
@noinline function _no_target_grid_error(src_lon, src_lat)
    error("""
    Land-use fractions need the simulation grid before they can be used.
    The source is $(round(_spacing(src_lon), digits=4))° × \
    $(round(_spacing(src_lat), digits=4))°, so a bare lookup would sample one \
    source cell per query rather than area-average over your grid cell.

    Pass the domain when building the system:
        DryDepositionGasFractional(domain)
    or set the grid directly:
        set_landuse_target_grid!(lon_centres_deg, lat_centres_deg;
                                 dlon = ..., dlat = ...)""")
end

function _populate_landuse_cache!()
    path = _LANDUSE_NC_PATH[]
    isfile(path) || error(
        "Land-use source NetCDF not found at $path. " *
            "See the NOTE at the top of src/landuse_fractions.jl for how it is built.",
    )
    frac, lon, lat = _read_landuse_nc(path)

    target = _LANDUSE_TARGET[]
    target === nothing && _no_target_grid_error(lon, lat)
    _check_landuse_coverage(target.lon, target.lat, lon, lat,
        target.dlon, target.dlat)
    frac = _remap_landuse(frac, lon, lat, target.lon, target.lat,
        target.dlon, target.dlat)
    lon, lat = target.lon, target.lat

    _LANDUSE_FRAC[]     = frac
    _LANDUSE_LON_MIN[]  = lon[1]
    _LANDUSE_LON_STEP[] = length(lon) > 1 ? lon[2] - lon[1] : 1.0
    _LANDUSE_LAT_MIN[]  = lat[1]
    _LANDUSE_LAT_STEP[] = length(lat) > 1 ? lat[2] - lat[1] : 1.0
    _LANDUSE_N_LON[]    = length(lon)
    _LANDUSE_N_LAT[]    = length(lat)
    _LANDUSE_READY[]    = true
    return nothing
end

# Populate on first use. The target grid is not known until
# `DryDepositionGasFractional(domain)` has run, so this cannot be eager.
# Threads may reach it concurrently under `SolverStrangThreads`.
@noinline function _populate_landuse_locked!()
    lock(_LANDUSE_LOCK) do
        _LANDUSE_READY[] || _populate_landuse_cache!()
    end
    return nothing
end

@inline _ensure_landuse!() = _LANDUSE_READY[] ? nothing : _populate_landuse_locked!()

"""
    reload_landuse_fractions!()

Re-read the NetCDF at the current path into the cache. Useful after
regenerating the source NetCDF.
"""
function reload_landuse_fractions!()
    lock(_LANDUSE_LOCK) do
        _LANDUSE_READY[] = false
        _populate_landuse_cache!()
    end
    return nothing
end

"""
    set_landuse_target_grid!(lon_deg, lat_deg; dlon, dlat)

Record the simulation grid (cell centres, DEGREES) that the source data should
be remapped onto, and invalidate any cached fractions. Called by
`DryDepositionGasFractional(domain)`; pass `nothing` to use the source grid
directly.

`dlon`/`dlat` give the cell size. They default to the centre spacing, which is
undefined for a single-cell axis — pass them explicitly in that case, as the
domain-aware constructor does.
"""
function set_landuse_target_grid!(lon_deg, lat_deg;
        dlon = _spacing(lon_deg), dlat = _spacing(lat_deg))
    lock(_LANDUSE_LOCK) do
        _LANDUSE_TARGET[] = (lon = collect(Float64, lon_deg),
            lat = collect(Float64, lat_deg),
            dlon = Float64(dlon), dlat = Float64(dlat))
        _LANDUSE_READY[] = false
    end
    return nothing
end

# Centre spacing of a uniform axis. NaN for a single cell, where the spacing
# cannot be recovered from the centres alone.
_spacing(c) = length(c) > 1 ? Float64(c[2]) - Float64(c[1]) : NaN

function set_landuse_target_grid!(::Nothing)
    lock(_LANDUSE_LOCK) do
        _LANDUSE_TARGET[] = nothing
        _LANDUSE_READY[] = false
    end
    return nothing
end

# Coverage of the source file, read without populating the cache. Used by
# `DryDepositionGasFractional(domain)` so an unsupported domain fails at
# construction rather than mid-solve.
function landuse_coverage()
    path = _LANDUSE_NC_PATH[]
    isfile(path) || error("Land-use source NetCDF not found at $path.")
    ds = NCDataset(path)
    try
        lon = Float64.(ds["lon"][:])
        lat = Float64.(ds["lat"][:])
        dlon = length(lon) > 1 ? lon[2] - lon[1] : 0.0
        dlat = length(lat) > 1 ? lat[2] - lat[1] : 0.0
        return (lon = (lon[1] - dlon / 2, lon[end] + dlon / 2),
            lat = (lat[1] - dlat / 2, lat[end] + dlat / 2))
    finally
        close(ds)
    end
end

"""
    record_landuse_grid!(domain)

Validate `domain` against the source coverage and record its lon/lat grid as
the remap target. Called by `DryDepositionGasFractional(domain)`.
"""
function record_landuse_grid!(domain)
    lon, lat, dlon, dlat = _domain_lonlat_deg(domain)
    cov = landuse_coverage()
    _coverage_error(_edges(lon, dlon), cov.lon, "lon")
    _coverage_error(_edges(lat, dlat), cov.lat, "lat")
    set_landuse_target_grid!(lon, lat; dlon = dlon, dlat = dlat)
    return nothing
end

# lon/lat cell centres and cell sizes of a DomainInfo, in degrees.
# `EarthSciMLBase.grid` returns one range per partial independent variable, in
# radians; taking the size from the range's `step` keeps it defined even when
# an axis has a single cell.
function _domain_lonlat_deg(domain)
    g = EarthSciMLBase.grid(domain)
    pvs = String.(Symbol.(EarthSciMLBase.pvars(domain)))
    ilon = findfirst(v -> occursin("lon", v), pvs)
    ilat = findfirst(v -> occursin("lat", v), pvs)
    (ilon === nothing || ilat === nothing) &&
        error("Domain has no lon/lat coordinates; got $(pvs).")
    return rad2deg.(collect(Float64, g[ilon])), rad2deg.(collect(Float64, g[ilat])),
    rad2deg(Float64(step(g[ilon]))), rad2deg(Float64(step(g[ilat])))
end

# Outer edges of a uniform axis given its cell centres and cell size.
function _edges(c, d = _spacing(c))
    isnan(d) && (d = 0.0)
    return (c[1] - d / 2, c[end] + d / 2)
end

function _check_landuse_coverage(tgt_lon, tgt_lat, src_lon, src_lat,
        tdlon = _spacing(tgt_lon), tdlat = _spacing(tgt_lat))
    _coverage_error(_edges(tgt_lon, tdlon), _edges(src_lon), "lon")
    _coverage_error(_edges(tgt_lat, tdlat), _edges(src_lat), "lat")
    return nothing
end

function _coverage_error(t, s, axis)
    # Coordinates are stored as Float32, so allow ~1e-4° of round-trip slack.
    # A real overrun is at least half a source cell (0.025°).
    (t[1] >= s[1] - 1.0e-4 && t[2] <= s[2] + 1.0e-4) && return nothing
    error("""
    Domain $axis range $(round(t[1], digits=4))°…$(round(t[2], digits=4))° extends beyond \
    the land-use coverage $(round(s[1], digits=4))°…$(round(s[2], digits=4))°.
    Regenerate $(basename(_LANDUSE_NC_PATH[])) with a wider window \
    (see the NOTE in src/landuse_fractions.jl), or restrict the domain.""")
end

@noinline function _landuse_oob(lon_deg, lat_deg)
    cov = try
        landuse_coverage()
    catch
        nothing
    end
    where = cov === nothing ? "" :
            " Coverage is $(round(cov.lon[1], digits=4))°…$(round(cov.lon[2], digits=4))° \
lon, $(round(cov.lat[1], digits=4))°…$(round(cov.lat[2], digits=4))° lat."
    error("Land-use lookup at lon=$(round(lon_deg, digits=4))°, \
lat=$(round(lat_deg, digits=4))° falls outside the available data.$where")
end

"""
    set_landuse_data_path!(path)

Override the NetCDF path used by `landuse_frac_at` and reload immediately.
Primarily used by tests to point at a synthetic fixture.
"""
function set_landuse_data_path!(path::AbstractString)
    lock(_LANDUSE_LOCK) do
        _LANDUSE_NC_PATH[] = String(path)
        _LANDUSE_READY[] = false
        _populate_landuse_cache!()
    end
    return nothing
end

# Area fraction of Wesely (1989) land-use class `class_index` at the
# given longitude/latitude (RADIANS, matching `DomainInfo`). Indices
# follow the `wesleyLandUse` enum:
#   1 = urban, 2 = agricultural, 3 = range, 4 = deciduous, 5 = coniferous,
#   6 = mixed forest, 7 = water, 8 = barren, 9 = wetland, 10 = rangeag,
#   11 = rocky shrubs.
#
# Nearest-cell lookup on the cached uniform grid — the simulation grid when
# `DryDepositionGasFractional(domain)` recorded one, otherwise the source
# grid. Query points outside the grid raise rather than return a boundary
# cell, so an unsupported domain cannot pass silently.
function landuse_frac_at(lon_rad, lat_rad, class_index)::Float64
    _ensure_landuse!()
    lon_deg = rad2deg(lon_rad)
    lat_deg = rad2deg(lat_rad)
    i_lon = round(Int, (lon_deg - _LANDUSE_LON_MIN[]) / _LANDUSE_LON_STEP[]) + 1
    i_lat = round(Int, (lat_deg - _LANDUSE_LAT_MIN[]) / _LANDUSE_LAT_STEP[]) + 1
    (1 <= i_lon <= _LANDUSE_N_LON[]) || _landuse_oob(lon_deg, lat_deg)
    (1 <= i_lat <= _LANDUSE_N_LAT[]) || _landuse_oob(lon_deg, lat_deg)
    @inbounds Float64(_LANDUSE_FRAC[][i_lon, i_lat, class_index])
end
# Registered as symbolic so that `lon`/`lat` are preserved as parameter
# references in the System equations (otherwise the symbolic engine
# might fold based on the lookup logic).
@register_symbolic landuse_frac_at(lon, lat, class_index)
landuse_frac_at(::DynamicQuantities.Quantity, ::DynamicQuantities.Quantity, ::Any) = 0.0
ModelingToolkit.get_unit(::typeof(landuse_frac_at)) = 1.0

# Wesely (1989) season index from absolute Unix time, assuming
# Northern-Hemisphere mid-latitudes:
#   Jun–Aug → 1 (Midsummer)
#   Sep–Oct → 2 (Autumn)
#   Nov     → 3 (LateAutumn)
#   Dec–Feb → 4 (Winter)
#   Mar–May → 5 (Transitional)
#
# Called from `couple2(::DryDepositionGasCoupler, ::GEOSFPCoupler)` in
# `ext/EarthSciDataExt.jl` as `season_at(gp.t_ref + t)`. For
# Southern-Hemisphere or tropical work, override `season` manually instead of
# binding to this function.
#
# TODO: lat-aware seasons for global use.
function season_at(t_abs_seconds)
    m = month(unix2datetime(t_abs_seconds))
    return m == 6 || m == 7 || m == 8  ? Int(wesleyMidsummer)    :
           m == 9 || m == 10           ? Int(wesleyAutumn)       :
           m == 11                     ? Int(wesleyLateAutumn)   :
           m == 12 || m == 1 || m == 2 ? Int(wesleyWinter)       :
                                         Int(wesleyTransitional)
end
@register_symbolic season_at(t_abs_seconds)
season_at(::DynamicQuantities.Quantity) = 1
ModelingToolkit.get_unit(::typeof(season_at)) = 1.0
