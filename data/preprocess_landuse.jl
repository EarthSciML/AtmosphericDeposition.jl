#!/usr/bin/env julia
# Downloads MODIS MCD12C1 IGBP land cover for the configured year,
# re-classes IGBP (17 classes) to Wesely (1989) (11 classes), extracts the
# CONUS-plus-buffer window at the native 0.05° resolution, and writes
# `data/landuse_wesely_conus_0p05.nc`.
#
# Binning to a simulation grid happens at run time in
# `src/landuse_fractions.jl`, which conservatively remaps this file onto
# whatever grid the domain specifies. This script is dev-only: run it to
# regenerate the bundled source, never as part of a simulation.
#
# Prerequisites
# -------------
# 1. NASA Earthdata Login account — register (free) at
#    https://urs.earthdata.nasa.gov/
# 2. Credentials in ~/.netrc:
#       machine urs.earthdata.nasa.gov
#         login YOUR_USERNAME
#         password YOUR_PASSWORD
#    Then: chmod 600 ~/.netrc
# 3. Julia 1.10+, environment instantiated:
#       julia --project=data -e 'using Pkg; Pkg.instantiate()'
#
# Run
# ---
#   julia --project=data data/preprocess_landuse.jl
#
# Notes
# -----
# MCD12C1 is a Climate Modeling Grid (CMG) product at 0.05° resolution
# (7200 × 3600 cells, global). The `Land_Cover_Type_1_Percent` subdataset
# provides per-cell area fractions of all 17 IGBP classes in 0–100 percent.
# We re-class the 17 IGBP percentages to 11 Wesely (1989) classes using the
# `IGBP_TO_WESELY` table below, keeping the native 0.05° resolution.
# Percentages are integers summing to 100 per cell, so the re-class is exact
# and the result stores losslessly as UInt8.

using ArchGDAL
using NCDatasets
using Downloads
using Printf

# ─────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────

const MCD12C1_YEAR = 2016
# NASA migrated MCD12C1 v6.1 to Earthdata Cloud (March 2024). The legacy
# e4ftl01.cr.usgs.gov path returns 404. The current cloud URL embeds the
# full granule ID as a directory.
# If the URL below 404s, query the CMR catalogue for the current granule:
#   curl -s "https://cmr.earthdata.nasa.gov/search/granules.json?short_name=MCD12C1&version=061&temporal=$(MCD12C1_YEAR)-01-01,$(MCD12C1_YEAR)-12-31&page_size=1"
const MCD12C1_GRANULE = "MCD12C1.A$(MCD12C1_YEAR)001.061.2022168010533"
const MCD12C1_URL = "https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/" *
    "MCD12C1.061/$(MCD12C1_GRANULE)/$(MCD12C1_GRANULE).hdf"
const MCD12C1_PATH = joinpath(@__DIR__, "mcd12c1_$(MCD12C1_YEAR).hdf")

const OUTPUT_NC = joinpath(@__DIR__, "landuse_wesely_conus_0p05.nc")

# Extracted window: CONUS plus a generous buffer. Bounds are cell *edges* and
# are exact multiples of 0.05°, so the window aligns with MCD12C1 cell
# boundaries and no source cell is split.
const WIN_LON_MIN = -130.0
const WIN_LON_MAX =  -60.0
const WIN_LAT_MIN =   20.0
const WIN_LAT_MAX =   55.0

# MCD12C1 CMG grid metadata.
const MCD_LON_STEP = 0.05
const MCD_LAT_STEP = 0.05
const MCD_NX = 7200          # lon -180..180
const MCD_NY = 3600          # lat  90..-90 (decreasing j → decreasing lat)
const MCD_N_CLASSES = 17     # IGBP 0..16

# Wesely (1989) class indices (1-based, matching the `wesleyLandUse` enum).
const WESELY_CLASSES = [
    "urban", "agricultural", "range", "deciduous", "coniferous",
    "mixedforest", "water", "barren", "wetland", "rangeag", "rockyshrubs",
]
const N_WESELY = length(WESELY_CLASSES)

# IGBP (0-based) → Wesely (1-based).
const IGBP_TO_WESELY = Dict{Int, Int}(
    0  => 7,   # Water bodies                              → wesleyWater
    1  => 5,   # Evergreen needleleaf forest               → wesleyConiferous
    2  => 4,   # Evergreen broadleaf forest                → wesleyDeciduous
    3  => 5,   # Deciduous needleleaf forest               → wesleyConiferous
    4  => 4,   # Deciduous broadleaf forest                → wesleyDeciduous
    5  => 6,   # Mixed forest                              → wesleyMixedForest
    6  => 11,  # Closed shrubland                          → wesleyRockyShrubs
    7  => 11,  # Open shrubland                            → wesleyRockyShrubs
    8  => 3,   # Woody savanna                             → wesleyRange
    9  => 3,   # Savanna                                   → wesleyRange
    10 => 3,   # Grassland                                 → wesleyRange
    11 => 9,   # Permanent wetland                         → wesleyWetland
    12 => 2,   # Cropland                                  → wesleyAgricultural
    13 => 1,   # Urban and built-up                        → wesleyUrban
    14 => 10,  # Cropland / natural-vegetation mosaic      → wesleyRangeAg
    15 => 8,   # Permanent snow and ice                    → wesleyBarren
    16 => 8,   # Barren or sparsely vegetated              → wesleyBarren
)

# ─────────────────────────────────────────────────────────────────────
# Download
# ─────────────────────────────────────────────────────────────────────

function ensure_mcd12c1()
    if isfile(MCD12C1_PATH) && filesize(MCD12C1_PATH) > 10_000_000
        @info "MCD12C1 already at $MCD12C1_PATH (skipping download)"
        return MCD12C1_PATH
    end
    @info "Downloading MCD12C1 from\n    $MCD12C1_URL\nto\n    $MCD12C1_PATH"
    @info "Using NASA Earthdata credentials from ~/.netrc."
    @info "(File is ~1.2 GB; expect a few minutes on typical broadband.)"
    try
        Downloads.download(MCD12C1_URL, MCD12C1_PATH; verbose = false)
    catch e
        @error """
        Download failed. Common causes:
          1. ~/.netrc missing the urs.earthdata.nasa.gov entry, or
             credentials wrong. Register at https://urs.earthdata.nasa.gov/.
          2. URL outdated — NASA reprocesses periodically. Query CMR for
             the current granule:
                curl -s "https://cmr.earthdata.nasa.gov/search/granules.json?short_name=MCD12C1&version=061&temporal=$(MCD12C1_YEAR)-01-01,$(MCD12C1_YEAR)-12-31&page_size=1"
             and update MCD12C1_GRANULE in this script.
          3. NASA host transiently down — retry.
        """
        rethrow(e)
    end
    @info "Downloaded $(round(Int, filesize(MCD12C1_PATH) / 1024 / 1024)) MB."
    return MCD12C1_PATH
end

# ─────────────────────────────────────────────────────────────────────
# Read MCD12C1 (HDF4)
# ─────────────────────────────────────────────────────────────────────

# Returns a (MCD_NX, MCD_NY, MCD_N_CLASSES) UInt8 array of percent values 0..100.
function read_mcd12c1(hdf_path)
    # HDF4-EOS subdataset path syntax.
    subdataset = "HDF4_EOS:EOS_GRID:\"$hdf_path\":MOD12C1:Land_Cover_Type_1_Percent"
    @info "Opening MCD12C1 subdataset:\n    $subdataset"
    local ds
    try
        ds = ArchGDAL.read(subdataset)
    catch e
        @error """
        Could not open MCD12C1 HDF4 subdataset. The GDAL_jll bundled with
        ArchGDAL may not have HDF4 driver support.
        """
        rethrow(e)
    end
    nx = ArchGDAL.width(ds)
    ny = ArchGDAL.height(ds)
    nb = ArchGDAL.nraster(ds)
    @info "MCD12C1 raster: $nx × $ny × $nb"
    @assert nx == MCD_NX  "expected $MCD_NX lon cells, got $nx"
    @assert ny == MCD_NY  "expected $MCD_NY lat cells, got $ny"
    @assert nb == MCD_N_CLASSES  "expected $MCD_N_CLASSES IGBP bands, got $nb"

    igbp = Array{UInt8, 3}(undef, nx, ny, nb)
    for b in 1:nb
        band = ArchGDAL.getband(ds, b)
        igbp[:, :, b] = ArchGDAL.read(band)
    end
    return igbp
end

# ─────────────────────────────────────────────────────────────────────
# Extract window + re-class IGBP → Wesely
# ─────────────────────────────────────────────────────────────────────

function extract_window(igbp)
    # Source cell i spans [-180 + 0.05(i-1), -180 + 0.05i]; cell j spans
    # [90 - 0.05j, 90 - 0.05(j-1)], i.e. lat decreases with j.
    i_lo = Int(round((WIN_LON_MIN + 180.0) / MCD_LON_STEP)) + 1
    i_hi = Int(round((WIN_LON_MAX + 180.0) / MCD_LON_STEP))
    j_lo = Int(round((90.0 - WIN_LAT_MAX) / MCD_LAT_STEP)) + 1
    j_hi = Int(round((90.0 - WIN_LAT_MIN) / MCD_LAT_STEP))
    @assert 1 <= i_lo <= i_hi <= MCD_NX "lon window out of range: $i_lo..$i_hi"
    @assert 1 <= j_lo <= j_hi <= MCD_NY "lat window out of range: $j_lo..$j_hi"

    n_lon_out = i_hi - i_lo + 1
    n_lat_out = j_hi - j_lo + 1
    @info "Extracting window at native $(MCD_LON_STEP)°: $n_lon_out × $n_lat_out"

    # Cell centres. Latitude is emitted increasing, so the source rows are
    # reversed on write; `landuse_frac_at` assumes an increasing axis.
    lon_centers = [WIN_LON_MIN + MCD_LON_STEP * (k - 0.5) for k in 1:n_lon_out]
    lat_centers = [WIN_LAT_MIN + MCD_LAT_STEP * (k - 0.5) for k in 1:n_lat_out]

    # Percent per Wesely class, 0..100. Several IGBP classes map to one Wesely
    # class, so accumulate; the total per cell stays ≤ 100.
    output = zeros(UInt8, n_lon_out, n_lat_out, N_WESELY)
    acc = zeros(UInt16, n_lon_out, n_lat_out, N_WESELY)
    for igbp_class in 0:MCD_N_CLASSES-1
        wes_idx = get(IGBP_TO_WESELY, igbp_class, 0)
        wes_idx == 0 && continue  # IGBP fill / unmapped
        b = igbp_class + 1
        @inbounds for j_out in 1:n_lat_out
            # Reverse: output row 1 is the southernmost, source row j_hi.
            jj = j_hi - j_out + 1
            for i_out in 1:n_lon_out
                acc[i_out, j_out, wes_idx] += UInt16(igbp[i_lo + i_out - 1, jj, b])
            end
        end
    end
    over = 0
    @inbounds for k in eachindex(acc)
        v = acc[k]
        v > 100 && (over += 1; v = UInt16(100))
        output[k] = UInt8(v)
    end
    over > 0 && @warn "$over class totals exceeded 100% and were capped"
    return output, lon_centers, lat_centers
end

# ─────────────────────────────────────────────────────────────────────
# Write NetCDF
# ─────────────────────────────────────────────────────────────────────

function write_netcdf(output, lon_centers, lat_centers; source_desc)
    n_lon, n_lat, n_class = size(output)
    isfile(OUTPUT_NC) && rm(OUTPUT_NC)
    NCDataset(OUTPUT_NC, "c") do ds
        defDim(ds, "lon", n_lon)
        defDim(ds, "lat", n_lat)
        defDim(ds, "class", n_class)

        lon_v = defVar(ds, "lon", Float32, ("lon",))
        lon_v[:] = Float32.(lon_centers)
        lon_v.attrib["units"] = "degrees_east"
        lon_v.attrib["long_name"] = "longitude (cell center)"

        lat_v = defVar(ds, "lat", Float32, ("lat",))
        lat_v[:] = Float32.(lat_centers)
        lat_v.attrib["units"] = "degrees_north"
        lat_v.attrib["long_name"] = "latitude (cell center)"

        class_v = defVar(ds, "class", Int32, ("class",))
        class_v[:] = Int32.(1:n_class)
        class_v.attrib["long_name"] = "Wesely (1989) land-use class index"

        class_name_v = defVar(ds, "class_name", String, ("class",))
        class_name_v[:] = WESELY_CLASSES
        class_name_v.attrib["long_name"] = "Wesely (1989) class name"

        pct_v = defVar(ds, "percent", output, ("lon", "lat", "class");
            deflatelevel = 9, shuffle = true,
            chunksizes = [min(n_lon, 256), min(n_lat, 256), 1])
        pct_v.attrib["long_name"] = "Wesely land-use class area percentage"
        pct_v.attrib["units"] = "percent"
        pct_v.attrib["valid_min"] = UInt8(0)
        pct_v.attrib["valid_max"] = UInt8(100)

        ds.attrib["title"] = "Wesely (1989) land-use percentages, CONUS + buffer, native 0.05°"
        ds.attrib["source"] = source_desc
        ds.attrib["window_lon_min"] = WIN_LON_MIN
        ds.attrib["window_lon_max"] = WIN_LON_MAX
        ds.attrib["window_lat_min"] = WIN_LAT_MIN
        ds.attrib["window_lat_max"] = WIN_LAT_MAX
        ds.attrib["cell_size"] = MCD_LON_STEP
    end
    @info "Wrote $OUTPUT_NC ($(round(filesize(OUTPUT_NC) / 1024 / 1024, digits = 1)) MB)"
end

# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────

function main()
    hdf_path = ensure_mcd12c1()
    igbp = read_mcd12c1(hdf_path)
    output, lon_centers, lat_centers = extract_window(igbp)

    # Sanity: percentages should sum to ≈100 per cell. Warn if not.
    n_lon, n_lat, _ = size(output)
    bad_cells = 0
    min_sum, max_sum = typemax(Int), typemin(Int)
    for j in 1:n_lat, i in 1:n_lon
        s = sum(Int, @view output[i, j, :])
        min_sum = min(min_sum, s)
        max_sum = max(max_sum, s)
        (99 ≤ s ≤ 101) || (bad_cells += 1)
    end
    @info @sprintf("Percent-sum sanity: min=%d, max=%d, off-by->1%% cells: %d / %d",
        min_sum, max_sum, bad_cells, n_lon * n_lat)

    source_desc = "MODIS MCD12C1 v6.1 IGBP Land_Cover_Type_1_Percent, year " *
        "$MCD12C1_YEAR; re-classed to Wesely (1989) classes; native 0.05° " *
        "window $WIN_LON_MIN..$WIN_LON_MAX E, $WIN_LAT_MIN..$WIN_LAT_MAX N."
    write_netcdf(output, lon_centers, lat_centers; source_desc = source_desc)
end

main()
