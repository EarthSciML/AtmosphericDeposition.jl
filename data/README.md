# Land-use fractions data

This directory holds the bundled MODIS-derived Wesely (1989) land-use data
used by `DryDepositionGasFractional`, plus the offline script that produces
it.

The bundled file is stored at the native MCD12C1 resolution and is **not**
tied to any simulation grid. `src/landuse_fractions.jl` conservatively
area-averages it onto whatever grid the domain specifies, on first use.

## Files

- `landuse_wesely_conus_0p05.nc` — the bundled NetCDF read at runtime by
  `AtmosphericDeposition.landuse_frac_at`. Native 0.05° source window:
  - longitude: -130° to -60° (1400 cells)
  - latitude: 20° to 55° (700 cells)
  - class dim: 11 (Wesely 1989 land-use classes; see below)

  Bounds are cell *edges* on exact 0.05° multiples. A domain extending past
  them raises rather than silently clamping, so widen the window here if you
  need a larger region.
- `preprocess_landuse.jl` — full pipeline. Downloads MCD12C1 from NASA
  Earthdata, re-classes IGBP → Wesely, extracts the window, writes the
  NetCDF. **Dev-only — run once to regenerate the bundled file**, never as
  part of a simulation. Requires NASA Earthdata credentials in `~/.netrc`.
- `Project.toml` — environment for the preprocessing script. Heavier
  deps (`ArchGDAL`, `Downloads`) live here, not in the main package.

## Wesely (1989) land-use classes

| index | name              |
|-------|-------------------|
| 1     | urban             |
| 2     | agricultural      |
| 3     | range             |
| 4     | deciduous forest  |
| 5     | coniferous forest |
| 6     | mixed forest      |
| 7     | water             |
| 8     | barren            |
| 9     | wetland           |
| 10    | range/agricultural mixed |
| 11    | rocky shrubs      |

## Regenerating the bundled file

Set up NASA Earthdata credentials at https://urs.earthdata.nasa.gov/
and add to `~/.netrc`:

```
machine urs.earthdata.nasa.gov
  login YOUR_USERNAME
  password YOUR_PASSWORD
```

Then run:

```bash
julia --project=data -e 'using Pkg; Pkg.instantiate()'
julia --project=data data/preprocess_landuse.jl
```

This downloads MCD12C1 (1.24 GB HDF4), re-classes IGBP → Wesely, extracts
the 0.05° window, and overwrites `landuse_wesely_conus_0p05.nc` (~1.3 MB).

Re-running against granule `MCD12C1.A2016001.061.2022168010533` reproduces
the committed file byte for byte (SHA-256
`e3a9730ababae391b93830277f6d455db39eb576dcb9cb6275500d0c8505a51f`, checked
2026-10-22).

## NetCDF schema

```
dimensions:
  lon = 1400
  lat = 700
  class = 11

variables:
  float lon(lon)                 ; units = "degrees_east"
  float lat(lat)                 ; units = "degrees_north"
  int class(class)               ; long_name = "Wesely (1989) land-use class index"
  string class_name(class)       ; long_name = "Wesely (1989) class name"
  ubyte percent(lon, lat, class) ; units = "percent" ; valid_range = 0, 100
```

Percentages are integers summing to exactly 100 per cell, so the IGBP →
Wesely re-class is lossless and stores in one byte per class.
