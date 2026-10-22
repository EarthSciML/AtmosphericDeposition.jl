# Dry Deposition

## Model

`DryDepositionGasFractional` calculates dry deposition velocities and first-order
loss rates for gas-phase species. It uses fractional (mosaic) land use: rather
than one land-use class per grid cell, the deposition velocity is the
area-weighted sum over the 11 Wesely classes,

```
Vd = Σᵢ fᵢ / (Ra + Rb + Rcᵢ)
```

where the aerodynamic resistance `Ra` and quasi-laminar sublayer resistance `Rb`
follow Seinfeld and Pandis (2006) eqs. 19.13/19.14 and 19.17 and are the same for
every class, and the surface resistance `Rcᵢ` of class `i` follows Wesely (1989).

Deposition is applied only in the lowest model level; the rate is
`k = Vd · g · ρA / ΔP`.

## Building the model

```julia
using AtmosphericDeposition
using ModelingToolkit
using ModelingToolkit: t

model = DryDepositionGasFractional()
sys = mtkcompile(model)
```

The system is purely algebraic — it produces deposition velocities and rates, not
concentrations. After `mtkcompile` there are no unknowns and 264 observed
quantities: a velocity `v_X(t)` [m/s] and a rate `k_X(t)` [1/s] for each of the
132 supported species.

```julia
observed(sys)   # v_ACET(t), v_ACTA(t), ..., k_RCOOH(t)
```

## Parameters

```julia
parameters(sys)
```

Meteorology and surface state:

| parameter | description |
|-----------|-------------|
| `Ts`      | surface air temperature [K] |
| `z`       | height to the mid-point of level 1 [m] |
| `del_P`   | pressure thickness of level 1 [Pa] |
| `z₀`      | roughness length [m] |
| `u_star`  | friction velocity [m/s] |
| `L`       | Monin-Obukhov length [m] |
| `ρA`      | air density [kg/m³] |
| `G`       | solar irradiation [W/m²] |
| `θ`       | slope of the local terrain [radians] |
| `lev`     | level index; deposition is applied only where `lev == 1` |
| `season`  | Wesely (1989) season index, 1–5 |

Land-use area fractions, one per Wesely (1989) class, expected to sum to ≈ 1:

`f_urban`, `f_agricultural`, `f_range`, `f_deciduous`, `f_coniferous`,
`f_mixedforest`, `f_water`, `f_barren`, `f_wetland`, `f_rangeag`,
`f_rockyshrubs`.

By default `f_mixedforest = 1.0` and the rest are zero.

The remaining parameters (`g`, `κ`, `Rc_unit`, `unit_m`, …) are physical
constants and unit-handling scalars; they are not intended to be varied.

## Coupling to chemistry

On its own the model only reports `v_X` and `k_X`. To see concentrations change,
couple it to a chemical mechanism — the `k_X` are then applied as loss terms to
the matching species:

```julia
using AtmosphericDeposition, GasChem, EarthSciMLBase, ModelingToolkit

model = couple(SuperFast(), DryDepositionGasFractional())
sys = convert(System, model)
```

`GEOSChemGasPhase` and `Pollu` are wired the same way; see
`ext/GasChemExt.jl` for the species maps.

## Coupling to GEOS-FP meteorology

Coupling with `EarthSciData.GEOSFP` binds the surface meteorology, drives
`season` from the simulation date, and looks the 11 land-use fractions up from
the bundled MODIS-derived CONUS dataset at each cell's longitude and latitude:

```julia
using AtmosphericDeposition, EarthSciData, EarthSciMLBase, ModelingToolkit, Dates

domain = DomainInfo(
    DateTime(2022, 5, 1), DateTime(2022, 5, 2);
    latrange = deg2rad(24.0f0):deg2rad(2):deg2rad(50.0f0),
    lonrange = deg2rad(-125.0f0):deg2rad(2.5):deg2rad(-66.0f0),
    levrange = 1:10,
)

model = couple(GEOSFP("4x5", domain), DryDepositionGasFractional(domain))
sys = convert(System, model)
```

Passing `domain` to `DryDepositionGasFractional` declares the simulation grid, and
the land-use source is conservatively area-averaged onto it. Without it the
system still builds, but the first land-use lookup raises.

Note the bundled land-use dataset covers CONUS only (-130…-60 °E, 20…55 °N); a
domain extending outside it raises rather than being clamped to the boundary.
`season` assumes Northern-Hemisphere mid-latitudes.
