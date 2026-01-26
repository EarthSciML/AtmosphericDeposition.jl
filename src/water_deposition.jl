export DryDepositionWater

# ============================================================================
# Two-Film Model for Water Surfaces (Seinfeld & Pandis 2006, Section 19.5.1)
# Equations 19.31-19.49
# ============================================================================

# Unit conversion constant: the kL formulas give results in cm/h
# To convert to m/s: divide by (100 cm/m * 3600 s/h) = 360000
@constants kL_unit_convert = 360000 [description = "Conversion factor from cm/h to m/s (100*3600)"]
@constants unit_kL = 1 [unit = u"m/s", description = "Unit for mass transfer coefficient"]
@constants unit_rc = 1 [unit = u"s/m", description = "Unit for surface resistance"]
@constants T_K_unit = 1 [unit = u"K^-1", description = "Unit to convert temperature"]

"""
    kG_gas_transfer(u10)

Gas-phase mass transfer coefficient [m/s] for air-water interface.
From Hicks and Liss (1976), as cited in Seinfeld & Pandis Section 19.5.1.

The gas-phase transfer coefficient represents the resistance to gas transfer
in the atmospheric boundary layer above the water surface.

# Arguments
- `u10`: Wind speed at 10m (dimensionless, representing m/s value)

# Returns
- Gas-phase mass transfer coefficient [m/s]
"""
function kG_gas_transfer(u10)
    return 0.0013 * u10 * unit_kL  # Approximately 0.13% of wind speed, returns [m/s]
end

"""
    kL_liss_merlivat(u10, Sc_ratio)

Liquid-phase mass transfer coefficient [m/s] using Liss & Merlivat (1986).
From Seinfeld & Pandis Eq. 19.43.

This parameterization divides the wind speed range into three regimes:
- Smooth surface (u10 <= 3.6 m/s): kL proportional to u10
- Rough surface (3.6 < u10 <= 13 m/s): Enhanced transfer due to waves
- Breaking waves (u10 > 13 m/s): Further enhanced transfer due to bubble injection

# Arguments
- `u10`: Wind speed at 10m (dimensionless, representing m/s value)
- `Sc_ratio`: Schmidt number ratio = Sc_CO2(293K) / Sc_species(T) (dimensionless)

# Returns
- Liquid-phase mass transfer coefficient [m/s]

# Note
Original equations give kL in cm/h, converted here to m/s by dividing by 360000.
"""
function kL_liss_merlivat(u10, Sc_ratio)
    # Smooth regime: u10 <= 3.6 m/s
    kL_smooth = (0.17 * u10 * Sc_ratio^(2 / 3)) / kL_unit_convert * unit_kL

    # Rough regime: 3.6 < u10 <= 13 m/s
    kL_rough = (0.612 * Sc_ratio^(2 / 3) + (2.85 * u10 - 10.26) * Sc_ratio^0.5) / kL_unit_convert * unit_kL

    # Breaking wave regime: u10 > 13 m/s
    kL_break = (0.612 * Sc_ratio^(2 / 3) + (5.9 * u10 - 49.9) * Sc_ratio^0.5) / kL_unit_convert * unit_kL

    # Use ifelse for symbolic compatibility
    kL = ifelse(u10 <= 3.6, kL_smooth,
        ifelse(u10 <= 13.0, kL_rough, kL_break))

    return kL
end

"""
    kL_wanninkhof(u10, Sc_ratio)

Liquid-phase mass transfer coefficient [m/s] using Wanninkhof (1992).
From Seinfeld & Pandis Eq. 19.46.

This is a simpler parameterization that uses a quadratic dependence on wind speed,
which is generally preferred for global applications.

# Arguments
- `u10`: Wind speed at 10m (dimensionless, representing m/s value)
- `Sc_ratio`: Schmidt number ratio = Sc_CO2(293K) / Sc_species(T) (dimensionless)

# Returns
- Liquid-phase mass transfer coefficient [m/s]

# Note
Original equation gives kL in cm/h, converted here to m/s.
"""
function kL_wanninkhof(u10, Sc_ratio)
    return (0.31 * u10^2 * Sc_ratio^0.5) / kL_unit_convert * unit_kL
end

"""
    schmidt_CO2_seawater(T_C)

Schmidt number for CO2 in seawater at temperature T [Celsius].
From Seinfeld & Pandis Eq. 19.45.

The Schmidt number is the ratio of kinematic viscosity to molecular diffusivity.
For CO2 in seawater, it varies strongly with temperature.

# Arguments
- `T_C`: Temperature [Celsius] (dimensionless number)

# Returns
- Schmidt number for CO2 in seawater (dimensionless)
"""
function schmidt_CO2_seawater(T_C)
    return 2073.1 - 125.62 * T_C + 3.6276 * T_C^2 - 0.043219 * T_C^3
end

"""
    rc_water_wanninkhof(H_star, u10, T)

Surface resistance for gas deposition to water [s/m] using Wanninkhof (1992).
From Seinfeld & Pandis Eq. 19.47-19.49.

The two-film model treats gas exchange across the air-water interface as
transport through two stagnant films: one on the air side (characterized by kG)
and one on the water side (characterized by kL).

# Arguments
- `H_star`: Effective Henry's law constant [M/atm] (dimensionless number)
- `u10`: Wind speed at 10m [m/s] (dimensionless number)
- `T`: Temperature [K]

# Returns
- Surface resistance [s/m]

# Physical interpretation
- For highly soluble gases (H_star >> 1e3 M/atm), rc ≈ 1/kG (Eq. 19.49)
  Gas-phase resistance dominates; transfer is limited by atmospheric transport
- For slightly soluble gases (H_star << 1 M/atm), rc ≈ 1/(kL * H̃) (Eq. 19.48)
  Liquid-phase resistance dominates; transfer is limited by dissolution

# Note on units
The effective Henry's law constant H_star has units [M/atm] = [mol/L/atm].
To convert to the dimensionless Henry's law constant H̃ (Eq. 19.42):
H̃ = H* × R × T (with R in appropriate units)
"""
function rc_water_wanninkhof(H_star, u10, T)
    T_C = T * T_K_unit - 273.15  # Convert to Celsius (dimensionless)

    # Gas-phase transfer coefficient [m/s]
    kG = kG_gas_transfer(u10)

    # Dimensionless Henry's law constant (Eq. 19.42)
    # H̃ = H* × R × T where H* is in mol/L/atm
    # For dimensional consistency: H̃ = H* [M/atm] × R [L·atm/mol/K] × T [K]
    # Using R = 0.08206 L·atm/(mol·K) for this conversion
    H_tilde = H_star * 0.08206 * T * T_K_unit  # Dimensionless

    # Schmidt number for CO2 at reference temperature (293 K = 20°C)
    Sc_CO2_ref = 660.0  # Approximate value at 20°C

    # Schmidt number for CO2 at actual temperature
    Sc_CO2 = schmidt_CO2_seawater(T_C)

    # Schmidt number ratio (assuming species has similar Sc to CO2)
    # For other species, this would need to be adjusted based on molecular diffusivity
    Sc_ratio = Sc_CO2_ref / Sc_CO2

    # Liquid-phase mass transfer coefficient using Wanninkhof (1992) [m/s]
    kL = kL_wanninkhof(u10, Sc_ratio)

    # Total surface resistance (Eq. 19.47) [s/m]
    # rc = 1/kG + 1/(kL × H̃)
    # The second term represents the liquid-phase resistance scaled by the
    # dimensionless Henry's law constant
    rc = 1 / kG + 1 / (kL * H_tilde)

    return rc
end

"""
    rc_water_liss_merlivat(H_star, u10, T)

Surface resistance for gas deposition to water [s/m] using Liss-Merlivat (1986).
From Seinfeld & Pandis Eq. 19.47-19.49.

Same as rc_water_wanninkhof but uses the Liss-Merlivat parameterization for kL.
"""
function rc_water_liss_merlivat(H_star, u10, T)
    T_C = T * T_K_unit - 273.15

    kG = kG_gas_transfer(u10)
    H_tilde = H_star * 0.08206 * T * T_K_unit
    Sc_CO2_ref = 660.0
    Sc_CO2 = schmidt_CO2_seawater(T_C)
    Sc_ratio = Sc_CO2_ref / Sc_CO2

    kL = kL_liss_merlivat(u10, Sc_ratio)

    rc = 1 / kG + 1 / (kL * H_tilde)

    return rc
end

# Default values for water deposition model (used for documentation purposes)
# The actual constants are declared at the module level with @constants

struct DryDepositionWaterCoupler
    sys::Any
end

"""
    DryDepositionWater(; name=:DryDepositionWater, use_liss_merlivat=false)

Create a ModelingToolkit System for gas dry deposition to water surfaces
using the two-film model from Seinfeld & Pandis (2006) Section 19.5.1.

The two-film model treats gas exchange across the air-water interface as
transport through two stagnant films: one on the air side and one on the
water side. This is appropriate for deposition to oceans, lakes, and other
water bodies.

# Keyword Arguments
- `name`: Name of the system (default: :DryDepositionWater)
- `use_liss_merlivat`: If true, use Liss-Merlivat (1986) for kL; otherwise use Wanninkhof (1992) (default: false)

# Parameters
- `z`: Height of surface layer [m] (default: 50)
- `z₀`: Roughness length [m] (default: 0.0002 for water)
- `u_star`: Friction velocity [m/s] (default: 0.1)
- `L`: Monin-Obukhov length [m] (default: 0)
- `ρA`: Air density [kg/m³] (default: 1.2)
- `T`: Temperature [K] (default: 293)
- `u10`: Wind speed at 10m [m/s] (default: 5)
- `H_star`: Effective Henry's law constant [M/atm] (default: 1e3 for moderately soluble gas)
- `lev`: Atmospheric level (default: 1)

# Variables
- `v_dep`: Deposition velocity [m/s]
- `k_dep`: Deposition rate [1/s]

# Example

```julia
using AtmosphericDeposition
using ModelingToolkit

# Use Wanninkhof (1992) parameterization (default)
model = DryDepositionWater()

# Use Liss-Merlivat (1986) parameterization
model_lm = DryDepositionWater(use_liss_merlivat=true)

sys = structural_simplify(model)
```

# References
- Seinfeld, J.H. and Pandis, S.N. (2006) Atmospheric Chemistry and Physics:
  From Air Pollution to Climate Change. 2nd Edition, John Wiley & Sons, New York.
  Chapter 19, Section 19.5.1, Equations 19.31-19.49.
- Liss, P.S. and Merlivat, L. (1986) Air-sea gas exchange rates: Introduction
  and synthesis. In The Role of Air-Sea Exchange in Geochemical Cycling.
- Wanninkhof, R. (1992) Relationship between wind speed and gas exchange over
  the ocean. Journal of Geophysical Research, 97(C5), 7373-7382.
"""
function DryDepositionWater(; name = :DryDepositionWater, use_liss_merlivat = false)
    params = @parameters begin
        z = 50, [unit = u"m", description = "Height of surface layer"]
        z₀ = 0.0002, [unit = u"m", description = "Roughness length for water"]
        u_star = 0.1, [unit = u"m/s", description = "Friction velocity"]
        L = 0, [unit = u"m", description = "Monin-Obukhov length"]
        ρA = 1.2, [unit = u"kg*m^-3", description = "Air density"]
        T = 293, [unit = u"K", description = "Temperature"]
        u10 = 5, [description = "Wind speed at 10m height (m/s value, dimensionless parameter)"]
        H_star = 1e3, [description = "Effective Henry's law constant (M/atm value, dimensionless parameter)"]
        lev = 1, [description = "Atmospheric level"]
    end

    vars = @variables begin
        v_dep(t), [unit = u"m/s", description = "Deposition velocity to water"]
        k_dep(t), [unit = u"1/s", description = "Deposition rate to water"]
    end

    # Calculate aerodynamic resistance [s/m]
    Ra = ra(z, z₀, u_star, L)

    # Calculate quasi-laminar sublayer resistance for gas [s/m]
    μ = mu(T)
    Dg = dH2O(T)  # Use water vapor diffusivity as proxy
    Sc = sc(μ, ρA, Dg)
    Rb = RbGas(Sc, u_star)

    # Calculate surface resistance using two-film model [s/m]
    # Selection of kL parameterization is done at model creation time
    if use_liss_merlivat
        Rc = rc_water_liss_merlivat(H_star, u10, T)
    else
        Rc = rc_water_wanninkhof(H_star, u10, T)
    end

    # Deposition velocity (only apply at surface level)
    i = ifelse(lev == 1, 1, 0)

    eqs = [
        v_dep ~ i / (Ra + Rb + Rc),
        k_dep ~ v_dep / z
    ]

    System(
        eqs,
        t,
        vars,
        [params; [kL_unit_convert, unit_kL, unit_rc, T_K_unit, unit_T, unit_convert_mu, T_unitless, unit_dH2O, unit_m, κ]];
        name = name,
        metadata = Dict(CoupleType => DryDepositionWaterCoupler)
    )
end
