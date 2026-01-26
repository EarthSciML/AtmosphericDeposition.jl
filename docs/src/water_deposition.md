# Dry Deposition to Water Surfaces

## Overview

This module implements the **two-film model** for gas dry deposition to water surfaces, based on Seinfeld & Pandis (2006) Chapter 19, Section 19.5.1, Equations 19.31-19.49.

The two-film model is appropriate for calculating gas exchange between the atmosphere and water bodies such as oceans, lakes, and rivers. Unlike deposition to land surfaces (which uses surface resistance parameterizations like Wesely 1989), deposition to water is controlled by diffusive transport through thin boundary layers on both sides of the air-water interface.

## Physical Background

### The Two-Film Concept

Gas exchange across the air-water interface is modeled as transport through two stagnant films:

1. **Gas-phase film**: A thin layer of air just above the water surface where molecular diffusion dominates
2. **Liquid-phase film**: A thin layer of water just below the surface where molecular diffusion dominates

The total resistance to gas transfer is the sum of resistances in these two films.

### Governing Equations

The **total surface resistance** (Eq. 19.47) is:

```math
r_c = \frac{1}{k_G} + \frac{1}{k_L \cdot \tilde{H}}
```

where:
- ``k_G`` is the gas-phase transfer coefficient [m/s]
- ``k_L`` is the liquid-phase transfer coefficient [m/s]
- ``\tilde{H}`` is the dimensionless Henry's law constant

### Henry's Law Constant

The **dimensionless Henry's law constant** (Eq. 19.42) relates gas and liquid phase concentrations:

```math
\tilde{H} = H^* \cdot R \cdot T
```

where:
- ``H^*`` is the effective Henry's law constant [M/atm]
- ``R`` is the gas constant
- ``T`` is temperature [K]

### Limiting Cases

The relative importance of gas-phase and liquid-phase resistances depends on the gas solubility:

**Highly soluble gases** (``H^* >> 10^3`` M/atm):
```math
r_c \approx \frac{1}{k_G}
```
Gas-phase resistance dominates. Examples: HNO₃, NH₃, SO₂

**Slightly soluble gases** (``H^* << 1`` M/atm):
```math
r_c \approx \frac{1}{k_L \cdot \tilde{H}}
```
Liquid-phase resistance dominates. Examples: O₃, NO, CO

## Parameterizations

### Gas-Phase Transfer Coefficient

The gas-phase transfer coefficient (Hicks & Liss, 1976) is proportional to wind speed:

```math
k_G = 0.0013 \cdot u_{10}
```

where ``u_{10}`` is the wind speed at 10 m height [m/s].

### Liquid-Phase Transfer Coefficient

Two parameterizations are available:

#### Liss-Merlivat (1986)

This parameterization divides wind speed into three regimes (Eq. 19.43):

- **Smooth surface** (``u_{10} \leq 3.6`` m/s):
  ```math
  k_L = 0.17 \cdot u_{10} \cdot S_A^{2/3}
  ```

- **Rough surface** (``3.6 < u_{10} \leq 13`` m/s):
  ```math
  k_L = 0.612 \cdot S_A^{2/3} + (2.85 \cdot u_{10} - 10.26) \cdot S_A^{0.5}
  ```

- **Breaking waves** (``u_{10} > 13`` m/s):
  ```math
  k_L = 0.612 \cdot S_A^{2/3} + (5.9 \cdot u_{10} - 49.9) \cdot S_A^{0.5}
  ```

where ``S_A = Sc_{CO2}(293K) / Sc_{species}(T)`` is the Schmidt number ratio.

#### Wanninkhof (1992)

A simpler quadratic parameterization (Eq. 19.46):

```math
k_L = 0.31 \cdot u_{10}^2 \cdot S_A^{0.5}
```

This is generally preferred for global applications.

### Schmidt Number for CO₂

The Schmidt number for CO₂ in seawater (Eq. 19.45):

```math
Sc_{CO2} = 2073.1 - 125.62T + 3.6276T^2 - 0.043219T^3
```

where ``T`` is temperature in Celsius.

## Usage

### Basic Usage

```julia
using AtmosphericDeposition
using ModelingToolkit

# Create the model with default parameters
model = DryDepositionWater()

# Simplify the system
sys = structural_simplify(model)
```

### Setting Parameters

```julia
using AtmosphericDeposition
using ModelingToolkit

model = DryDepositionWater()
sys = structural_simplify(model)

# Access parameters
@unpack u10, H_star, T, z = sys

# Example: SO₂ deposition (highly soluble)
# H* for SO₂ ≈ 10⁵ M/atm
params_SO2 = [H_star => 1e5, u10 => 8.0, T => 293.0, z => 50.0]

# Example: O₃ deposition (slightly soluble)
# H* for O₃ ≈ 0.01 M/atm
params_O3 = [H_star => 0.01, u10 => 8.0, T => 293.0, z => 50.0]
```

### Available Parameters

| Parameter | Description | Default | Units |
|-----------|-------------|---------|-------|
| `z` | Height of surface layer | 50 | m |
| `z₀` | Roughness length | 0.0002 | m |
| `u_star` | Friction velocity | 0.1 | m/s |
| `L` | Monin-Obukhov length | 0 | m |
| `ρA` | Air density | 1.2 | kg/m³ |
| `T` | Temperature | 293 | K |
| `u10` | Wind speed at 10m | 5 | m/s |
| `H_star` | Henry's law constant | 1000 | M/atm |
| `use_wanninkhof` | Use Wanninkhof parameterization | true | - |
| `lev` | Atmospheric level | 1 | - |

### Output Variables

| Variable | Description | Units |
|----------|-------------|-------|
| `v_dep` | Deposition velocity | m/s |
| `k_dep` | Deposition rate | 1/s |

## Henry's Law Constants for Common Gases

| Gas | H* [M/atm] | Solubility |
|-----|------------|------------|
| HNO₃ | 10¹⁴ | Very high |
| NH₃ | 2×10⁴ | High |
| SO₂ | 10⁵ | High |
| H₂O₂ | 10⁵ | High |
| HCHO | 6×10³ | Moderate |
| O₃ | 0.01 | Low |
| NO | 0.003 | Very low |
| CO | 0.001 | Very low |

## Coupling with Other Models

The `DryDepositionWater` model can be coupled with other EarthSciML models using the standard coupling interface:

```julia
using EarthSciMLBase

# Example: Couple with atmospheric chemistry model
# coupled_model = couple(chemistry_model, DryDepositionWater())
```

## References

- Seinfeld, J.H. and Pandis, S.N. (2006) Atmospheric Chemistry and Physics: From Air Pollution to Climate Change. 2nd Edition, John Wiley & Sons, New York.
- Hicks, B.B. and Liss, P.S. (1976) Transfer of SO₂ and other reactive gases across the air-sea interface. Tellus, 28, 348-354.
- Liss, P.S. and Merlivat, L. (1986) Air-sea gas exchange rates: Introduction and synthesis. In The Role of Air-Sea Exchange in Geochemical Cycling (P. Buat-Menard, ed.), pp. 113-127, D. Reidel.
- Wanninkhof, R. (1992) Relationship between wind speed and gas exchange over the ocean. Journal of Geophysical Research, 97(C5), 7373-7382.
