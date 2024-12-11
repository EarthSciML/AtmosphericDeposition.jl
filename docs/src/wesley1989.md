# Wesley 1989 Dry Deposition Surface Resistance
## Introduction
This is the Wesely 1989 algorithm for surface resistance to dry deposition.

Citation for the original article, followed by citation for an article with some corrections which have been 
incorporated here:

1. M. L. Wesely, Parameterization of surface resistances to gaseous dry deposition in regional-scale numerical models, 
Atmos. Environ. 23, 1293–1304 (1989), http://dx.doi.org/10.1016/0004-6981(89)90153-4.

2. J. Walmsley, and M. L. Wesely, Modification of coded parametrizations of surface resistances to gaseous dry deposition, 
Atmos. Environ. 30(7), 1181–1188 (1996), http://dx.doi.org/10.1016/1352-2310(95)00403-3.

The abstract of the original article:

Methods for estimating the dry deposition velocities of atmospheric gases in the U.S. and surrounding areas have been improved and incorporated into a revised computer code module for use in numerical models of atmospheric transport and deposition of pollutants over regional scales. The key improvement is the computation of bulk surface resistances along three distinct pathways of mass transfer to sites of deposition at the upper portions of vegetative canopies or structures, the lower portions, and the ground (or water surface). This approach replaces the previous technique of providing simple look-up tables of bulk surface resistances. With the surface resistances divided explicitly into distinct pathways, the bulk surface resistances for a large number of gases in addition to those usually addressed in acid deposition models (SO₂, O₃, NOₓ, and HNO₃) can be computed, if estimates of the effective Henry’s Law constants and appropriate measures of the chemical reactivity of the various substances are known. This has been accomplished successfully for H₂O₂, HCHO, CH₃O₂H (to represent organic peroxides), CH₃C(O)O₂H, HCOOH (to represent organic acids), NH₃, CH₃C(O)O₂NO₂, and HNO₂. Other factors considered include surface temperature, stomatal response to environmental parameters, the wetting of surfaces by dew and rain, and the covering of surfaces by snow. Surface emission of gases and variations of uptake characteristics by individual plant species within the land use types are not considered explicitly.

## Main functions
Function ```WesleySurfaceResistance``` is used to calculate surface resistance to dry depostion [s/m] based on Wesely (1989) equation 2.
The inputs of the function are information on the chemical of interest ```gasData```, solar irradiation ```G``` [W/m²], the surface air temperature ```Ts``` [°C], the slope of the local terrain ```Θ``` [radians], the season index ```iSeason```, the land use index ```iLandUse```, whether there is currently rain or dew ```rain``` or ```dew```, and whether the chemical of interest is either SO₂ ```isSO2``` or O₃ ```isO3```.
Here's an example:

```julia @example 1
gasData::GasData = AtmosphericDeposition.So2Data
G, Ts, θ = [1.0, 20.0, 1.0] # [W/m², °C, radians]
iSeason, iLandUse = [1, 1] # the season index and land use index need to be integers
rain::Bool, dew::Bool, isSO2::Bool, isO3::Bool = [true, true, true, false]

WesleySurfaceResistance(gasData, G, Ts, θ, iSeason, iLandUse, rain, dew, isSO2, isO3) # [s/m]
```
This will return you the value of surface resistance to dry deposition of SO₂ with given solar irradiation, temperature and local terrain during midsummer with lush vegetation (season) in areas of evergreen needleleaf trees (land).

## Default parameters
a. For season index ```iSeason```, there're five seasonal categories

    1. Midsummer with lush vegetation
    2. Autumn with cropland not harvested
    3. Late autumn after frost, no snow
    4. Winter, snow on ground
    5. Transitional 

b. For land use index ```iLandUse```, there're eleven land use categories

    1. Urban land 
    2. agricultural land
    3. range land
    4. deciduous forest
    5. coniferous forest 
    6. mixed forest including wetland 
    7. water, both salt and fresh 
    8. barren land, mostly desert
    9. nonforested wetland
    10. mixed agricultural and range land
    11. rocky open areas with low-growing shrubs

c. For gasData ```gasData```, we include the following species: 

| Variable name | Species |
| -----------| --- |
|So2Data | SO2|
|O3Data | O3|
|No2Data | NO2|
|NoData | NO|
|Hno3Data | HNO3|
|H2o2Data | H2O2|
|AldData  | Acetaldehyde (aldehyde class)|
|HchoData  | Formaldehyde|
|OpData | Methyl hydroperoxide (organic peroxide class)|
|PaaData | Peroxyacetyl nitrate|
|OraData | Formic acid (organic acid class)|
|Nh3Data | NH3|
|PanData | Peroxyacetyl nitrate|
|Hno2Data | Nitrous acid|

Wesely (1989) suggests that, in general, the sum of NO and NO₂ should be considered rather than NO₂ alone because rapid in-air chemical reactions can cause a significant change of NO and NO₂ vertical fluxes between the surface and the point at which deposition velocities are applied, but the sum of NO and NO₂ fluxes should be practically unchanged.