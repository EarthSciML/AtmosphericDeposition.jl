export WesleySurfaceResistance

const inf = 1.0e25

@enum wesleyLandUse begin # Land use categories:
    wesleyUrban = 1 # Urban land
    wesleyAgricultural # agricultural land
    wesleyRange # range land
    wesleyDeciduous # deciduous forest
    wesleyConiferous # coniferous forest
    wesleyMixedForest # mixed forest including wetland
    wesleyWater # water, both salt and fresh
    wesleyBarren # barren land, mostly desert
    wesleyWetland # nonforested wetland
    wesleyRangeAg # mixed agricultural and range land
    wesleyRockyShrubs # rocky open areas with low-growing shrubs
end

@enum wesleySeason begin # Season categories:
    wesleyMidsummer = 1 # Midsummer with lush vegetation
    wesleyAutumn    # Autumn with unharvested cropland
    wesleyLateAutumn # Late autumn after frost, no snow
    wesleyWinter     # Winter, snow on ground and subfreezing
    wesleyTransitional  # Transitional spring with partially green short annuals
end

# r_i represents the minimum bulk canopy resistances for water vapor.
const r_i = SA_F32[
    inf 60 120 70 130 100 inf inf 80 100 150
    inf inf inf inf 250 500 inf inf inf inf inf
    inf inf inf inf 250 500 inf inf inf inf inf
    inf inf inf inf 400 800 inf inf inf inf inf
    inf 120 240 140 250 190 inf inf 160 200 300
]

# r_lu signifies leaf cuticles in healthy vegetation and otherwise the outer surfaces in the upper canopy.
const r_lu = SA_F32[
    inf 2000 2000 2000 2000 2000 inf inf 2500 2000 4000
    inf 9000 9000 9000 4000 8000 inf inf 9000 9000 9000
    inf inf 9000 9000 4000 8000 inf inf 9000 9000 9000
    inf inf inf inf 6000 9000 inf inf 9000 9000 9000
    inf 4000 4000 4000 2000 3000 inf inf 4000 4000 8000
]

# r_ac signifies transfer that depends only on canopy height and density.
const r_ac = SA_F32[
    100.0 200 100 2000 2000 2000 0 0 300 150 200
    100 150 100 1500 2000 1700 0 0 200 120 140
    100 10 100 1000 2000 1500 0 0 100 50 120
    100 10 10 1000 2000 1500 0 0 50 10 50
    100 50 80 1200 2000 1500 0 0 200 60 120
]

# r_gs signifies uptake at the "ground" by soil, leaf litter, snow, water etc. 'S' and 'O' stand for SO2 and O3 respectively.
const r_gsS = SA_F32[
    400.0 150 350 500 500 100 0 1000 0 220 40
    400 200 350 500 500 100 0 1000 0 300 400
    400 150 350 500 500 200 0 1000 0 200 400
    100 100 100 100 100 100 0 1000 100 100 50
    500 150 350 500 500 200 0 1000 0 250 40
]

const r_gsO = SA_F32[
    300.0 150 200 200 200 300 2000 400 1000 180 200
    300 150 200 200 200 300 2000 400 800 180 200
    300 150 200 200 200 300 2000 400 1000 180 20
    600 3500 3500 3500 3500 3500 2000 400 3500 3500 3500
    300 150 200 200 200 300 2000 400 1000 180 200
]

# r_cl is meant to account for uptake pathways at the leaves, bark, etc. 'S' and 'O' stand for SO2 and O3 respectively.
const r_clS = SA_F32[
    inf 2000 2000 2000 2000 2000 inf inf 2500 2000 4000
    inf 9000 9000 9000 2000 4000 inf inf 9000 9000 9000
    inf inf 9000 9000 3000 6000 inf inf 9000 9000 9000
    inf inf inf 9000 200 400 inf inf 9000 inf 9000
    inf 4000 4000 4000 2000 3000 inf inf 4000 4000 8000
]

const r_clO = SA_F32[
    inf 1000 1000 1000 1000 1000 inf inf 1000 1000 1000
    inf 400 400 400 1000 600 inf inf 400 400 400
    inf 1000 400 400 1000 600 inf inf 800 600 600
    inf 1000 1000 400 1500 600 inf inf 800 1000 800
    inf 1000 500 500 1500 700 inf inf 600 800 800
]

# Holder for gas properties from Wesely (1989) Table 2.'
struct GasData
    Dh2oPerDx::AbstractFloat
    Hstar::AbstractFloat
    Fo::AbstractFloat
end

const So2Data = GasData(1.9, 1.0e5, 0)
const No2Data = GasData(1.6, 0.01, 0.1) # Wesely (1989) suggests that,
# in general, the sum of NO and NO2 should be considered rather
# than NO2 alone because rapid in-air chemical reactions can cause
# a significant change of NO and NO2 vertical fluxes between the
# surface and the point at which deposition velocities are applied,
# but the sum of NO and NO2 fluxes should be practically unchanged.
const NoData = GasData(1.3, 3.0e-3, 0) # Changed according to Walmsley (1996)
const Hno3Data = GasData(1.9, 1.0e14, 0)
const H2o2Data = GasData(1.4, 1.0e5, 1)
const AldData = GasData(1.6, 15, 0)     # Acetaldehyde (aldehyde class)
const HchoData = GasData(1.3, 6.0e3, 0)   # Formaldehyde
const OpData = GasData(1.6, 240, 0.1) # Methyl hydroperoxide (organic peroxide class)
const PaaData = GasData(2.0, 540, 0.1)  # Peroxyacetyl nitrate
const OraData = GasData(1.6, 4.0e6, 0)   # Formic acid (organic acid class)
const Nh3Data = GasData(0.97, 2.0e4, 0)  # Changed according to Walmsley (1996)
const PanData = GasData(2.6, 3.6, 0.1)  # Peroxyacetyl nitrate
const Hno2Data = GasData(1.6, 1.0e5, 0.1) # Nitrous acid
const ACETData = GasData(1.795685104, 100000.0, 1.0)
const ACTAData = GasData(1.825879681, 4100.0, 1.0)
const ALD2Data = GasData(1.563873924, 11.0, 1.0)
const AROMP4Data = GasData(1.943968688, 410000.0, 1.0)
const AROMP5Data = GasData(2.333533261, 2.0e6, 1.0)
const ATOOHData = GasData(2.236236775, 294.0, 1.0)
const BALDData = GasData(2.427046586, 38.0, 1.0)
const BENZPData = GasData(2.472252775, 2900.0, 1.0)
const Br2Data = GasData(2.978296144, 0.76, 0.0)
const BrClData = GasData(2.531491424, 0.97, 0.0)
const BrNO3Data = GasData(2.806635356, 1.0e20, 0.0)
const BZCO3HData = GasData(2.768903222, 24000.0, 1.0)
const BZPANData = GasData(3.188213391, 70.0, 1.0)
const CH2OData = GasData(1.291091904, 3000.0, 1.0)
const Cl2Data = GasData(1.983821576, 0.092, 0.0)
const ClNO2Data = GasData(2.126302433, 0.045, 0.0)
const ClNO3Data = GasData(2.325789543, 1.0e20, 0.0)
const ClOData = GasData(1.689943485, 0.7, 0.0)
const ClOOData = GasData(1.934953215, 1.0, 0.0)
const CSLData = GasData(2.450037177, 420.0, 1.0)
const EOHData = GasData(1.599147774, 190.0, 0.0)
const ETHLNData = GasData(2.414894654, 2.0e6, 1.0)
const ETHNData = GasData(2.437885999, 39000.0, 0.1)
const ETHPData = GasData(2.081716485, 650000.0, 0.1)
const ETNO3Data = GasData(2.248490219, 1.6, 0.1)
const ETPData = GasData(1.856330695, 294.0, 1.0)
const GLYCData = GasData(1.825879681, 41000.0, 1.0)
const GLYXData = GasData(1.794912135, 360000.0, 1.0)
const H2O2Data = GasData(1.374189565, 5.0e7, 1.0)
const HACData = GasData(2.027822692, 1.4e6, 1.0)
const HBrData = GasData(2.119242195, 7.1e15, 0.0)
const HC5AData = GasData(2.357553733, 7800.0, 0.0)
const HClData = GasData(1.422421336, 2.0e13, 0.0)
const HCOOHData = GasData(1.598453398, 8900.0, 1.0)
const HIData = GasData(2.664598268, 2.35e16, 0.0)
const HMHPData = GasData(1.885554366, 1.3e6, 1.0)
const HMMLData = GasData(2.380632525, 120000.0, 1.0)
const HNO3Data = GasData(1.870183545, 1.0e14, 0.0)
const HOBrData = GasData(2.319336638, 1300.0, 0.0)
const HOClData = GasData(1.706287613, 650.0, 0.0)
const HOIData = GasData(2.826147328, 15400.0, 0.0)
const HONITData = GasData(3.454607581, 2.0e6, 1.0)
const HPALD1Data = GasData(2.538935715, 40000.0, 0.0)
const HPALD2Data = GasData(2.538935715, 40000.0, 0.0)
const HPALD3Data = GasData(2.538935715, 40000.0, 0.0)
const HPALD4Data = GasData(2.538935715, 40000.0, 0.0)
const HPETHNLData = GasData(2.054743675, 41000.0, 1.0)
const I2Data = GasData(3.753403898, 2.7, 0.0)
const I2O2Data = GasData(3.983002729, 1.0e20, 0.0)
const I2O3Data = GasData(4.092975165, 1.0e20, 0.0)
const I2O4Data = GasData(4.200069126, 1.0e20, 0.0)
const IBrData = GasData(3.388907673, 24.3, 0.0)
const ICHEData = GasData(2.538935715, 8.0e7, 1.0)
const IClData = GasData(3.002889488, 111.0, 0.0)
const ICNData = GasData(2.838298642, 2.0e6, 1.0)
const ICPDHData = GasData(2.886969263, 1.0e8, 1.0)
const IDCData = GasData(2.333652194, 40000.0, 0.0)
const IDCHPData = GasData(2.867483999, 1.0e8, 1.0)
const IDHDPData = GasData(3.055299253, 1.0e8, 1.0)
const IDHPEData = GasData(2.886969263, 1.0e8, 1.0)
const IDNData = GasData(3.265875962, 1.0e8, 1.0)
const IEPOXAData = GasData(2.427275283, 8.0e7, 1.0)
const IEPOXBData = GasData(2.427275283, 8.0e7, 1.0)
const IEPOXDData = GasData(2.427275283, 8.0e7, 1.0)
const IHN1Data = GasData(2.857982893, 2.0e6, 1.0)
const IHN2Data = GasData(2.857982893, 2.0e6, 1.0)
const IHN3Data = GasData(2.857982893, 2.0e6, 1.0)
const IHN4Data = GasData(2.857982893, 2.0e6, 1.0)
const INPBData = GasData(3.009352286, 2.0e6, 1.0)
const INPDData = GasData(3.009352286, 2.0e6, 1.0)
const IONOData = GasData(3.098058022, 0.3, 0.0)
const IONO2Data = GasData(3.238224586, 1.0e20, 0.0)
const IPRNO3Data = GasData(2.415469232, 0.79, 0.1)
const ITCNData = GasData(3.291271958, 1.0e8, 1.0)
const ITHNData = GasData(3.308262103, 1.0e8, 1.0)
const LIMOData = GasData(2.750196241, 0.07, 0.0)
const LVOCData = GasData(2.925550478, 1.0e8, 1.0)
const MACRData = GasData(1.972597602, 6.5, 1.0)
const MACR1OOHData = GasData(2.380632525, 294.0, 1.0)
const MAPData = GasData(2.054743675, 840.0, 1.0)
const MCRDHData = GasData(2.404067025, 1.4e6, 1.0)
const MCRENOLData = GasData(2.186155589, 294.0, 1.0)
const MCRHNData = GasData(2.876953728, 2.0e6, 1.0)
const MCRHNBData = GasData(2.876953728, 2.0e6, 1.0)
const MCRHPData = GasData(2.582183808, 1.4e6, 1.0)
const MCTData = GasData(2.623555974, 420.0, 1.0)
const MENO3Data = GasData(2.068072755, 2.0, 0.1)
const MGLYData = GasData(2.000123225, 3700.0, 1.0)
const MOHData = GasData(1.333808586, 203.0, 1.0)
const MONITSData = GasData(3.456856361, 2.0e6, 1.0)
const MONITUData = GasData(3.456856361, 2.0e6, 1.0)
const MPANData = GasData(2.857497296, 1.72, 1.0)
const MTPAData = GasData(2.750196241, 0.049, 0.0)
const MTPOData = GasData(2.750196241, 0.049, 0.0)
const MVKData = GasData(1.972456898, 44.0, 1.0)
const MVKDHData = GasData(2.415699025, 1.4e6, 1.0)
const MVKHCData = GasData(2.380632525, 1.4e6, 1.0)
const MVKHCBData = GasData(2.380632525, 1.4e6, 1.0)
const MVKHPData = GasData(2.582183808, 1.4e6, 1.0)
const MVKNData = GasData(2.877050197, 2.0e6, 1.0)
const MVKPCData = GasData(2.560380085, 1.4e6, 1.0)
const N2O5Data = GasData(2.44867743, 1.0e14, 0.0)
const NO2Data = GasData(1.598106097, 0.01, 0.1)
const NPHENData = GasData(2.77880881, 2300.0, 1.0)
const NPRNO3Data = GasData(2.415469232, 1.1, 0.1)
const O3Data = GasData(1.632300488, 0.01, 1.0)
const PANData = GasData(2.592267569, 3.6, 1.0)
const PHENData = GasData(2.285585007, 2800.0, 1.0)
const PPData = GasData(2.261168272, 294.0, 1.0)
const PPNData = GasData(2.738262115, 3.6, 1.0)
const PROPNNData = GasData(2.570981223, 500000.0, 1.0)
const PRPNData = GasData(2.758760869, 294.0, 1.0)
const PYACData = GasData(2.211024169, 314000.0, 1.0)
const R4N2Data = GasData(2.571197117, 17000.0, 1.0)
const R4PData = GasData(2.236857245, 294.0, 1.0)
const RA3PData = GasData(2.055418934, 294.0, 1.0)
const RB3PData = GasData(2.055418934, 294.0, 1.0)
const RIPAData = GasData(2.560922022, 1.7e6, 1.0)
const RIPBData = GasData(2.560922022, 1.7e6, 1.0)
const RIPCData = GasData(2.560922022, 1.7e6, 1.0)
const RIPDData = GasData(2.560922022, 1.7e6, 1.0)
const RPData = GasData(2.236236775, 294.0, 1.0)
const SO2Data = GasData(1.885407166, 100000.0, 0.0)
const RCOOHData = GasData(2.027959554, 1520.0, 1.0)

# Obtain values from matrix using symbolic parameter iSeason and iLandUse
function obtain_value(iSeason, iLandUse, matrix)
    index = (iLandUse - 1) * 5 + iSeason
    interpolate_r_i = DataInterpolations.LinearInterpolation(vec(matrix), 1:55)
    return interpolate_r_i(index)
end

# Calculate bulk canopy stomatal resistance [s m-1] based on Wesely (1989) equation 3 when given the solar irradiation (G [W m-2]), the surface air temperature (Ts [°C]), the season index (iSeason), the land use index (iLandUse), and whether there is currently rain or dew.
function r_s(G, Ts, iSeason, iLandUse, rainOrDew::Bool)
    rs = 0.0
    rs = ifelse(
        (Ts >= 39.9),
        inf,
        obtain_value(iSeason, iLandUse, r_i) *
            (1 + (200.0 * 1.0 / (G + 1.0))^2) *
            (400.0 * 1.0 / (Ts * (40.0 - Ts)))
    )
    rs = ifelse(
        (Ts <= 0.1),
        inf,
        obtain_value(iSeason, iLandUse, r_i) *
            (1 + (200.0 * 1.0 / (G + 1.0))^2) *
            (400.0 * 1.0 / (Ts * (40.0 - Ts)))
    )
    # Adjust for dew and rain (from "Effects of dew and rain" section).
    if rainOrDew
        rs *= 3
    end
    return rs
end

# Calculate the resistance from the effects of mixing forced by buoyant convection when sunlight heats the ground or lower canopy and by penetration of wind into canopies on the sides of hills [s m-1] when given the solar irradiation (G [W m-2]) and the slope of the local terrain (θ [radians]). From Wesely (1989) equation 5.
function r_dc(G, θ)
    return 100.0 * (1.0 + 1000.0 / (G + 10.0)) / (1.0 + 1000.0 * θ)
end

#Calculate mesophyll resistance [s m-1] based on Wesely (1989) equation 6 when given the effective Henry's law coefficient (Hstar [M atm-1]) and the reactivity factor (fo [-]), both available in Wesely (1989) table 2.
function r_mx(Hstar::T, fo::T) where {T <: AbstractFloat}
    return 1.0 / (Hstar / 3000.0 + 100.0 * fo)
end

#Calculate combined minimum stomatal and mesophyll resistance [s m-1] based on Wesely (1989) equation 4 when given stomatal resistance (r_s [s m-1]), ratio of water to chemical-of-interest diffusivities (Dh2oPerDx [-]), and mesophyll resistance (r_mx [s m-1]).
function r_smx(r_s, Dh2oPerDx::T, r_mx::T) where {T <: AbstractFloat}
    return r_s * Dh2oPerDx + r_mx
end

# Calculate the resistance of the outer surfaces in the upper canopy (leaf cuticular resistance in healthy vegetation) based on Wesely (1989) equations 7 and 10-14 when given the effective Henry's law coefficient (Hstar [M atm-1]), the reactivity factor (fo [-]) (both available in Wesely (1989) table 2), the season index (iSeason), the land use index (iLandUse), whether there is currently rain or dew, and whether the chemical of interest is either SO2 or O3.
function r_lux(
        Hstar::T,
        fo::T,
        iSeason,
        iLandUse,
        rain::Bool,
        dew::Bool,
        isSO2::Bool,
        isO3::Bool
    ) where {T <: AbstractFloat}
    rlux = 0.0
    if dew && (iSeason != 4) # Dew doesn't have any effect in the winter
        if isSO2
            if iLandUse == 1
                rlux = 50.0 # equation 13 and a half
            else
                rlux = 100 # equation 10.
            end
        elseif isO3
            # equation 11
            rlux = 1.0 / (1.0 / 3000.0 + 1.0 / (3 * obtain_value(iSeason, iLandUse, r_lu)))
        else
            rluO = 1.0 / (1.0 / 3000.0 + 1.0 / (3 * obtain_value(iSeason, iLandUse, r_lu))) #equation 11
            rlux = 1.0 / (
                1.0 /
                    (3 * obtain_value(iSeason, iLandUse, r_lu) / (1.0e-5 * Hstar + fo)) +
                    1.0e-7 * Hstar +
                    fo / rluO
            ) # equation 14, modified to match Walmsley eq. 5g
        end
    elseif rain && (iSeason != 4)
        if isSO2
            if iLandUse == 1
                rlux = 50 #equation 13 and a half
            else
                # equation 12
                rlux = 1.0 /
                    (1.0 / 5000.0 + 1.0 / (3 * obtain_value(iSeason, iLandUse, r_lu)))
            end
        elseif isO3
            # equation 13
            rlux = 1.0 / (1.0 / 1000.0 + 1.0 / (3 * obtain_value(iSeason, iLandUse, r_lu)))
        else
            rluO = 1.0 / (1.0 / 1000.0 + 1.0 / (3 * obtain_value(iSeason, iLandUse, r_lu))) # equation 13
            rlux = 1.0 / (
                1.0 /
                    (3 * obtain_value(iSeason, iLandUse, r_lu) / (1.0e-5 * Hstar + fo)) +
                    1.0e-7 * Hstar +
                    fo / rluO
            ) # equation 14, modified to match Walmsley eq. 5g
        end
    else
        rlux = obtain_value(iSeason, iLandUse, r_lu) / (1.0e-5 * Hstar + fo)
    end
    return rlux
end

# Calculate the resistance of the exposed surfaces in the lower portions of structures (canopies, buildings) above the ground based on Wesely (1989) equation 8 when given the effective Henry's law coefficient (Hstar [M atm-1]), the reactivity factor (fo [-]) (both available in Wesely (1989) table 2), the season index (iSeason), and the land use index (iLandUse).
function r_clx(Hstar::T, fo::T, iSeason, iLandUse) where {T <: AbstractFloat}
    return 1.0 / (
        Hstar / (1.0e5 * obtain_value(iSeason, iLandUse, r_clS)) +
            fo / obtain_value(iSeason, iLandUse, r_clO)
    )
end

# Calculate the resistance to uptake at the 'ground' surface based on Wesely (1989) equation 9 when given the effective Henry's law coefficient (Hstar [M atm-1]), the reactivity factor (fo [-]) (both available in Wesely (1989) table 2), the season index (iSeason), and the land use index (iLandUse).
function r_gsx(Hstar::T, fo::T, iSeason, iLandUse) where {T <: AbstractFloat}
    return 1.0 / (
        Hstar / (1.0e5 * obtain_value(iSeason, iLandUse, r_gsS)) +
            fo / obtain_value(iSeason, iLandUse, r_gsO)
    )
end

# Calculates surface resistance to dry depostion [s m-1] based on Wesely (1989) equation 2 when given information on the chemical of interest (gasData), solar irradiation (G [W m-2]), the surface air temperature (Ts [°C]), the slope of the local terrain (Θ [radians]), the season index (iSeason), the land use index (iLandUse), whether there is currently rain or dew, and whether the chemical of interest is either SO2 (isSO2) or O3 (isO3).
# From Wesely (1989) regarding rain and dew inputs: "A direct computation of the surface wetness would be most desirable, e.g. by estimating the amount of free surface water accumulated and then evaporated. Alternatively, surface relative humidity might be a useful	index. After dewfall and rainfall events are completed, surface wetness	often disappears as a result of evaporation after approximately 2	hours of good atmospheric mixing, the period of time recommended earlier (Sheih et al., 1986)".
function WesleySurfaceResistance(
        gasData::GasData,
        G,
        Ts,
        θ,
        iSeason,
        iLandUse,
        rain::Bool,
        dew::Bool,
        isSO2::Bool,
        isO3::Bool
    )
    rs = r_s(G, Ts, iSeason, iLandUse, rain || dew)
    rmx = r_mx(gasData.Hstar, gasData.Fo)
    rsmx = r_smx(rs, gasData.Dh2oPerDx, rmx)
    rdc = r_dc(G, θ)
    rlux = r_lux(gasData.Hstar, gasData.Fo, iSeason, iLandUse, rain, dew, isSO2, isO3)
    rclx = 0.0
    rgsx = 0.0
    if isSO2
        rclx = obtain_value(iSeason, iLandUse, r_clS)
        rgsx = obtain_value(iSeason, iLandUse, r_gsS)
    elseif isO3
        rclx = obtain_value(iSeason, iLandUse, r_clO)
        rgsx = obtain_value(iSeason, iLandUse, r_gsO)
    else
        rclx = r_clx(gasData.Hstar, gasData.Fo, iSeason, iLandUse)
        rgsx = r_gsx(gasData.Hstar, gasData.Fo, iSeason, iLandUse)
    end

    rac = obtain_value(iSeason, iLandUse, r_ac)

    # Correction for cold temperatures from page 4 column 1.
    correction = ifelse((Ts < 0.0), 1000.0 * exp(-Ts - 4), 0) # [s m-1] #mark
    rlux += correction
    rclx += correction
    rgsx += correction

    r_c = 1.0 / (1.0 / (rsmx) + 1.0 / rlux + 1.0 / (rdc + rclx) + 1.0 / (rac + rgsx))
    r_c = max(r_c, 10.0) # From "Results and conclusions" section to avoid extremely high deposition velocities over extremely rough surfaces.
    r_c = min(r_c, 9999.0)
    return r_c
end
