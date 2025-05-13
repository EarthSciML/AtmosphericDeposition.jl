export WesleySurfaceResistance

const inf = 1.e25

# r_i represents the minimum bulk canopy resistances for water vapor.
const r_i = SA_F32[inf 60 120 70 130 100 inf inf 80 100 150
                   inf inf inf inf 250 500 inf inf inf inf inf
                   inf inf inf inf 250 500 inf inf inf inf inf
                   inf inf inf inf 400 800 inf inf inf inf inf
                   inf 120 240 140 250 190 inf inf 160 200 300]

# r_lu signifies leaf cuticles in healthy vegetation and otherwise the outer surfaces in the upper canopy.
const r_lu = SA_F32[inf 2000 2000 2000 2000 2000 inf inf 2500 2000 4000
                    inf 9000 9000 9000 4000 8000 inf inf 9000 9000 9000
                    inf inf 9000 9000 4000 8000 inf inf 9000 9000 9000
                    inf inf inf inf 6000 9000 inf inf 9000 9000 9000
                    inf 4000 4000 4000 2000 3000 inf inf 4000 4000 8000]

# r_ac signifies transfer that depends only on canopy height and density.
const r_ac = SA_F32[100.0 200 100 2000 2000 2000 0 0 300 150 200
                    100 150 100 1500 2000 1700 0 0 200 120 140
                    100 10 100 1000 2000 1500 0 0 100 50 120
                    100 10 10 1000 2000 1500 0 0 50 10 50
                    100 50 80 1200 2000 1500 0 0 200 60 120]

# r_gs signifies uptake at the "ground" by soil, leaf litter, snow, water etc. 'S' and 'O' stand for SO2 and O3 respectively.
const r_gsS = SA_F32[400.0 150 350 500 500 100 0 1000 0 220 40
                     400 200 350 500 500 100 0 1000 0 300 400
                     400 150 350 500 500 200 0 1000 0 200 400
                     100 100 100 100 100 100 0 1000 100 100 50
                     500 150 350 500 500 200 0 1000 0 250 40]

const r_gsO = SA_F32[300.0 150 200 200 200 300 2000 400 1000 180 200
                     300 150 200 200 200 300 2000 400 800 180 200
                     300 150 200 200 200 300 2000 400 1000 180 20
                     600 3500 3500 3500 3500 3500 2000 400 3500 3500 3500
                     300 150 200 200 200 300 2000 400 1000 180 200]

# r_cl is meant to account for uptake pathways at the leaves, bark, etc. 'S' and 'O' stand for SO2 and O3 respectively.
const r_clS = SA_F32[inf 2000 2000 2000 2000 2000 inf inf 2500 2000 4000
                     inf 9000 9000 9000 2000 4000 inf inf 9000 9000 9000
                     inf inf 9000 9000 3000 6000 inf inf 9000 9000 9000
                     inf inf inf 9000 200 400 inf inf 9000 inf 9000
                     inf 4000 4000 4000 2000 3000 inf inf 4000 4000 8000]

const r_clO = SA_F32[inf 1000 1000 1000 1000 1000 inf inf 1000 1000 1000
                     inf 400 400 400 1000 600 inf inf 400 400 400
                     inf 1000 400 400 1000 600 inf inf 800 600 600
                     inf 1000 1000 400 1500 600 inf inf 800 1000 800
                     inf 1000 500 500 1500 700 inf inf 600 800 800]

# Holder for gas properties from Wesely (1989) Table 2.'
struct GasData
    Dh2oPerDx::AbstractFloat
    Hstar::AbstractFloat
    Fo::AbstractFloat
end

const So2Data = GasData(1.9, 1.e5, 0)
const O3Data = GasData(1.6, 0.01, 1)
const No2Data = GasData(1.6, 0.01, 0.1) # Wesely (1989) suggests that,
# in general, the sum of NO and NO2 should be considered rather
# than NO2 alone because rapid in-air chemical reactions can cause
# a significant change of NO and NO2 vertical fluxes between the
# surface and the point at which deposition velocities are applied,
# but the sum of NO and NO2 fluxes should be practically unchanged.
const NoData = GasData(1.3, 3.e-3, 0) # Changed according to Walmsley (1996)
const Hno3Data = GasData(1.9, 1.e14, 0)
const H2o2Data = GasData(1.4, 1.e5, 1)
const AldData = GasData(1.6, 15, 0)     # Acetaldehyde (aldehyde class)
const HchoData = GasData(1.3, 6.e3, 0)   # Formaldehyde
const OpData = GasData(1.6, 240, 0.1) # Methyl hydroperoxide (organic peroxide class)
const PaaData = GasData(2.0, 540, 0.1)  # Peroxyacetyl nitrate
const OraData = GasData(1.6, 4.e6, 0)   # Formic acid (organic acid class)
const Nh3Data = GasData(0.97, 2.e4, 0)  # Changed according to Walmsley (1996)
const PanData = GasData(2.6, 3.6, 0.1)  # Peroxyacetyl nitrate
const Hno2Data = GasData(1.6, 1.e5, 0.1) # Nitrous acid

# Obtain values from matrix using symbolic parameter iSeason and iLandUse
function obtain_value(iSeason, iLandUse, matrix)
    index = (iLandUse - 1) * 5 + iSeason
    interpolate_r_i = DataInterpolations.LinearInterpolation(vec(matrix), 1:55)
    interpolate_r_i(index)
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
                (3 * obtain_value(iSeason, iLandUse, r_lu) / (1.e-5 * Hstar + fo)) +
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
