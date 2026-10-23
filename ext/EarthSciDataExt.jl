module EarthSciDataExt

using AtmosphericDeposition,
    EarthSciData, EarthSciMLBase, DynamicQuantities, ModelingToolkit
using ModelingToolkit: t

@constants(MW_air=0.029,
    [unit=u"kg/mol", description="Dry air molar mass"],
    vK=0.4,
    [description="von Karman's constant"],
    Cp=1000,
    [unit=u"W*s/kg/K", description="Specific heat at constant pressure"],
    R=8.31446261815324,
    [unit=u"m^3*Pa/mol/K", description="Ideal gas constant"],
    g=9.81,
    [unit=u"m*s^-2", description="Gravitational acceleration"],
    P_unit = 1,
    [unit=u"Pa", description="Unit for pressure"])

air_density(P, T) = P / (T * R) * MW_air

# Monin-Obhukov length = -Air density * Cp * T(surface air) * Ustar^3/（vK   * g  * Sensible Heat flux）
MoninObhukovLength(ρ_air, Ts, u_star, HFLUX) = -ρ_air * Cp * Ts * (u_star)^3 / (vK * g * HFLUX)

# First level pressure thickness using the first 2 values of Ap and Bp.
# Ap(2) = 0.04804826 hPa = 4.804826 Pa; the hPa value was inserted as Pa (100× too small).
first_level_pressure_thickness(P) = -4.804826 * P_unit + P * 0.015048

# First level midpoint height above ground from the hypsometric relation, using the
# mean of the level-1 and 2 m virtual temperatures
@constants(Rd = 287.05,
    [unit = u"J/(kg*K)", description = "Dry-air gas constant"],
    g₀ = 9.80665,
    [unit = u"m*s^-2", description = "Standard gravity"])
first_level_midpoint_height(PS, T, QV, T2M, QV2M) =
    Rd / g₀ * 0.5 * (T * (1 + 0.61 * QV) + T2M * (1 + 0.61 * QV2M)) *
    log(PS / (PS - first_level_pressure_thickness(PS) / 2))

# Gas dry deposition (fractional / mosaic) bound to GEOS-FP. Binds surface
# meteorology, the date-driven `season`, and the 11 land-use area fractions
# (looked up from the bundled CONUS NetCDF via
# `landuse_frac_at(gp.lon, gp.lat, i)`).
function EarthSciMLBase.couple2(
        d::AtmosphericDeposition.DryDepositionGasCoupler,
        gp::EarthSciData.GEOSFPCoupler
    )
    d, gp = d.sys, gp.sys

    d = param_to_var(d,
        :Ts, :z, :del_P, :z₀, :u_star, :G, :ρA, :L, :lev,
        :season,
        :f_urban, :f_agricultural, :f_range, :f_deciduous, :f_coniferous,
        :f_mixedforest, :f_water, :f_barren, :f_wetland, :f_rangeag,
        :f_rockyshrubs,
    )

    return ConnectorSystem(
        [
            d.Ts ~ gp.A1₊TS,
            d.z ~ first_level_midpoint_height(gp.I3₊PS, gp.I3₊T, gp.I3₊QV, gp.A1₊T2M, gp.A1₊QV2M),
            d.del_P ~ first_level_pressure_thickness(gp.I3₊PS),
            d.z₀ ~ gp.A1₊Z0M,
            d.u_star ~ gp.A1₊USTAR,
            d.G ~ gp.A1₊SWGDN,
            d.ρA ~ air_density(gp.P, gp.I3₊T),
            d.L ~ MoninObhukovLength(d.ρA, gp.A1₊TS, gp.A1₊USTAR, gp.A1₊HFLUX),
            d.lev ~ gp.lev,
            d.season ~ AtmosphericDeposition.season_at(gp.t_ref + t),
            d.f_urban        ~ AtmosphericDeposition.landuse_frac_at(gp.lon, gp.lat, 1),
            d.f_agricultural ~ AtmosphericDeposition.landuse_frac_at(gp.lon, gp.lat, 2),
            d.f_range        ~ AtmosphericDeposition.landuse_frac_at(gp.lon, gp.lat, 3),
            d.f_deciduous    ~ AtmosphericDeposition.landuse_frac_at(gp.lon, gp.lat, 4),
            d.f_coniferous   ~ AtmosphericDeposition.landuse_frac_at(gp.lon, gp.lat, 5),
            d.f_mixedforest  ~ AtmosphericDeposition.landuse_frac_at(gp.lon, gp.lat, 6),
            d.f_water        ~ AtmosphericDeposition.landuse_frac_at(gp.lon, gp.lat, 7),
            d.f_barren       ~ AtmosphericDeposition.landuse_frac_at(gp.lon, gp.lat, 8),
            d.f_wetland      ~ AtmosphericDeposition.landuse_frac_at(gp.lon, gp.lat, 9),
            d.f_rangeag      ~ AtmosphericDeposition.landuse_frac_at(gp.lon, gp.lat, 10),
            d.f_rockyshrubs  ~ AtmosphericDeposition.landuse_frac_at(gp.lon, gp.lat, 11),
        ],
        d,
        gp,
    )
end

function EarthSciMLBase.couple2(
        d::AtmosphericDeposition.DryDepositionAerosolCoupler,
        gp::EarthSciData.GEOSFPCoupler
    )
    d, gp = d.sys, gp.sys

    d = param_to_var(d, :Ts, :z, :del_P, :z₀, :u_star, :ρA, :L, :lev)

    return ConnectorSystem(
        [
            d.Ts ~ gp.A1₊TS,
            d.z ~ 0.1 * gp.A1₊PBLH, # Ra reference height only (surface layer = 10% of PBL)
            d.del_P ~ first_level_pressure_thickness(gp.I3₊PS), # level-1 box thickness for k = v·g·ρA/del_P
            d.z₀ ~ gp.A1₊Z0M,
            d.u_star ~ gp.A1₊USTAR,
            d.ρA ~ air_density(gp.P, gp.I3₊T),
            d.L ~ MoninObhukovLength(d.ρA, gp.A1₊TS, gp.A1₊USTAR, gp.A1₊HFLUX),
            d.lev ~ gp.lev,
        ],
        d,
        gp
    )
end

function EarthSciMLBase.couple2(
        d::AtmosphericDeposition.WetDepositionCoupler,
        g::EarthSciData.GEOSFPCoupler
    )
    d, g = d.sys, g.sys

    @constants(Vdr = 5.0, [unit = u"m/s", description = "droplet velocity"])

    # From EMEP algorithm: P = QRAIN * Vdr * ρgas => QRAIN = P / Vdr / ρgas
    # kg*m-2*s-1/(m/s)/(kg/m3)

    # lev is deliberately NOT promoted or linked: with the EMEP fixed-h_s rate
    # form the rates no longer divide by the local layer thickness, so Delta-z (and
    # lev, which only entered through get_lev_depth(lev) * Delta-z_unit) is unused
    # in the WetDeposition equations and is pruned at convert time; promoting
    # or linking it here throws "variable lev does not exist".
    d = param_to_var(d, :cloudFrac, :ρ_air, :qrain)
    return ConnectorSystem(
        [
            d.cloudFrac ~ g.A3cld₊CLOUD,
            d.ρ_air ~ air_density(g.P, g.I3₊T),
            d.qrain ~ (g.A3mstE₊PFLCU + g.A3mstE₊PFLLSAN) / Vdr / (g.P / (g.I3₊T * R) * MW_air),
        ],
        d,
        g
    )
end

end
