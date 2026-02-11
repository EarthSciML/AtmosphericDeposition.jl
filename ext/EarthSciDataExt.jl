module EarthSciDataExt

using AtmosphericDeposition,
    EarthSciData, EarthSciMLBase, DynamicQuantities, ModelingToolkit

@constants(
    MW_air = 0.029,
    [unit = u"kg/mol", description = "Dry air molar mass"],
    vK = 0.4,
    [description = "von Karman's constant"],
    Cp = 1000,
    [unit = u"W*s/kg/K", description = "Specific heat at constant pressure"],
    R = 8.31446261815324,
    [unit = u"m^3*Pa/mol/K", description = "Ideal gas constant"],
    g = 9.81,
    [unit = u"m*s^-2", description = "Gravitational acceleration"],
)

air_density(P, T) = P / (T * R) * MW_air

# Monin-Obhukov length = -Air density * Cp * T(surface air) * Ustar^3/（vK   * g  * Sensible Heat flux）
MoninObhukovLength(ρ_air, Ts, u_star, HFLUX) = -ρ_air * Cp * Ts * (u_star)^3 / (vK * g * HFLUX)

function EarthSciMLBase.couple2(
        d::AtmosphericDeposition.DryDepositionGasCoupler,
        gp::EarthSciData.GEOSFPCoupler
    )
    d, gp = d.sys, gp.sys

    d = param_to_var(d, :Ts, :z, :z₀, :u_star, :G, :ρA, :L, :lev)

    return ConnectorSystem(
        [
            d.Ts ~ gp.A1₊TS,
            d.z ~ 0.1 * gp.A1₊PBLH, # the surface layer height is 10% of the boundary layer height
            d.z₀ ~ gp.A1₊Z0M,
            d.u_star ~ gp.A1₊USTAR,
            d.G ~ gp.A1₊SWGDN,
            d.ρA ~ air_density(gp.P, gp.I3₊T),
            d.L ~ MoninObhukovLength(d.ρA, gp.A1₊TS, gp.A1₊USTAR, gp.A1₊HFLUX),
            d.lev ~ gp.lev,
        ],
        d,
        gp
    )
end

function EarthSciMLBase.couple2(
        d::AtmosphericDeposition.DryDepositionAerosolCoupler,
        gp::EarthSciData.GEOSFPCoupler
    )
    d, gp = d.sys, gp.sys

    d = param_to_var(d, :Ts, :z, :z₀, :u_star, :ρA, :L, :lev)

    return ConnectorSystem(
        [
            d.Ts ~ gp.A1₊TS,
            d.z ~ 0.1 * gp.A1₊PBLH, # the surface layer height is 10% of the boundary layer height
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

    d = param_to_var(d, :cloudFrac, :ρ_air, :qrain, :lev)
    return ConnectorSystem(
        [
            d.cloudFrac ~ g.A3cld₊CLOUD,
            d.ρ_air ~ air_density(g.P, g.I3₊T),
            d.qrain ~ (g.A3mstE₊PFLCU + g.A3mstE₊PFLLSAN) / Vdr / (g.P / (g.I3₊T * R) * MW_air),
            d.lev ~ g.lev,
        ],
        d,
        g
    )
end

end
