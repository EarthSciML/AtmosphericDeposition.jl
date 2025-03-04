module EarthSciDataExt

using AtmosphericDeposition, EarthSciData, EarthSciMLBase, DynamicQuantities, ModelingToolkit

function EarthSciMLBase.couple2(d::AtmosphericDeposition.DrydepositionGCoupler, g::EarthSciData.GEOSFPCoupler)
    d, g = d.sys, g.sys

    @constants(
        MW_air = 0.029, [unit = u"kg/mol", description="dry air molar mass"],
        R = 8.31446261815324, [unit = u"m^3*Pa/mol/K", description="Ideal gas constant"],
        vK = 0.4, [description = "von Karman's constant"],
        Cp = 1000, [unit = u"W*s/kg/K", description="specific heat at constant pressure"],
        gg = 9.81, [unit = u"m*s^-2", description="gravitational acceleration"],
    )

    # ρA in DrydepositionG() are in units of "kg/m3".
    # ρ = P*M/(R*T)= Pa*(g/mol)/(m3*Pa/mol/K*K) = g/m3
    # Overall, ρA = P*M/(R*T)*kgperg

    # Monin-Obhukov length = -Air density * Cp * T(surface air) * Ustar^3/（vK   * g  * Sensible Heat flux）

    d = param_to_var(d, :T, :z, :z₀, :u_star, :G, :ρA, :L, :lev)

    ConnectorSystem([
            d.T ~ g.I3₊T,
            d.z ~ 0.1 * g.A1₊PBLH,
            d.z₀ ~ g.A1₊Z0M,
            d.u_star ~ g.A1₊USTAR,
            d.G ~ g.A1₊SWGDN,
            d.ρA ~ g.P/(g.I3₊T*R)*MW_air,
            d.L ~ -g.P/(g.I3₊T*R)*MW_air * Cp * g.A1₊TS * (g.A1₊USTAR)^3/(vK*gg*g.A1₊HFLUX),
            d.lev ~ g.lev,
        ], d, g)
end

function EarthSciMLBase.couple2(d::AtmosphericDeposition.WetdepositionCoupler, g::EarthSciData.GEOSFPCoupler)
    d, g = d.sys, g.sys

    @constants(
        MW_air = 0.029, [unit = u"kg/mol", description="dry air molar mass"],
        R = 8.31446261815324, [unit = u"m^3*Pa/mol/K", description="Ideal gas constant"],
        Vdr = 5.0, [unit = u"m/s", description="droplet velocity"],
    )

    # ρ_air in Wetdeposition() are in units of "kg/m3".
    # ρ = P*M/(R*T)= Pa*(g/mol)/(m3*Pa/mol/K*K) = g/m3
    # Overall, ρ_air = P*M/(R*T)*kgperg

    # From EMEP algorithm: P = QRAIN * Vdr * ρgas => QRAIN = P / Vdr / ρgas
    # kg*m-2*s-1/(m/s)/(kg/m3)

    d = param_to_var(d, :cloudFrac, :ρ_air, :qrain, :lev)
    ConnectorSystem([
            d.cloudFrac ~ g.A3cld₊CLOUD,
            d.ρ_air ~ g.P/(g.I3₊T*R)*MW_air,
            d.qrain ~ (g.A3mstE₊PFLCU + g.A3mstE₊PFLLSAN) / Vdr / (g.P/(g.I3₊T*R)*MW_air),
            d.lev ~ g.lev,
        ], d, g)
end

end