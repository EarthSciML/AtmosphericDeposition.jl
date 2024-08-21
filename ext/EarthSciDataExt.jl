module EarthSciDataExt

using AtmosphericDeposition, EarthSciData, EarthSciMLBase, Unitful, ModelingToolkit

function EarthSciMLBase.couple2(d::AtmosphericDeposition.DrydepositionGCoupler, g::EarthSciData.GEOSFPCoupler)
    d, g = d.sys, g.sys

    @constants(
        PaPerhPa = 100, [unit = u"Pa/hPa", description="Conversion factor from hPa to Pa"],
        MW_air = 29, [unit = u"g/mol", description="dry air molar mass"],
        kgperg = 1e-3, [unit = u"kg/g", description="Conversion factor from g to kg"],
        R = 8.31446261815324, [unit = u"m^3*Pa/mol/K", description="Ideal gas constant"],
    )

    # ρA in DrydepositionG(t) are in units of "kg/m3".
    # ρ = P*M/(R*T)= Pa*(g/mol)/(m3*Pa/mol/K*K) = g/m3
    # Overall, ρA = P*M/(R*T)*kgperg

    d = param_to_var(d, :T, :z, :z₀, :u_star, :G, :ρA)
    ConnectorSystem([
            d.T ~ g.I3₊T,
            d.z ~ g.A1₊PBLH,
            d.z₀ ~ g.A1₊Z0M,
            d.u_star ~ g.A1₊USTAR,
            d.G ~ g. A1₊SWGDN,
            d.ρA ~ g.P * PaPerhPa/(g.I3₊T*R)*kgperg*MW_air,
        ], d, g)
end

function EarthSciMLBase.couple2(d::AtmosphericDeposition.WetdepositionCoupler, g::EarthSciData.GEOSFPCoupler)
    d, g = d.sys, g.sys

    @constants(
        PaPerhPa = 100, [unit = u"Pa/hPa", description="Conversion factor from hPa to Pa"],
        MW_air = 29, [unit = u"g/mol", description="dry air molar mass"],
        kgperg = 1e-3, [unit = u"kg/g", description="Conversion factor from g to kg"],
        R = 8.31446261815324, [unit = u"m^3*Pa/mol/K", description="Ideal gas constant"],
    )

    # ρ_air in Wetdeposition(t) are in units of "kg/m3".
    # ρ = P*M/(R*T)= Pa*(g/mol)/(m3*Pa/mol/K*K) = g/m3
    # Overall, ρ_air = P*M/(R*T)*kgperg

    d = param_to_var(d, :cloudFrac, :ρ_air)
    ConnectorSystem([
            d.cloudFrac ~ g.A3cld₊CLOUD,
            d.ρ_air ~ g.P * PaPerhPa/(g.I3₊T*R)*kgperg*MW_air,
        ], d, g)
end

end