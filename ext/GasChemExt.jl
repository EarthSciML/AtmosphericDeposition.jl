module GasChemExt

using AtmosphericDeposition, GasChem, EarthSciMLBase, DynamicQuantities, ModelingToolkit

function EarthSciMLBase.couple2(c::GasChem.SuperFastCoupler, d::AtmosphericDeposition.DrydepositionGCoupler)
    c, d = c.sys, d.sys

    operator_compose(convert(ODESystem, c), d, Dict(
        #c.SO2 => d.SO2, # SuperFast does not currently have SO2
        c.NO2 => d.NO2,
        c.O3 => d.O3,
        c.H2O2 => d.H2O2,
        c.CH2O => d.CH2O,
    ))
end

function EarthSciMLBase.couple2(c::GasChem.SuperFastCoupler, d::AtmosphericDeposition.WetdepositionCoupler)
    c, d = c.sys, d.sys

    operator_compose(convert(ODESystem, c), d, Dict(
        #c.SO2 => d.SO2, # SuperFast does not currently have SO2
        c.NO2 => d.NO2,
        c.O3 => d.O3,
        c.H2O2 => d.H2O2,
        c.CH2O => d.CH2O,
    ))
end

end
