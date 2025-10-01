module GasChemExt

using AtmosphericDeposition, GasChem, EarthSciMLBase, DynamicQuantities, ModelingToolkit

function EarthSciMLBase.couple2(
        c::GasChem.SuperFastCoupler,
        d::AtmosphericDeposition.DryDepositionGasCoupler
)
    c, d = c.sys, d.sys

    operator_compose(
        convert(ODESystem, c),
        d,
        Dict(
            #c.SO2 => d.SO2 => c.SO2, # SuperFast does not currently have SO2
            c.HNO3 => d.k_HNO3 => -c.HNO3,
            c.NO2 => d.k_NO2 => -c.NO2,
            c.NO => d.k_NO2 => -c.NO,
            c.O3 => d.k_O3 => -c.O3,
            c.H2O2 => d.k_H2O2 => -c.H2O2,
            c.CH2O => d.k_HCHO => -c.CH2O
        )
    )
end

function EarthSciMLBase.couple2(
        c::GasChem.SuperFastCoupler,
        d::AtmosphericDeposition.WetDepositionCoupler
)
    c, d = c.sys, d.sys

    operator_compose(
        convert(ODESystem, c),
        d,
        Dict(
            #c.SO2 => d.k_SO2 => c.SO2, # SuperFast does not currently have SO2
            c.HNO3 => d.k_othergas => -c.HNO3,
            c.H2O2 => d.k_othergas => -c.H2O2,
            c.CH2O => d.k_othergas => -c.CH2O,
            c.CH3OOH => d.k_othergas => -c.CH3OOH
        )
    )
end

end
