module AerosolExt
using AtmosphericDeposition, Aerosol
using EarthSciMLBase, ModelingToolkit

function EarthSciMLBase.couple2(a::Aerosol.ElementalCarbonCoupler,
    d::AtmosphericDeposition.DryDepositionAerosolCoupler)
    a, d = a.sys, d.sys

    operator_compose(a, d, Dict(
        a.EC => d.k => -a.EC,
    ))
end

function EarthSciMLBase.couple2(a::Aerosol.ElementalCarbonCoupler,
    d::AtmosphericDeposition.WetDepositionCoupler)
    a, d = a.sys, d.sys

    operator_compose(a, d, Dict(
        a.EC => d.k_particle => -a.EC,
    ))
end

end
