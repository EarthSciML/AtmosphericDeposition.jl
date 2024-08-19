module GasChemExt

using AtmosphericDeposition, GasChem, EarthSciMLBase, Unitful, ModelingToolkit

export register_couplings_ext

# Use a global flag to track initialization and ensure that the coupling between SuperFast and NEI Emission is only initialized once
const couplings_registered_ext = Ref(false)

function register_couplings_ext()
    if couplings_registered_ext[]
        println("Couplings have already been registered.")
        return
    end
    
    println("Registering couplings in ext")
    @parameters t [unit = u"s"]

    register_coupling(SuperFast(t), Wetdeposition(t)) do c, e
        operator_compose(convert(ODESystem, c), e, Dict(
            c.SO2 => e.SO2,
            c.NO2 => e.NO2,
            c.O3 => e.O3,
            c.CH4 => e.CH4,
            c.CO => e.CO,
            c.DMS => e.DMS,
            c.ISOP => e.ISOP))
    end

    register_coupling(SuperFast(t), DrydepositionG(t)) do c, e
        operator_compose(convert(ODESystem, c), e, Dict(
            c.SO2 => e.SO2,
            c.NO2 => e.NO2,
            c.O3 => e.O3,
            c.NO => e.NO,
            c.H2O2 => e.H2O2,
            c.CH2O => e.CH2O))
    end
    
    couplings_registered_ext[] = true
    println("Coupling registry after registration: ", EarthSciMLBase.coupling_registry)
end

function __init__()
    println("Initializing EmissionExt module")
    register_couplings_ext()
end

end