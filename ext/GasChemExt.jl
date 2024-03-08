module GasChemExt

using AtmosphericDeposition, GasChem
Base.:(+)(w::Wetdeposition, b::SuperFast) = operator_compose(b, w)
Base.:(+)(b::SuperFast, w::Wetdeposition) = w + b

Base.:(+)(d::DrydepositionG, b::SuperFast) = operator_compose(b, d)
Base.:(+)(b::SuperFast, d::DrydepositionG) = d + b

end