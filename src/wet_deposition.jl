export WetDeposition, Wetdeposition

"""
Calculate wet deposition based on formulas at
www.emep.int/UniDoc/node12.html.
Inputs are fraction of grid cell covered by clouds (cloudFrac),
rain mixing ratio (qrain), air density (ρ_air [kg/m3]),
and fall distance (Δz [m]).
Outputs are wet deposition rates for PM2.5, SO2, and other gases
(wdParticle, wdSO2, and wdOtherGas [1/s]).
"""
function WetDeposition(cloudFrac, qrain, ρ_air, Δz, output)
    A = 5.2u"m^3/kg/s"          # m3 kg-1 s-1; Empirical coefficient
	E = 0.1           # size-dependent collection efficiency of aerosols by the raindrops
	wSubSO2 = 0.15   # sub-cloud scavanging ratio
	wSubOther = 0.5  # sub-cloud scavanging ratio
	wInSO2 = 0.3     # in-cloud scavanging ratio
	wInParticle = 1.0 # in-cloud scavanging ratio
	wInOther = 1.4   # in-cloud scavanging ratio
	ρwater = 1000.0u"kg*m^-3"   # kg/m3
	Vdr = 5.0u"m/s"        # raindrop fall speed, m/s

    # precalculated constant combinations
	AE = A * E
	wSubSO2VdrPerρwater = wSubSO2 * Vdr / ρwater
	wSubOtherVdrPerρwater = wSubOther * Vdr / ρwater
	wInSO2VdrPerρwater = wInSO2 * Vdr / ρwater
	wInParticleVdrPerρwater = wInParticle * Vdr / ρwater
	wInOtherVdrPerρwater = wInOther * Vdr / ρwater

    # wdParticle (subcloud) = A * P / Vdr * E; P = QRAIN * Vdr * ρgas => wdParticle = A * QRAIN * ρgas * E
	# wdGas (subcloud) = wSub * P / Δz / ρwater = wSub * QRAIN * Vdr * ρgas / Δz / ρwater
	# wd (in-cloud) = wIn * P / Δz / ρwater = wIn * QRAIN * Vdr * ρgas / Δz / ρwater

    wdParticle = qrain * ρ_air * (AE + cloudFrac * (wInParticleVdrPerρwater / Δz))
	wdSO2 = (wSubSO2VdrPerρwater + cloudFrac*wSubSO2VdrPerρwater) * qrain * ρ_air / Δz
	wdOtherGas = (wSubOtherVdrPerρwater + cloudFrac*wSubOtherVdrPerρwater) * qrain * ρ_air / Δz

	if output == 1
		return wdParticle
	elseif output == 2 
		return wdSO2
	else 
		return wdOtherGas
	end
end

using EarthSciMLBase
using ModelingToolkit
using Unitful
# Add unit "ppb" to Unitful 
module MyUnits
using Unitful
@unit ppb "ppb" Number 1 / 1000000000 false
end
Unitful.register(MyUnits)

struct Wetdeposition <: EarthSciMLODESystem
    sys::ODESystem
    function Wetdeposition(t, cloudFrac, qrain, ρ_air, Δz)
		@parameters k1 = WetDeposition(cloudFrac, qrain, ρ_air, Δz, 2) * 1u"s" [unit = u"s^-1"]
		@parameters k2 = WetDeposition(cloudFrac, qrain, ρ_air, Δz, 3) * 1u"s" [unit = u"s^-1"]
        #@parameters t [unit = u"s"]

        D = Differential(t)

        @variables SO2(t) [unit = u"ppb"]
        @variables O3(t) [unit = u"ppb"]

        eqs = [
            D(SO2) ~  -k1 * SO2
            D(O3) ~ -k2 * O3
        ]

        new(ODESystem(eqs, t, [SO2, O3], [k1,k2]; name=:Wetdeposition))
    end
end 

