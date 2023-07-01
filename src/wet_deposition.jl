export WetDeposition, Wetdeposition, wd_defaults

"""
Calculate wet deposition based on formulas at
www.emep.int/UniDoc/node12.html.
Inputs are fraction of grid cell covered by clouds (cloudFrac),
rain mixing ratio (qrain), air density (ρ_air [kg/m3]),
and fall distance (Δz [m]).
Outputs are wet deposition rates for PM2.5, SO2, and other gases
(wdParticle, wdSO2, and wdOtherGas [1/s]).
"""

@constants A_wd = 5.2 [unit = u"m^3/kg/s"] # m3 kg-1 s-1; Empirical coefficient
@constants ρwater = 1000.0 [unit = u"kg*m^-3"]  # kg/m3
@constants Vdr = 5.0 [unit = u"m/s"] # raindrop fall speed, m/s

function WetDeposition(cloudFrac, qrain, ρ_air, Δz)
    #@constants A_wd = 5.2 [unit = u"m^3/kg/s"] # m3 kg-1 s-1; Empirical coefficient
	E = 0.1           # size-dependent collection efficiency of aerosols by the raindrops
	wSubSO2 = 0.15   # sub-cloud scavanging ratio
	wSubOther = 0.5  # sub-cloud scavanging ratio
	wInSO2 = 0.3     # in-cloud scavanging ratio
	wInParticle = 1.0 # in-cloud scavanging ratio
	wInOther = 1.4   # in-cloud scavanging ratio
	# @constants ρwater = 1000.0 [unit = u"kg*m^-3"]  # kg/m3
	# @constants Vdr = 5.0 [unit = u"m/s"] # raindrop fall speed, m/s

    # precalculated constant combinations
	AE = A_wd * E
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

	return wdParticle, wdSO2, wdOtherGas
end
wd_defaults = [A_wd => 5.2, ρwater => 1000.0, Vdr => 5.0]

# Add unit "ppb" to Unitful 
module myUnits
using Unitful
@unit ppb "ppb" Number 1 / 1000000000 false
end
Unitful.register(myUnits)

struct Wetdeposition <: EarthSciMLODESystem
    sys::ODESystem
    function Wetdeposition(t)
		@parameters cloudFrac = 0.5
		@parameters qrain = 0.5
		@parameters ρ_air = 1.204 [unit = u"kg*m^-3"]
		@parameters Δz = 200 [unit = u"m"]
        @parameters t [unit = u"s"]

        D = Differential(t)

        @variables SO2(t) [unit = u"ppb"]
        @variables O3(t) [unit = u"ppb"]

        eqs = [
            D(SO2) ~  -WetDeposition(cloudFrac, qrain, ρ_air, Δz)[2] * SO2
            D(O3) ~ -WetDeposition(cloudFrac, qrain, ρ_air, Δz)[3] * O3
        ]

        new(ODESystem(eqs, t, [SO2, O3], [cloudFrac, qrain, ρ_air, Δz]; name=:Wetdeposition))
    end
end 