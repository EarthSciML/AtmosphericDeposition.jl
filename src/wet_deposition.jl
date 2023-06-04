export WetDeposition

"""
Calculate wet deposition based on formulas at
www.emep.int/UniDoc/node12.html.
Inputs are fraction of grid cell covered by clouds (cloudFrac),
rain mixing ratio (qrain), air density (ρ_air [kg/m3]),
and fall distance (Δz [m]).
Outputs are wet deposition rates for PM2.5, SO2, and other gases
(wdParticle, wdSO2, and wdOtherGas [1/s]).
"""
function WetDeposition(cloudFrac, qrain, ρ_air, Δz)
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

    wdParticle = qrain * ρ_air * (AE + cloudFrac*(wInParticleVdrPerρwater/Δz))
	wdSO2 = (wSubSO2VdrPerρwater + cloudFrac*wSubSO2VdrPerρwater) * qrain * ρ_air / Δz
	wdOtherGas = (wSubOtherVdrPerρwater + cloudFrac*wSubOtherVdrPerρwater) * qrain * ρ_air / Δz

    return wdParticle, wdSO2, wdOtherGas
end
