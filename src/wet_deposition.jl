export WetDeposition, Wetdeposition, wd_defaults

@constants A_wd = 5.2 [unit = u"m^3/kg/s", description = "Empirical coefficient"]
@constants ρwater = 1000.0 [unit = u"kg*m^-3", description = "water density"]
@constants Vdr = 5.0 [unit = u"m/s", description = "raindrop fall speed"]

"""
Calculate wet deposition based on formulas at
https://www.emep.int/publ/reports/2003/emep_report_1_part1_2003.pdf.
Inputs are fraction of grid cell covered by clouds (cloudFrac),
rain mixing ratio (qrain), air density (ρ_air [kg/m3]),
and fall distance (Δz [m]).
Outputs are wet deposition rates for PM2.5, SO2, and other gases
(wdParticle, wdSO2, and wdOtherGas [1/s]).
"""
function WetDeposition(cloudFrac, qrain, ρ_air, Δz)
    E = 0.1           # size-dependent collection efficiency of aerosols by the raindrops
    wSubSO2 = 0.15   # sub-cloud scavanging ratio
    wSubOther = 0.5  # sub-cloud scavanging ratio
    wInSO2 = 0.3     # in-cloud scavanging ratio
    wInParticle = 1.0 # in-cloud scavanging ratio
    wInOther = 1.4   # in-cloud scavanging ratio

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
    wdSO2 = (wSubSO2VdrPerρwater + cloudFrac * wSubSO2VdrPerρwater) * qrain * ρ_air / Δz
    wdOtherGas = (wSubOtherVdrPerρwater + cloudFrac * wSubOtherVdrPerρwater) * qrain * ρ_air / Δz

    return wdParticle, wdSO2, wdOtherGas
end
wd_defaults = [A_wd => 5.2, ρwater => 1000.0, Vdr => 5.0]

# Add unit "ppb" to Unitful 
module myUnits
using Unitful
@unit ppb "ppb" Number 1 / 1000000000 false
end
Unitful.register(myUnits)

"""
Description: This is a box model used to calculate wet deposition based on formulas at EMEP model.
Build Wetdeposition model
# Example
``` julia
	@parameters t 
	wd = Wetdeposition(t)
```
"""
function Wetdeposition(t)
    @parameters cloudFrac = 0.5 [description = "fraction of grid cell covered by clouds"]
    @parameters qrain = 0.5 [description = "rain mixing ratio"]
    @parameters ρ_air = 1.204 [unit = u"kg*m^-3", description = "air density"]
    @parameters Δz = 200 [unit = u"m", description = "fall distance"]

    D = Differential(t)

    @variables SO2(t) [unit = u"ppb"]
    @variables O3(t) [unit = u"ppb"]
    @variables NO2(t) [unit = u"ppb"]
    @variables CH4(t) [unit = u"ppb"]
    @variables CO(t) [unit = u"ppb"]
    @variables DMS(t) [unit = u"ppb"]
    @variables ISOP(t) [unit = u"ppb"]

    eqs = [
        D(SO2) ~ -WetDeposition(cloudFrac, qrain, ρ_air, Δz)[2] * SO2
        D(O3) ~ -WetDeposition(cloudFrac, qrain, ρ_air, Δz)[3] * O3
        D(NO2) ~ -WetDeposition(cloudFrac, qrain, ρ_air, Δz)[3] * NO2
        D(CH4) ~ -WetDeposition(cloudFrac, qrain, ρ_air, Δz)[3] * CH4
        D(CO) ~ -WetDeposition(cloudFrac, qrain, ρ_air, Δz)[3] * CO
        D(DMS) ~ -WetDeposition(cloudFrac, qrain, ρ_air, Δz)[3] * DMS
        D(ISOP) ~ -WetDeposition(cloudFrac, qrain, ρ_air, Δz)[3] * ISOP
    ]

    ODESystem(eqs, t, [SO2, O3, NO2, CH4, CO, DMS, ISOP], [cloudFrac, qrain, ρ_air, Δz]; name=:WetDeposition)
end
