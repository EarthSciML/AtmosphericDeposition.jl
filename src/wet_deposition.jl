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

struct WetdepositionCoupler
    sys
end

"""
Description: This is a box model used to calculate wet deposition based on formulas at EMEP model.
Build Wetdeposition model
# Example
``` julia
	@parameters t 
	wd = Wetdeposition(t)
```
"""
function Wetdeposition(t; name=:Wetdeposition)
    params = @parameters(
        cloudFrac = 0.5, [description = "fraction of grid cell covered by clouds"],
        qrain = 0.5, [description = "rain mixing ratio"],
        ρ_air = 1.204, [unit = u"kg*m^-3", description = "air density"],
        Δz = 200, [unit = u"m", description = "fall distance"],
    )

    D = Differential(t)

    vars = @variables(
        SO2(t), [unit = u"nmol/mol"],
        O3(t), [unit = u"nmol/mol"],
        NO2(t), [unit = u"nmol/mol"],
        CH4(t), [unit = u"nmol/mol"],
        CO(t), [unit = u"nmol/mol"],
        DMS(t), [unit = u"nmol/mol"],
        ISOP(t), [unit = u"nmol/mol"],
    )

    eqs = [
        D(SO2) ~ -WetDeposition(cloudFrac, qrain, ρ_air, Δz)[2] * SO2
        D(O3) ~ -WetDeposition(cloudFrac, qrain, ρ_air, Δz)[3] * O3
        D(NO2) ~ -WetDeposition(cloudFrac, qrain, ρ_air, Δz)[3] * NO2
        D(CH4) ~ -WetDeposition(cloudFrac, qrain, ρ_air, Δz)[3] * CH4
        D(CO) ~ -WetDeposition(cloudFrac, qrain, ρ_air, Δz)[3] * CO
        D(DMS) ~ -WetDeposition(cloudFrac, qrain, ρ_air, Δz)[3] * DMS
        D(ISOP) ~ -WetDeposition(cloudFrac, qrain, ρ_air, Δz)[3] * ISOP
    ]

    ODESystem(eqs, t, vars, params; name=name,
        metadata=Dict(:coupletype => WetdepositionCoupler))
end
