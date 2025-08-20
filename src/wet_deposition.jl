export WetDeposition

@constants A_wd = 5.2 [unit = u"m^3/kg/s", description = "Empirical coefficient"]
@constants ρwater = 1000.0 [unit = u"kg*m^-3", description = "water density"]
@constants Vdr = 5.0 [unit = u"m/s", description = "raindrop fall speed"]

const lev_depth = [
    123.66503647136699,
    126.1276322862293,
    127.81793001768432,
    129.5320284164788,
    131.3096296124803,
    133.16517799136489,
    135.08768602272573,
    137.10796269057562,
    139.21228482592755,
    141.44376989951024,
    143.7574669454459,
    146.11933204200864,
    198.54224783503264,
    254.42969509717796,
    261.67340274464914,
    269.231033059717,
    277.3039413513152,
    285.9069077030176,
    446.4439568700527,
    469.32492993415644,
    494.6124439883338,
    522.6842502004765,
    554.6091144685597,
    591.0927240282208,
    633.2768512023067,
    682.3924743496791,
    741.7294428383548,
    1090.491426716162,
    1091.224127316289,
    1063.1347680566032,
    1038.8233086979671,
    1024.2528747829,
    1011.555297383451,
    1001.2550381718829,
    993.8747259925203,
    988.6668674973989,
    1000.6772919922369,
    1018.1853036083521,
    1038.4077259678052,
    1060.279636055635,
    1080.0011699235038,
    1105.2749220476362,
    1128.8308686418932,
    1150.9605352645412,
    1174.281709164923,
    1196.5658056986103,
    1220.477376106457,
    1248.055139801716,
    1289.3240338666546,
    1330.1699211514533,
    1371.1351955508508,
    1418.9073193940494,
    1469.6721949144776,
    1524.1720211194188,
    1589.3067557895556,
    1669.0982349095284,
    1748.6026372218985,
    1819.7984799666592,
    1873.109532132221,
    1911.520986089934,
    1927.4499053974941,
    1929.2166042272438,
    1929.6506325142327,
    1932.4995109991505,
    1946.1510743334002,
    1971.4849133077078,
    1995.6526762815556,
    2040.0593226684432,
    2158.2645757768187,
    2421.761343639795,
    3116.2739840471622,
    4343.478215407085
] # unit: m. The depth of each level.

function get_lev_depth(lev)
    j = 0
    for i in 1:length(lev_depth)
        j += ifelse(i <= lev, 1, 0)
    end
    result = lev_depth[j]
    return result
end
@register_symbolic get_lev_depth(lev)

"""
Calculate wet deposition based on formulas at
https://www.emep.int/publ/reports/2003/emep_report_1_part1_2003.pdf.
Inputs are fraction of grid cell covered by clouds (cloudFrac),
rain mixing ratio (qrain), air density (ρ_air [kg/m3]),
and fall distance (Δz [m]).
Outputs are wet deposition rates for PM2.5, SO2, and other gases
(wdParticle, wdSO2, and wdOtherGas [1/s]).
"""
function _WetDeposition(cloudFrac, qrain, ρ_air, Δz)
    E = 0.1           # size-dependent collection efficiency of aerosols by the raindrops
    wSubSO2 = 0.15   # sub-cloud scavanging ratio
    wSubOther = 0.5  # sub-cloud scavanging ratio
    wInSO2 = 0.3     # in-cloud scavanging ratio
    wInParticle = 1.0 # in-cloud scavanging ratio
    wInOther = 1.4   # in-cloud scavanging ratio

    i = ifelse(qrain > 0, 1, 0) # index to check if rain is present, avoid negative qrain

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

    wdParticle = i * qrain * ρ_air * (AE + cloudFrac * (wInParticleVdrPerρwater / Δz))
    wdSO2 = i * (wSubSO2VdrPerρwater + cloudFrac * wSubSO2VdrPerρwater) * qrain * ρ_air / Δz
    wdOtherGas = i * (wSubOtherVdrPerρwater + cloudFrac * wSubOtherVdrPerρwater) * qrain *
                 ρ_air / Δz

    return wdParticle, wdSO2, wdOtherGas
end
wd_defaults = [A_wd => 5.2, ρwater => 1000.0, Vdr => 5.0]

struct WetDepositionCoupler
    sys::Any
end

"""
Description: This is a box model used to calculate wet deposition based on formulas at EMEP model.
Build WetDeposition model

# Example

```julia
@parameters t
wd = WetDeposition(t)
```
"""
function WetDeposition(; name = :WetDeposition)
    params = @parameters(cloudFrac=0.5,
        [description="fraction of grid cell covered by clouds"],
        qrain=0.5,
        [description="rain mixing ratio"],
        ρ_air=1.204,
        [unit=u"kg*m^-3", description="air density"],
        lev=1,
        [description="level of the grid cell"],)

    @constants Δz_unit = 1 [unit = u"m", description = "unit depth"]

    vars = @variables begin
        k_particle(t), [unit=u"1/s"]
        k_SO2(t), [unit=u"1/s"]
        k_othergas(t), [unit=u"1/s"]
    end

    wdParticle, wdSO2,
    wdOtherGas = _WetDeposition(cloudFrac, qrain, ρ_air, get_lev_depth(lev) * Δz_unit)

    eqs = [k_particle ~ wdParticle, k_SO2 ~ wdSO2, k_othergas ~ wdOtherGas]

    ODESystem(
        eqs,
        t,
        vars,
        params;
        name = name,
        metadata = Dict(:coupletype => WetDepositionCoupler)
    )
end
