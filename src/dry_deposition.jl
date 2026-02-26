export DryDepositionGas, DryDepositionAerosol

@constants g = 9.81 [unit = u"m*s^-2", description = "gravitational acceleration"]
@constants κ = 0.4 [description = "von Karman constant"]
@constants k = 1.3806488e-23 [unit = u"m^2*kg*s^-2/K", description = "Boltzmann constant"]
@constants M_air = 28.97e-3 [unit = u"kg/mol", description = "molecular weight of air"]
@constants R = 8.3144621 [unit = u"kg*m^2*s^−2*K^-1*mol^-1", description = "Gas constant"]

@constants unit_m = 1 [unit = u"m"]
"""
Function Ra calculates aerodynamic resistance to dry deposition
where z is the top of the surface layer [m], z₀ is the roughness length [m], u_star is friction velocity [m/s], and L is Monin-Obukhov length [m]
Based on Seinfeld and Pandis (2006) [Seinfeld, J.H. and Pandis, S.N. (2006) Atmospheric Chemistry and Physics: From Air Pollution to Climate Change. 2nd Edition, John Wiley & Sons, New York.]
equation 19.13 & 19.14.
"""
function ra(z, z₀, u_star, L)
    ζ = ifelse((L / unit_m == 0), 0, z / L)
    ζ₀ = ifelse((L / unit_m == 0), 0, z₀ / L)
    # Whether stable or neutral
    rₐ_1 = ifelse(
        (0 < ζ),
        1 / (κ * u_star) * (log(z / z₀) + 4.7 * (ζ - ζ₀)),
        1 / (κ * u_star) * log(z / z₀)
    )

    i = ifelse((ζ < 0), 1, 0) # Use index i to avoid a DomainError when calculating (1 - 15 * ζ) ^ 1/4 calculation for ζ > 0. η₀ and η are only needed when ζ < 0 but both sides are evaluated when using the ifelse function.
    η₀ = (1 - 15 * ζ₀ * i)^(1 / 4)
    η = (1 - 15 * ζ * i)^(1 / 4)

    # Whether unstable
    rₐ = ifelse(
        (ζ < 0),
        1 / (κ * u_star) * [
            log(z / z₀) +
                log(((η₀^2 + 1) * (η₀ + 1)^2) / ((η^2 + 1) * (η + 1)^2)) +
                2 * (atan(η) - atan(η₀)),
        ][1],
        rₐ_1
    ) #the [1] is to pass the ModelingToolkit unit check
    return rₐ
end

@constants unit_T = 1 [unit = u"K", description = "unit one for temperature"]
@constants unit_convert_mu = 1 [unit = u"kg/m/s", description = "unit one for mu"]
"""
Function mu calculates the dynamic viscosity of air [kg m-1 s-1] where T is temperature [K].
"""
function mu(T)
    return (1.458 * 10^-6 * (T / unit_T)^(3 / 2) / ((T / unit_T) + 110.4)) * unit_convert_mu
end

"""
Function mfp calculates the mean free path of air [m]
where T is temperature [K] P is pressure [Pa], and Mu is dynamic viscosity [kg/(m s)].
From Seinfeld and Pandis (2006) equation 9.6
"""
function mfp(T, P, μ)
    return 2 * μ / (P * (8 * M_air / (pi * R * T))^0.5)
end

"""
Function cc calculates the Cunnningham slip correction factor
where Dp is particle diameter [m], T is temperature [K], and P is pressure [Pa].
From Seinfeld and Pandis (2006) equation 9.34.
"""
function cc(Dₚ, T, P, μ)
    λ = mfp(T, P, μ)
    return 1 + 2 * λ / Dₚ * (1.257 + 0.4 * exp(-1.1 * Dₚ / (2 * λ)))
end

@constants unit_v = 1 [unit = u"m/s", description = "unit one for speed"]
@constants v_zero = 0 [unit = u"m/s", description = "zero velocity"]
"""
Function vs calculates the terminal setting velocity of a
particle where Dp is particle diameter [m], ρₚ is particle density [kg/m3], Cc is the Cunningham slip correction factor, and μ is air dynamic viscosity [kg/(s m)].
From equation 9.42 in Seinfeld and Pandis (2006).
"""
function vs(Dₚ, ρₚ, Cc, μ)
    return ifelse((Dₚ > 20.0e-6 * unit_m), 99999999 * unit_v, Dₚ^2 * ρₚ * g * Cc / (18 * μ))
    # Particle diameter Dₚ greater than 20um; Stokes settling no longer applies.
end

"""
Function dParticle calculates the brownian diffusivity of a particle [m2/s] using the Stokes-Einstein-Sutherland relation
(Seinfeld and Pandis eq. 9.73)
where T is air temperature [K], P is pressure [Pa], Dp is particle diameter [m], and μ is air dynamic viscosity [kg/(s m)]
"""
function dParticle(T, P, Dₚ, Cc, μ)
    return k * T * Cc / (3 * pi * μ * Dₚ)
end

@constants T_unitless = 1 [unit = u"K^-1", description = "used to offset temperature unit"]
@constants unit_dH2O = 1 [unit = u"m^2/s", description = "unit for molecular diffusivity"]
"""
Function dH2O calculates molecular diffusivity of water vapor in air [m2/s] where T is temperature [K]
using a regression fit to data in Bolz and Tuve (1976) found here: http://www.cambridge.org/us/engineering/author/nellisandklein/downloads/examples/EXAMPLE_9.2-1.pdf
"""
function dH2O(T)
    return (-2.775e-6 + 4.479e-8 * T * T_unitless + 1.656e-10 * (T * T_unitless)^2) *
        unit_dH2O
end

"""
Function sc computes the dimensionless Schmidt number,
where μ is dynamic viscosity of air [kg/(s m)], ρ is air density [kg/m3], and D is the molecular diffusivity of the gas speciesof interest [m2/s]
"""
function sc(μ, ρ, D)
    return μ / (ρ * D)
end

"""
Function stSmooth computes the dimensionless Stokes number for dry deposition of particles on smooth surfaces or surfaces with bluff roughness elements,
where vs is settling velocity [m/s], u_star is friction velocity [m/s], μ is dynamic viscosity of air [kg/(s m)], and ρ is air density [kg/m3],
based on Seinfeld and Pandis (2006) equation 19.23.
"""
function stSmooth(vₛ, u_star, μ, ρ)
    return vₛ * u_star^2 * ρ / (g * μ)
end

"""
Function stVeg computes the dimensionless Stokes number for dry deposition of particles on vegetated surfaces,
where vs is settling velocity [m/s], u_star is friction velocity [m/s], and A is the characteristic collector radius [m],
based on Seinfeld and Pandis (2006) equation 19.24.
"""
function stVeg(vₛ, u_star, A)
    return vₛ * u_star / (g * A)
end

"""
Function RbGas calculates the quasi-laminar sublayer resistance to dry deposition for a gas species [s/m],
where Sc is the dimensionless Schmidt number and u_star is the friction velocity [m/s].
From Seinfeld and Pandis (2006) equation 19.17.
"""
function RbGas(Sc, u_star)
    return 5 * Sc^(2 / 3) / u_star
end

@enum seinfeldLandUse begin
    seinfeldEvergreen = 1 #	0. Evergreen-needleleaf trees
    seinfeldDeciduous # Deciduous broadleaf trees
    seinfeldGrass # Grass
    seinfeldDesert # Desert
    seinfeldShrubs # Shrubs and interrupted woodlands
end

@enum seinfeldSeason begin
    seinfeldMidsummer = 1 # Midsummer with lush vegetation
    seinfeldAutumn # Autumn with cropland not harvested
    seinfeldLateAutumn # Late autumn after frost, no snow
    seinfeldWinter # Winter, snow on ground
    seinfeldTransitional # transitional
end

"""
Values for the characteristic radii of collectors [m]
where the columns are land use categories and the rows are seasonal categories.
Given in Seinfeld and Pandis Table 19.2

Land use options are given in SeinfeldLandUse and season options are given in SeinfeldSeason.
"""
z₀_table = [
    0.8 1.05 0.1 0.04 0.1
    0.9 1.05 0.1 0.04 0.1
    0.9 0.95 0.05 0.04 0.1
    0.9 0.55 0.02 0.04 0.1
    0.8 0.75 0.05 0.04 0.1
] # unit:[m]

function A_table(iSeason, iLandUse)
    return SA_F32[
        2.0 5.0 2.0 Inf 10.0
        2.0 5.0 2.0 Inf 10.0
        2.0 10.0 5.0 Inf 10.0
        2.0 10.0 2.0 Inf 10.0
        2.0 5.0 2.0 Inf 10.0
    ][
        iSeason,
        iLandUse,
    ] * 1.0e-3 # unit:[mm]
end # unit:[mm]
@register_symbolic A_table(iSeason, iLandUse)
A_table(::DynamicQuantities.Quantity, ::DynamicQuantities.Quantity) = 1.0
ModelingToolkit.get_unit(::typeof(A_table)) = 1.0 # TODO(CT): Can't figure out how to set units to meters.

α_table(iLandUse) = SA_F32[1.0 0.8 1.2 50.0 1.3][iLandUse]
@register_symbolic α_table(iLandUse)
ModelingToolkit.get_unit(::typeof(α_table)) = 1.0
α_table(::DynamicQuantities.Quantity) = 1.0

γ_table(iLandUse) = SA_F32[0.56 0.56 0.54 0.54 0.54][iLandUse]
@register_symbolic γ_table(iLandUse)
ModelingToolkit.get_unit(::typeof(γ_table)) = 1.0
γ_table(::DynamicQuantities.Quantity) = 1.0

"""
Function RbParticle calculates the quasi-laminar sublayer resistance to dry deposition for a particles [s/m],
where Sc is the dimensionless Schmidt number, u_star is the friction velocity [m/s], St is the dimensionless Stokes number,
Dp is particle diameter [m], and iSeason and iLandUse are season and land use indexes, respectively.
From Seinfeld and Pandis (2006) equation 19.27.
"""
function RbParticle(Sc, u_star, St, Dₚ, iSeason, iLandUse)
    α = α_table(iLandUse)
    γ = γ_table(iLandUse)
    A = A_table(iSeason, iLandUse) * unit_m
    R1 = exp(-St^0.5)
    term_1 = Sc^(-γ)
    term_2 = (St / (α + St))^2
    term_3 = 1 / 2 * (Dₚ / A)^2
    return 1 / (3 * u_star * (term_1 + term_2 + term_3) * R1)
end

@constants G_unitless = 1 [
    unit = u"m^2/W",
    description = "used to offset the unit of irradiation",
]
@constants Rc_unit = 1 [unit = u"s/m", description = "unit for surface resistance"]
"""
Function DryDepGas calculates dry deposition velocity [m/s] for a gas species,
where z is the height of the surface layer [m], zo is roughness length [m], u_star is friction velocity [m/s],
L is Monin-Obukhov length [m], T is surface air temperature [K], ρA is air density [kg/m3]
gasData is data about the gas species for surface resistance calculations, G is solar
irradiation [W m-2], Θ is the slope of the local terrain [radians], iSeason and iLandUse are indexes for the season and land use,
dew and rain indicate whether there is dew or rain on the ground, and isSO2 and isO3 indicate whether the gas species of interest is either SO2 or O3, respectively.
Based on Seinfeld and Pandis (2006) equation 19.2.
"""
function DryDepGas(
        lev,
        z,
        z₀,
        u_star,
        L,
        ρA,
        gasData::GasData,
        G,
        Ts,
        θ,
        iwesleySeason,
        iwesleyLandUse,
        rain::Bool,
        dew::Bool,
        isSO2::Bool,
        isO3::Bool
    )
    Ra = ra(z, z₀, u_star, L)
    μ = mu(Ts)
    Dg = dH2O(Ts) / gasData.Dh2oPerDx # Diffusivity of gas of interest [m2/s]
    Sc = sc(μ, ρA, Dg)
    Rb = RbGas(Sc, u_star)
    Rc = WesleySurfaceResistance(
        gasData,
        G * G_unitless,
        (Ts * T_unitless - 273),
        θ,
        iwesleySeason,
        iwesleyLandUse,
        rain::Bool,
        dew::Bool,
        isSO2::Bool,
        isO3::Bool
    ) * Rc_unit
    i = ifelse(lev == 1, 1, 0)
    result = i / (Ra + Rb + Rc)
    return result
end

"""
Function DryDepParticle calculates particle dry deposition velocity [m/s]
where z is the height of the surface layer [m], zo is roughness length [m], u_star is friction velocity [m/s], L is Monin-Obukhov length [m],
Dp is particle diameter [m], Ts is surface air temperature [K], P is pressure [Pa], ρParticle is particle density [kg/m3], ρAir is air density [kg/m3],
and iSeason and iLandUse are indexes for the season and land use.
Based on Seinfeld and Pandis (2006) equation 19.7.
"""
function DryDepParticle(
        lev, z, z₀, u_star, L, Dp, Ts, P, ρParticle, ρA,
        iSeinfeldSeason, iWesleySeason, iSeinfeldLandUse, iWesleyLandUse
    )
    Ra = ra(z, z₀, u_star, L)
    μ = mu(Ts)
    Cc = cc(Dp, Ts, P, μ)
    Vs = vs(Dp, ρParticle, Cc, μ)
    St = ifelse(
        iSeinfeldLandUse == Int(seinfeldDesert),
        stSmooth(Vs, u_star, μ, ρA),
        stVeg(Vs, u_star, A_table(iSeinfeldSeason, iSeinfeldLandUse) * unit_m)
    )
    D = dParticle(Ts, P, Dp, Cc, μ)
    Sc = sc(μ, ρA, D)
    Rb = RbParticle(Sc, u_star, St, Dp, iSeinfeldSeason, iSeinfeldLandUse)
    return ifelse(lev == 1, 1 / (Ra + Rb + Ra * Rb * Vs) + Vs, v_zero)
end

defaults = [
    g => 9.81,
    κ => 0.4,
    k => 1.3806488e-23,
    M_air => 28.97e-3,
    R => 8.3144621,
    unit_T => 1,
    unit_convert_mu => 1,
    T_unitless => 1,
    unit_dH2O => 1,
    unit_m => 1,
    G_unitless => 1,
    Rc_unit => 1,
    unit_v => 1,
    v_zero => 0,
]

struct DryDepositionGasCoupler
    sys::Any
end

"""
DescriptionGas: This is a box model used to calculate the gas species concentration rate changed by dry deposition.
Build Dry deposition model (gas)

# Example

```julia
@parameters t
d = DrydepositionGas(t)
```
"""
function DryDepositionGas(; name = :DryDepositionGas)
    rain = false
    dew = false
    params = @parameters begin
        season = Int(wesleyMidsummer), [description = "Index for season from Wesley (1989)"]
        landuse = Int(wesleyUrban), [description = "Index for land-use from Wesley (1989)"]
        z = 60, [unit = u"m", description = "Height from the ground to the mid-point of level 1"]
        del_P = 1520, [unit = u"Pa", description = "Pressure thinkness of level 1"]
        z₀ = 0.04, [unit = u"m", description = "Roughness length"]
        u_star = 0.44, [unit = u"m/s", description = "Friction velocity"]
        L = 0, [unit = u"m", description = "Monin-Obukhov length"]
        ρA = 1.2, [unit = u"kg*m^-3", description = "Air density"]
        G = 300, [unit = u"W*m^-2", description = "Solar irradiation"]
        Ts = 298, [unit = u"K", description = "Surface air temperature"]
        θ = 0, [description = "Slope of the local terrain, in unit radians"]
        lev = 1, [description = "Level of the atmospheric layer"]
    end

    depvel = @variables begin
        v_NO(t), [unit = u"m/s", description = "NO dry deposition velocity"]
        v_Ald(t),
            [
                unit = u"m/s", description = "Acetaldehyde (aldehyde class) dry deposition velocity",
            ]
        v_HCHO(t), [unit = u"m/s", description = "Formaldehyde dry deposition velocity"]
        v_OP(t),
            [
                unit = u"m/s",
                description = "Methyl hydroperoxide (organic peroxide class) dry deposition velocity",
            ]
        v_PAA(t),
            [unit = u"m/s", description = "Peroxyacetyl nitrate dry deposition velocity"]
        v_ORA(t),
            [
                unit = u"m/s", description = "Formic acid (organic acid class) dry deposition velocity",
            ]
        v_NH3(t), [unit = u"m/s", description = "NH3 dry deposition velocity"]
        v_HNO2(t), [unit = u"m/s", description = "Nitrous acid dry deposition velocity"]
        v_ACET(t), [unit = u"m/s", description = "Acetone dry deposition velocity"]
        v_ACTA(t), [unit = u"m/s", description = "Acetic acid dry deposition velocity"]
        v_ALD2(t), [unit = u"m/s", description = "Acetaldehyde dry deposition velocity"]
        v_AROMP4(t),
            [
                unit = u"m/s", description = "Generic C4 product of aromatics dry deposition velocity",
            ]
        v_AROMP5(t),
            [unit = u"m/s", description = "C5 unsaturated dicarbonyl dry deposition velocity"]
        v_ATOOH(t), [unit = u"m/s", description = "ATO2 peroxide dry deposition velocity"]
        v_BALD(t), [unit = u"m/s", description = "Benzaldehyde dry deposition velocity"]
        v_BENZP(t),
            [unit = u"m/s", description = "Phenyl hydroperoxide dry deposition velocity"]
        v_Br2(t), [unit = u"m/s", description = "Molecular Bromine dry deposition velocity"]
        v_BrCl(t), [unit = u"m/s", description = "Bromine chloride dry deposition velocity"]
        v_BrNO3(t), [unit = u"m/s", description = "Bromine nitrate dry deposition velocity"]
        v_BZCO3H(t),
            [unit = u"m/s", description = "Perbenzoic acid dry deposition velocity"]
        v_BZPAN(t),
            [unit = u"m/s", description = "Peroxybenzoylnitrate dry deposition velocity"]
        v_CH2O(t), [unit = u"m/s", description = "Formaldehyde dry deposition velocity"]
        v_Cl2(t),
            [unit = u"m/s", description = "Molecular chlorine dry deposition velocity"]
        v_ClNO2(t), [unit = u"m/s", description = "Nitryl chloride dry deposition velocity"]
        v_ClNO3(t),
            [unit = u"m/s", description = "Chlorine nitrate dry deposition velocity"]
        v_ClO(t), [unit = u"m/s", description = "Chlorine monoxide dry deposition velocity"]
        v_ClOO(t), [unit = u"m/s", description = "Chlorine dioxide dry deposition velocity"]
        v_CSL(t), [unit = u"m/s", description = "Cresols dry deposition velocity"]
        v_EOH(t), [unit = u"m/s", description = "Ethanol dry deposition velocity"]
        v_ETHLN(t), [unit = u"m/s", description = "Ethanol nitrate dry deposition velocity"]
        v_ETHN(t),
            [unit = u"m/s", description = "hydroxy-nitrooxy-ethane dry deposition velocity"]
        v_ETHP(t),
            [unit = u"m/s", description = "hydroxy-hydroperoxy-ethane dry deposition velocity"]
        v_ETNO3(t), [unit = u"m/s", description = "Ethyl nitrate dry deposition velocity"]
        v_ETP(t),
            [unit = u"m/s", description = "Ethylhydroperoxide dry deposition velocity"]
        v_GLYC(t), [unit = u"m/s", description = "Glycoaldehyde dry deposition velocity"]
        v_GLYX(t), [unit = u"m/s", description = "Glyoxal dry deposition velocity"]
        v_H2O2(t),
            [unit = u"m/s", description = "Hydrogen peroxide dry deposition velocity"]
        v_HAC(t), [unit = u"m/s", description = "Hydroxyacetone dry deposition velocity"]
        v_HBr(t), [unit = u"m/s", description = "Hypobromic acid dry deposition velocity"]
        v_HC5A(t),
            [
                unit = u"m/s", description = "isoprene-4,1-hydroxyaldehyde dry deposition velocity",
            ]
        v_HCl(t), [unit = u"m/s", description = "Hydrochloric acid dry deposition velocity"]
        v_HCOOH(t), [unit = u"m/s", description = "Formic acid dry deposition velocity"]
        v_HI(t), [unit = u"m/s", description = "Hydrogen iodide dry deposition velocity"]
        v_HMHP(t),
            [unit = u"m/s", description = "Hydroxymethyl hydroperoxide dry deposition velocity"]
        v_HMML(t),
            [
                unit = u"m/s", description = "hydroxymethyl-methyl-a-lactone dry deposition velocity",
            ]
        v_HNO3(t), [unit = u"m/s", description = "Nitric acid dry deposition velocity"]
        v_HOBr(t), [unit = u"m/s", description = "Hypobromous acid dry deposition velocity"]
        v_HOCl(t),
            [unit = u"m/s", description = "Hypochlorous acid dry deposition velocity"]
        v_HOI(t), [unit = u"m/s", description = "Hypoiodous acid dry deposition velocity"]
        v_HONIT(t),
            [
                unit = u"m/s",
                description = "2nd gen monoterpene organic nitrate dry deposition velocity",
            ]
        v_HPALD1(t),
            [
                unit = u"m/s", description = "d-4,1-C5-hydroperoxyaldehyde dry deposition velocity",
            ]
        v_HPALD2(t),
            [
                unit = u"m/s", description = "d-1,4-C5-hydroperoxyaldehyde dry deposition velocity",
            ]
        v_HPALD3(t),
            [
                unit = u"m/s", description = "b-2,1-C5-hydroperoxyaldehyde dry deposition velocity",
            ]
        v_HPALD4(t),
            [
                unit = u"m/s", description = "b-3,4-C5-hydroperoxyaldehyde dry deposition velocity",
            ]
        v_HPETHNL(t),
            [unit = u"m/s", description = "Hydroperoxy ethanal dry deposition velocity"]
        v_I2(t), [unit = u"m/s", description = "Molecular iodine dry deposition velocity"]
        v_I2O2(t), [unit = u"m/s", description = "Diiodine dioxide dry deposition velocity"]
        v_I2O3(t),
            [
                unit = u"m/s",
                description = "Diiodine trioxide dry deposition velocity",
            ]
        v_I2O4(t),
            [unit = u"m/s", description = "Diiodine tetraoxide dry deposition velocity"]
        v_IBr(t),
            [unit = u"m/s", description = "Iodine monobromide dry deposition velocity"]
        v_ICHE(t),
            [
                unit = u"m/s",
                description = "Isoprene hydroxy-carbonyl-epoxides dry deposition velocity",
            ]
        v_ICl(t),
            [unit = u"m/s", description = "Iodine monochloride dry deposition velocity"]
        v_ICN(t),
            [
                unit = u"m/s", description = "Lumped isoprene carbonyl-nitrates dry deposition velocity",
            ]
        v_ICPDH(t),
            [
                unit = u"m/s",
                description = "Isoprene dihydroxy hydroperoxycarbonyl dry deposition velocity",
            ]
        v_IDC(t),
            [unit = u"m/s", description = "Lumped isoprene dicarbonyls dry deposition velocity"]
        v_IDCHP(t),
            [
                unit = u"m/s",
                description = "Isoprene dicarbonyl hydroxy dihydroperoxide dry deposition velocity",
            ]
        v_IDHDP(t),
            [
                unit = u"m/s",
                description = "Isoprene dihydroxy dihydroperoxide dry deposition velocity",
            ]
        v_IDHPE(t),
            [
                unit = u"m/s",
                description = "Isoprene dihydroxy hydroperoxy epoxide dry deposition velocity",
            ]
        v_IDN(t),
            [unit = u"m/s", description = "Lumped isoprene dinitrates dry deposition velocity"]
        v_IEPOXA(t),
            [
                unit = u"m/s", description = "trans-Beta isoprene epoxydiol dry deposition velocity",
            ]
        v_IEPOXB(t),
            [unit = u"m/s", description = "cis-Beta isoprene epoxydiol dry deposition velocity"]
        v_IEPOXD(t),
            [unit = u"m/s", description = "Delta isoprene epoxydiol dry deposition velocity"]
        v_IHN1(t),
            [
                unit = u"m/s", description = "Isoprene-d-4,1-hydroxynitrate dry deposition velocity",
            ]
        v_IHN2(t),
            [
                unit = u"m/s", description = "Isoprene-b-1,2-hydroxynitrate dry deposition velocity",
            ]
        v_IHN3(t),
            [
                unit = u"m/s", description = "Isoprene-b-4,3-hydroxynitrate dry deposition velocity",
            ]
        v_IHN4(t),
            [
                unit = u"m/s", description = "Isoprene-d-4,1-hydroxynitrate dry deposition velocity",
            ]
        v_INPB(t),
            [
                unit = u"m/s",
                description = "Lumped b-hydroperoxy isoprene nitrates dry deposition velocity",
            ]
        v_INPD(t),
            [
                unit = u"m/s",
                description = "Lumped d-hydroperoxy isoprene nitrates dry deposition velocity",
            ]
        v_IONO(t), [unit = u"m/s", description = "Nitryl iodide dry deposition velocity"]
        v_IONO2(t), [unit = u"m/s", description = "Iodine nitrate dry deposition velocity"]
        v_IPRNO3(t),
            [unit = u"m/s", description = "Isopropyl nitrate dry deposition velocity"]
        v_ITCN(t),
            [
                unit = u"m/s",
                description = "lumped isoprene tetrafunctional carbonylnitrates dry deposition velocity",
            ]
        v_ITHN(t),
            [
                unit = u"m/s",
                description = "Lumped isoprene tetrafunctional hydroxynitrates dry deposition velocity",
            ]
        v_LIMO(t), [unit = u"m/s", description = "Limonene dry deposition velocity"]
        v_LVOC(t),
            [
                unit = u"m/s",
                description = "Gas-phase low-volatility non-IEPOX product of RIP ox dry deposition velocity",
            ]
        v_MACR(t), [unit = u"m/s", description = "Methacrolein dry deposition velocity"]
        v_MACR1OOH(t),
            [unit = u"m/s", description = "Peracid from MACR dry deposition velocity"]
        v_MAP(t), [unit = u"m/s", description = "Peroxyacetic acid dry deposition velocity"]
        v_MCRDH(t),
            [unit = u"m/s", description = "Dihydroxy-methacrolein dry deposition velocity"]
        v_MCRENOL(t),
            [unit = u"m/s", description = "Lumped enols from MVK/MACR dry deposition velocity"]
        v_MCRHN(t),
            [unit = u"m/s", description = "Nitrate from MACR dry deposition velocity"]
        v_MCRHNB(t),
            [unit = u"m/s", description = "Nitrate from MACR dry deposition velocity"]
        v_MCRHP(t),
            [
                unit = u"m/s", description = "Hydroxy-hydroperoxy-methacrolein dry deposition velocity",
            ]
        v_MCT(t),
            [
                unit = u"m/s", description = "Catechols and methyl catechols dry deposition velocity",
            ]
        v_MENO3(t), [unit = u"m/s", description = "Methyl nitrate dry deposition velocity"]
        v_MGLY(t), [unit = u"m/s", description = "Methylglyoxal dry deposition velocity"]
        v_MOH(t), [unit = u"m/s", description = "Methanol dry deposition velocity"]
        v_MONITS(t),
            [
                unit = u"m/s",
                description = "Saturated 1st gen monoterpene organic nitrate dry deposition velocity",
            ]
        v_MONITU(t),
            [
                unit = u"m/s",
                description = "Unsaturated 1st gen monoterpene organic nitrate dry deposition velocity",
            ]
        v_MPAN(t),
            [
                unit = u"m/s", description = "Peroxymethacroyl nitrate (PMN) dry deposition velocity",
            ]
        v_MTPA(t),
            [
                unit = u"m/s",
                description = "a-pinene, b-pinene, sabinene, carene dry deposition velocity",
            ]
        v_MTPO(t),
            [
                unit = u"m/s",
                description = "Terpinene, terpinolene, myrcene, ocimene, other monoterpenes dry deposition velocity",
            ]
        v_MVK(t),
            [unit = u"m/s", description = "Methyl vinyl ketone dry deposition velocity"]
        v_MVKDH(t), [unit = u"m/s", description = "dihydroxy-MVK dry deposition velocity"]
        v_MVKHC(t),
            [unit = u"m/s", description = "MVK hydroxy-carbonyl dry deposition velocity"]
        v_MVKHCB(t),
            [unit = u"m/s", description = "MVK hydroxy-carbonyl dry deposition velocity"]
        v_MVKHP(t),
            [unit = u"m/s", description = "MVK hydroxy-hydroperoxide dry deposition velocity"]
        v_MVKN(t), [unit = u"m/s", description = "Nitrate from MVK dry deposition velocity"]
        v_MVKPC(t),
            [unit = u"m/s", description = "MVK hydroperoxy-carbonyl dry deposition velocity"]
        v_N2O5(t),
            [unit = u"m/s", description = "Dinitrogen pentoxide dry deposition velocity"]
        v_NO2(t), [unit = u"m/s", description = "Nitrogen dioxide dry deposition velocity"]
        v_NPHEN(t), [unit = u"m/s", description = "Nitrophenols dry deposition velocity"]
        v_NPRNO3(t),
            [unit = u"m/s", description = "n-propyl nitrate dry deposition velocity"]
        v_O3(t), [unit = u"m/s", description = "Ozone dry deposition velocity"]
        v_PAN(t),
            [unit = u"m/s", description = "Peroxyacetyl nitrate dry deposition velocity"]
        v_PHEN(t), [unit = u"m/s", description = "Phenol dry deposition velocity"]
        v_PP(t), [unit = u"m/s", description = "Peroxide from PO2 dry deposition velocity"]
        v_PPN(t),
            [
                unit = u"m/s", description = "Lumped peroxypropionyl nitrate dry deposition velocity",
            ]
        v_PROPNN(t),
            [unit = u"m/s", description = "Propanone nitrate dry deposition velocity"]
        v_PRPN(t),
            [unit = u"m/s", description = "Peroxide from PRN1 dry deposition velocity"]
        v_PYAC(t), [unit = u"m/s", description = "Pyruvic acid dry deposition velocity"]
        v_R4N2(t),
            [unit = u"m/s", description = "Lumped alkyl nitrate dry deposition velocity"]
        v_R4P(t),
            [unit = u"m/s", description = "Peroxide from R4O2 dry deposition velocity"]
        v_RA3P(t),
            [unit = u"m/s", description = "Peroxide from A3O2 dry deposition velocity"]
        v_RB3P(t),
            [unit = u"m/s", description = "Peroxide from B3O2 dry deposition velocity"]
        v_RIPA(t), [unit = u"m/s", description = "1,2-ISOPOOH dry deposition velocity"]
        v_RIPB(t), [unit = u"m/s", description = "4,3-ISOPOOH dry deposition velocity"]
        v_RIPC(t), [unit = u"m/s", description = "d-1,4-ISOPOOH dry deposition velocity"]
        v_RIPD(t), [unit = u"m/s", description = "d-4,1-ISOPOOH dry deposition velocity"]
        v_RP(t), [unit = u"m/s", description = "Peroxide from RCO3 dry deposition velocity"]
        v_SO2(t), [unit = u"m/s", description = "Sulfur dioxide dry deposition velocity"]
        v_RCOOH(t),
            [unit = u"m/s", description = "> C2 organic acids dry deposition velocity"]
    end

    deprate = @variables begin
        k_NO(t), [unit = u"1/s", description = "NO dry deposition rate"]
        k_Ald(t),
            [unit = u"1/s", description = "Acetaldehyde (aldehyde class) dry deposition rate"]
        k_HCHO(t), [unit = u"1/s", description = "Formaldehyde dry deposition rate"]
        k_OP(t),
            [
                unit = u"1/s",
                description = "Methyl hydroperoxide (organic peroxide class) dry deposition rate",
            ]
        k_PAA(t), [unit = u"1/s", description = "Peroxyacetyl nitrate dry deposition rate"]
        k_ORA(t),
            [
                unit = u"1/s", description = "Formic acid (organic acid class) dry deposition rate",
            ]
        k_NH3(t), [unit = u"1/s", description = "NH3 dry deposition rate"]
        k_HNO2(t), [unit = u"1/s", description = "Nitrous acid dry deposition rate"]
        k_ACET(t), [unit = u"1/s", description = "Acetone dry deposition rate"]
        k_ACTA(t), [unit = u"1/s", description = "Acetic acid dry deposition rate"]
        k_ALD2(t), [unit = u"1/s", description = "Acetaldehyde dry deposition rate"]
        k_AROMP4(t),
            [unit = u"1/s", description = "Generic C4 product of aromatics dry deposition rate"]
        k_AROMP5(t),
            [unit = u"1/s", description = "C5 unsaturated dicarbonyl dry deposition rate"]
        k_ATOOH(t), [unit = u"1/s", description = "ATO2 peroxide dry deposition rate"]
        k_BALD(t), [unit = u"1/s", description = "Benzaldehyde dry deposition rate"]
        k_BENZP(t),
            [unit = u"1/s", description = "Phenyl hydroperoxide dry deposition rate"]
        k_Br2(t), [unit = u"1/s", description = "Molecular Bromine dry deposition rate"]
        k_BrCl(t), [unit = u"1/s", description = "Bromine chloride dry deposition rate"]
        k_BrNO3(t), [unit = u"1/s", description = "Bromine nitrate dry deposition rate"]
        k_BZCO3H(t), [unit = u"1/s", description = "Perbenzoic acid dry deposition rate"]
        k_BZPAN(t),
            [unit = u"1/s", description = "Peroxybenzoylnitrate dry deposition rate"]
        k_CH2O(t), [unit = u"1/s", description = "Formaldehyde dry deposition rate"]
        k_Cl2(t), [unit = u"1/s", description = "Molecular chlorine dry deposition rate"]
        k_ClNO2(t), [unit = u"1/s", description = "Nitryl chloride dry deposition rate"]
        k_ClNO3(t), [unit = u"1/s", description = "Chlorine nitrate dry deposition rate"]
        k_ClO(t), [unit = u"1/s", description = "Chlorine monoxide dry deposition rate"]
        k_ClOO(t), [unit = u"1/s", description = "Chlorine dioxide dry deposition rate"]
        k_CSL(t), [unit = u"1/s", description = "Cresols dry deposition rate"]
        k_EOH(t), [unit = u"1/s", description = "Ethanol dry deposition rate"]
        k_ETHLN(t), [unit = u"1/s", description = "Ethanol nitrate dry deposition rate"]
        k_ETHN(t),
            [unit = u"1/s", description = "hydroxy-nitrooxy-ethane dry deposition rate"]
        k_ETHP(t),
            [unit = u"1/s", description = "hydroxy-hydroperoxy-ethane dry deposition rate"]
        k_ETNO3(t), [unit = u"1/s", description = "Ethyl nitrate dry deposition rate"]
        k_ETP(t), [unit = u"1/s", description = "Ethylhydroperoxide dry deposition rate"]
        k_GLYC(t), [unit = u"1/s", description = "Glycoaldehyde dry deposition rate"]
        k_GLYX(t), [unit = u"1/s", description = "Glyoxal dry deposition rate"]
        k_H2O2(t), [unit = u"1/s", description = "Hydrogen peroxide dry deposition rate"]
        k_HAC(t), [unit = u"1/s", description = "Hydroxyacetone dry deposition rate"]
        k_HBr(t), [unit = u"1/s", description = "Hypobromic acid dry deposition rate"]
        k_HC5A(t),
            [unit = u"1/s", description = "isoprene-4,1-hydroxyaldehyde dry deposition rate"]
        k_HCl(t), [unit = u"1/s", description = "Hydrochloric acid dry deposition rate"]
        k_HCOOH(t), [unit = u"1/s", description = "Formic acid dry deposition rate"]
        k_HI(t), [unit = u"1/s", description = "Hydrogen iodide dry deposition rate"]
        k_HMHP(t),
            [unit = u"1/s", description = "Hydroxymethyl hydroperoxide dry deposition rate"]
        k_HMML(t),
            [unit = u"1/s", description = "hydroxymethyl-methyl-a-lactone dry deposition rate"]
        k_HNO3(t), [unit = u"1/s", description = "Nitric acid dry deposition rate"]
        k_HOBr(t), [unit = u"1/s", description = "Hypobromous acid dry deposition rate"]
        k_HOCl(t), [unit = u"1/s", description = "Hypochlorous acid dry deposition rate"]
        k_HOI(t), [unit = u"1/s", description = "Hypoiodous acid dry deposition rate"]
        k_HONIT(t),
            [
                unit = u"1/s", description = "2nd gen monoterpene organic nitrate dry deposition rate",
            ]
        k_HPALD1(t),
            [unit = u"1/s", description = "d-4,1-C5-hydroperoxyaldehyde dry deposition rate"]
        k_HPALD2(t),
            [unit = u"1/s", description = "d-1,4-C5-hydroperoxyaldehyde dry deposition rate"]
        k_HPALD3(t),
            [unit = u"1/s", description = "b-2,1-C5-hydroperoxyaldehyde dry deposition rate"]
        k_HPALD4(t),
            [unit = u"1/s", description = "b-3,4-C5-hydroperoxyaldehyde dry deposition rate"]
        k_HPETHNL(t),
            [unit = u"1/s", description = "Hydroperoxy ethanal dry deposition rate"]
        k_I2(t), [unit = u"1/s", description = "Molecular iodine dry deposition rate"]
        k_I2O2(t), [unit = u"1/s", description = "Diiodine dioxide dry deposition rate"]
        k_I2O3(t), [unit = u"1/s", description = "Diiodine trioxide dry deposition rate"]
        k_I2O4(t), [unit = u"1/s", description = "Diiodine tetraoxide dry deposition rate"]
        k_IBr(t), [unit = u"1/s", description = "Iodine monobromide dry deposition rate"]
        k_ICHE(t),
            [
                unit = u"1/s", description = "Isoprene hydroxy-carbonyl-epoxides dry deposition rate",
            ]
        k_ICl(t), [unit = u"1/s", description = "Iodine monochloride dry deposition rate"]
        k_ICN(t),
            [
                unit = u"1/s", description = "Lumped isoprene carbonyl-nitrates dry deposition rate",
            ]
        k_ICPDH(t),
            [
                unit = u"1/s",
                description = "Isoprene dihydroxy hydroperoxycarbonyl dry deposition rate",
            ]
        k_IDC(t),
            [unit = u"1/s", description = "Lumped isoprene dicarbonyls dry deposition rate"]
        k_IDCHP(t),
            [
                unit = u"1/s",
                description = "Isoprene dicarbonyl hydroxy dihydroperoxide dry deposition rate",
            ]
        k_IDHDP(t),
            [
                unit = u"1/s", description = "Isoprene dihydroxy dihydroperoxide dry deposition rate",
            ]
        k_IDHPE(t),
            [
                unit = u"1/s",
                description = "Isoprene dihydroxy hydroperoxy epoxide dry deposition rate",
            ]
        k_IDN(t),
            [unit = u"1/s", description = "Lumped isoprene dinitrates dry deposition rate"]
        k_IEPOXA(t),
            [unit = u"1/s", description = "trans-Beta isoprene epoxydiol dry deposition rate"]
        k_IEPOXB(t),
            [unit = u"1/s", description = "cis-Beta isoprene epoxydiol dry deposition rate"]
        k_IEPOXD(t),
            [unit = u"1/s", description = "Delta isoprene epoxydiol dry deposition rate"]
        k_IHN1(t),
            [unit = u"1/s", description = "Isoprene-d-4,1-hydroxynitrate dry deposition rate"]
        k_IHN2(t),
            [unit = u"1/s", description = "Isoprene-b-1,2-hydroxynitrate dry deposition rate"]
        k_IHN3(t),
            [unit = u"1/s", description = "Isoprene-b-4,3-hydroxynitrate dry deposition rate"]
        k_IHN4(t),
            [unit = u"1/s", description = "Isoprene-d-4,1-hydroxynitrate dry deposition rate"]
        k_INPB(t),
            [
                unit = u"1/s",
                description = "Lumped b-hydroperoxy isoprene nitrates dry deposition rate",
            ]
        k_INPD(t),
            [
                unit = u"1/s",
                description = "Lumped d-hydroperoxy isoprene nitrates dry deposition rate",
            ]
        k_IONO(t), [unit = u"1/s", description = "Nitryl iodide dry deposition rate"]
        k_IONO2(t), [unit = u"1/s", description = "Iodine nitrate dry deposition rate"]
        k_IPRNO3(t), [unit = u"1/s", description = "Isopropyl nitrate dry deposition rate"]
        k_ITCN(t),
            [
                unit = u"1/s",
                description = "lumped isoprene tetrafunctional carbonylnitrates dry deposition rate",
            ]
        k_ITHN(t),
            [
                unit = u"1/s",
                description = "Lumped isoprene tetrafunctional hydroxynitrates dry deposition rate",
            ]
        k_LIMO(t), [unit = u"1/s", description = "Limonene dry deposition rate"]
        k_LVOC(t),
            [
                unit = u"1/s",
                description = "Gas-phase low-volatility non-IEPOX product of RIP ox dry deposition rate",
            ]
        k_MACR(t), [unit = u"1/s", description = "Methacrolein dry deposition rate"]
        k_MACR1OOH(t),
            [unit = u"1/s", description = "Peracid from MACR dry deposition rate"]
        k_MAP(t), [unit = u"1/s", description = "Peroxyacetic acid dry deposition rate"]
        k_MCRDH(t),
            [unit = u"1/s", description = "Dihydroxy-methacrolein dry deposition rate"]
        k_MCRENOL(t),
            [unit = u"1/s", description = "Lumped enols from MVK/MACR dry deposition rate"]
        k_MCRHN(t), [unit = u"1/s", description = "Nitrate from MACR dry deposition rate"]
        k_MCRHNB(t), [unit = u"1/s", description = "Nitrate from MACR dry deposition rate"]
        k_MCRHP(t),
            [
                unit = u"1/s", description = "Hydroxy-hydroperoxy-methacrolein dry deposition rate",
            ]
        k_MCT(t),
            [unit = u"1/s", description = "Catechols and methyl catechols dry deposition rate"]
        k_MENO3(t), [unit = u"1/s", description = "Methyl nitrate dry deposition rate"]
        k_MGLY(t), [unit = u"1/s", description = "Methylglyoxal dry deposition rate"]
        k_MOH(t), [unit = u"1/s", description = "Methanol dry deposition rate"]
        k_MONITS(t),
            [
                unit = u"1/s",
                description = "Saturated 1st gen monoterpene organic nitrate dry deposition rate",
            ]
        k_MONITU(t),
            [
                unit = u"1/s",
                description = "Unsaturated 1st gen monoterpene organic nitrate dry deposition rate",
            ]
        k_MPAN(t),
            [unit = u"1/s", description = "Peroxymethacroyl nitrate (PMN) dry deposition rate"]
        k_MTPA(t),
            [
                unit = u"1/s", description = "a-pinene, b-pinene, sabinene, carene dry deposition rate",
            ]
        k_MTPO(t),
            [
                unit = u"1/s",
                description = "Terpinene, terpinolene, myrcene, ocimene, other monoterpenes dry deposition rate",
            ]
        k_MVK(t), [unit = u"1/s", description = "Methyl vinyl ketone dry deposition rate"]
        k_MVKDH(t), [unit = u"1/s", description = "dihydroxy-MVK dry deposition rate"]
        k_MVKHC(t),
            [unit = u"1/s", description = "MVK hydroxy-carbonyl dry deposition rate"]
        k_MVKHCB(t),
            [unit = u"1/s", description = "MVK hydroxy-carbonyl dry deposition rate"]
        k_MVKHP(t),
            [unit = u"1/s", description = "MVK hydroxy-hydroperoxide dry deposition rate"]
        k_MVKN(t), [unit = u"1/s", description = "Nitrate from MVK dry deposition rate"]
        k_MVKPC(t),
            [unit = u"1/s", description = "MVK hydroperoxy-carbonyl dry deposition rate"]
        k_N2O5(t), [unit = u"1/s", description = "Dinitrogen pentoxide dry deposition rate"]
        k_NO2(t), [unit = u"1/s", description = "Nitrogen dioxide dry deposition rate"]
        k_NPHEN(t), [unit = u"1/s", description = "Nitrophenols dry deposition rate"]
        k_NPRNO3(t), [unit = u"1/s", description = "n-propyl nitrate dry deposition rate"]
        k_O3(t), [unit = u"1/s", description = "Ozone dry deposition rate"]
        k_PAN(t), [unit = u"1/s", description = "Peroxyacetyl nitrate dry deposition rate"]
        k_PHEN(t), [unit = u"1/s", description = "Phenol dry deposition rate"]
        k_PP(t), [unit = u"1/s", description = "Peroxide from PO2 dry deposition rate"]
        k_PPN(t),
            [unit = u"1/s", description = "Lumped peroxypropionyl nitrate dry deposition rate"]
        k_PROPNN(t), [unit = u"1/s", description = "Propanone nitrate dry deposition rate"]
        k_PRPN(t), [unit = u"1/s", description = "Peroxide from PRN1 dry deposition rate"]
        k_PYAC(t), [unit = u"1/s", description = "Pyruvic acid dry deposition rate"]
        k_R4N2(t), [unit = u"1/s", description = "Lumped alkyl nitrate dry deposition rate"]
        k_R4P(t), [unit = u"1/s", description = "Peroxide from R4O2 dry deposition rate"]
        k_RA3P(t), [unit = u"1/s", description = "Peroxide from A3O2 dry deposition rate"]
        k_RB3P(t), [unit = u"1/s", description = "Peroxide from B3O2 dry deposition rate"]
        k_RIPA(t), [unit = u"1/s", description = "1,2-ISOPOOH dry deposition rate"]
        k_RIPB(t), [unit = u"1/s", description = "4,3-ISOPOOH dry deposition rate"]
        k_RIPC(t), [unit = u"1/s", description = "d-1,4-ISOPOOH dry deposition rate"]
        k_RIPD(t), [unit = u"1/s", description = "d-4,1-ISOPOOH dry deposition rate"]
        k_RP(t), [unit = u"1/s", description = "Peroxide from RCO3 dry deposition rate"]
        k_SO2(t), [unit = u"1/s", description = "Sulfur dioxide dry deposition rate"]
        k_RCOOH(t), [unit = u"1/s", description = "> C2 organic acids dry deposition rate"]
    end

    datas = [
        NoData,
        AldData,
        HchoData,
        OpData,
        PaaData,
        OraData,
        Nh3Data,
        Hno2Data,
        ACETData,
        ACTAData,
        ALD2Data,
        AROMP4Data,
        AROMP5Data,
        ATOOHData,
        BALDData,
        BENZPData,
        Br2Data,
        BrClData,
        BrNO3Data,
        BZCO3HData,
        BZPANData,
        CH2OData,
        Cl2Data,
        ClNO2Data,
        ClNO3Data,
        ClOData,
        ClOOData,
        CSLData,
        EOHData,
        ETHLNData,
        ETHNData,
        ETHPData,
        ETNO3Data,
        ETPData,
        GLYCData,
        GLYXData,
        H2O2Data,
        HACData,
        HBrData,
        HC5AData,
        HClData,
        HCOOHData,
        HIData,
        HMHPData,
        HMMLData,
        HNO3Data,
        HOBrData,
        HOClData,
        HOIData,
        HONITData,
        HPALD1Data,
        HPALD2Data,
        HPALD3Data,
        HPALD4Data,
        HPETHNLData,
        I2Data,
        I2O2Data,
        I2O3Data,
        I2O4Data,
        IBrData,
        ICHEData,
        IClData,
        ICNData,
        ICPDHData,
        IDCData,
        IDCHPData,
        IDHDPData,
        IDHPEData,
        IDNData,
        IEPOXAData,
        IEPOXBData,
        IEPOXDData,
        IHN1Data,
        IHN2Data,
        IHN3Data,
        IHN4Data,
        INPBData,
        INPDData,
        IONOData,
        IONO2Data,
        IPRNO3Data,
        ITCNData,
        ITHNData,
        LIMOData,
        LVOCData,
        MACRData,
        MACR1OOHData,
        MAPData,
        MCRDHData,
        MCRENOLData,
        MCRHNData,
        MCRHNBData,
        MCRHPData,
        MCTData,
        MENO3Data,
        MGLYData,
        MOHData,
        MONITSData,
        MONITUData,
        MPANData,
        MTPAData,
        MTPOData,
        MVKData,
        MVKDHData,
        MVKHCData,
        MVKHCBData,
        MVKHPData,
        MVKNData,
        MVKPCData,
        N2O5Data,
        NO2Data,
        NPHENData,
        NPRNO3Data,
        O3Data,
        PANData,
        PHENData,
        PPData,
        PPNData,
        PROPNNData,
        PRPNData,
        PYACData,
        R4N2Data,
        R4PData,
        RA3PData,
        RB3PData,
        RIPAData,
        RIPBData,
        RIPCData,
        RIPDData,
        RPData,
        SO2Data,
        RCOOHData,
    ]

    isSO2 = repeat([false], size(datas)[1])
    isSO2[131] = true
    isO3 = repeat([false], size(datas)[1])
    isO3[114] = true
    eqs = [
        depvel .~ DryDepGas.(
            lev, z, z₀, u_star, L, ρA, datas, G, Ts, θ,
            season, landuse, rain, dew, isSO2, isO3
        );
        deprate .~ depvel * g * ρA / del_P
    ]

    return System(
        eqs,
        t,
        [depvel; deprate],
        [
            params;
            [G_unitless, T_unitless, unit_dH2O, Rc_unit, unit_T, unit_m, unit_convert_mu, κ, g, k, M_air, R, unit_v]
        ];
        name = name,
        metadata = Dict(CoupleType => DryDepositionGasCoupler)
    )
end

struct DryDepositionAerosolCoupler
    sys::Any
end

"""
Aerosol dry deposition based on Seinfeld and Pandis (2006) equation 19.7.
"""
function DryDepositionAerosol(; name = :DryDepositionAerosol)
    params = @parameters begin
        SeinfeldSeason::Int = Int(seinfeldMidsummer),
            [description = "Index for season from Seinfeld and Pandis (2006)"]
        WesleySeason::Int = Int(wesleyMidsummer),
            [description = "Index for season from Wesley (1989)"]
        SeinfeldLandUse::Int = Int(seinfeldGrass),
            [description = "Index for land use from Seinfeld and Pandis (2006)"]
        WesleyLandUse::Int = Int(wesleyRangeAg),
            [description = "Index for land use from Wesley (1989)"]
        z = 50, [unit = u"m", description = "Top of the surface layer"]
        z₀ = 0.04, [unit = u"m", description = "Roughness length"]
        u_star = 0.44, [unit = u"m/s", description = "Friction velocity"]
        L = 0, [unit = u"m", description = "Monin-Obukhov length"]
        ρA = 1.2, [unit = u"kg*m^-3", description = "Air density"]
        Ts = 298, [unit = u"K", description = "Surface air temperature"]
        lev = 1, [description = "Level of the atmospheric layer"]
        Dp = 0.8e-6, [unit = u"m", description = "Particle diameter"]
        P = 101325, [unit = u"Pa", description = "Pressure"]
        ρParticle = 1000.0, [unit = u"kg*m^-3", description = "Particle density"]
    end

    @variables begin
        v(t), [unit = u"m/s", description = "Particle dry deposition velocity"]
        k(t), [unit = u"1/s", description = "Particle dry deposition rate"]
    end
    eqs = [
        v ~ DryDepParticle(
            lev, z, z₀, u_star, L, Dp, Ts, P, ρParticle, ρA,
            SeinfeldSeason, WesleySeason, SeinfeldLandUse, WesleyLandUse
        ),
        k ~ v / z,
    ]

    return System(
        eqs,
        t,
        [v; k],
        params;
        name = name,
        metadata = Dict(CoupleType => DryDepositionAerosolCoupler)
    )
end
