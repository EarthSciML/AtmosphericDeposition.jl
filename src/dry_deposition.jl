export DryDepositionGas, DryDepositionAerosol

@constants g=9.81 [unit = u"m*s^-2", description = "gravitational acceleration"]
@constants κ=0.4 [description = "von Karman constant"]
@constants k=1.3806488e-23 [unit = u"m^2*kg*s^-2/K", description = "Boltzmann constant"]
@constants M_air=28.97e-3 [unit = u"kg/mol", description = "molecular weight of air"]
@constants R=8.3144621 [unit = u"kg*m^2*s^−2*K^-1*mol^-1", description = "Gas constant"]

@constants unit_m=1 [unit = u"m"]
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
            2 * (atan(η) - atan(η₀))
        ][1],
        rₐ_1
    ) #the [1] is to pass the ModelingToolkit unit check
    return rₐ
end

@constants unit_T=1 [unit = u"K", description = "unit one for temperature"]
@constants unit_convert_mu=1 [unit = u"kg/m/s", description = "unit one for mu"]
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

@constants unit_v=1 [unit = u"m/s", description = "unit one for speed"]
"""
Function vs calculates the terminal setting velocity of a
particle where Dp is particle diameter [m], ρₚ is particle density [kg/m3], Cc is the Cunningham slip correction factor, and μ is air dynamic viscosity [kg/(s m)].
From equation 9.42 in Seinfeld and Pandis (2006).
"""
function vs(Dₚ, ρₚ, Cc, μ)
    ifelse((Dₚ > 20.e-6 * unit_m), 99999999 * unit_v, Dₚ^2 * ρₚ * g * Cc / (18 * μ))
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

@constants T_unitless=1 [unit = u"K^-1", description = "used to offset temperature unit"]
@constants unit_dH2O=1 [unit = u"m^2/s", description = "unit for molecular diffusivity"]
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
z₀_table = [0.8 1.05 0.1 0.04 0.1
            0.9 1.05 0.1 0.04 0.1
            0.9 0.95 0.05 0.04 0.1
            0.9 0.55 0.02 0.04 0.1
            0.8 0.75 0.05 0.04 0.1] # unit:[m]

function A_table(iSeason, iLandUse)
    SA_F32[2.0 5.0 2.0 Inf 10.0
           2.0 5.0 2.0 Inf 10.0
           2.0 10.0 5.0 Inf 10.0
           2.0 10.0 2.0 Inf 10.0
           2.0 5.0 2.0 Inf 10.0][
        iSeason,
        iLandUse
    ] * 1e-3 # unit:[mm]
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

@constants G_unitless=1 [
    unit = u"m^2/W",
    description = "used to offset the unit of irradiation"
]
@constants Rc_unit=1 [unit = u"s/m", description = "unit for surface resistance"]
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
function DryDepParticle(lev, z, z₀, u_star, L, Dp, Ts, P, ρParticle, ρA,
        iSeinfeldSeason, iWesleySeason, iSeinfeldLandUse, iWesleyLandUse)
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
    @constants v_zero=0 [unit = u"m/s", description = "zero velocity"]
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
    unit_v => 1
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
        season=Int(wesleyMidsummer), [description = "Index for season from Wesley (1989)"]
        landuse=Int(wesleyUrban), [description = "Index for land-use from Wesley (1989)"]
        z=60, [unit = u"m", description = "Height from the ground to the mid-point of level 1"]
        del_P=1520, [unit = u"Pa", description = "Pressure thinkness of level 1"]
        z₀=0.04, [unit = u"m", description = "Roughness length"]
        u_star=0.44, [unit = u"m/s", description = "Friction velocity"]
        L=0, [unit = u"m", description = "Monin-Obukhov length"]
        ρA=1.2, [unit = u"kg*m^-3", description = "Air density"]
        G=300, [unit = u"W*m^-2", description = "Solar irradiation"]
        Ts=298, [unit = u"K", description = "Surface air temperature"]
        θ=0, [description = "Slope of the local terrain, in unit radians"]
        lev=1, [description = "Level of the atmospheric layer"]
    end

    depvel = @variables begin
        v_SO2(t), [unit = u"m/s", description = "SO2 dry deposition velocity"]
        v_O3(t), [unit = u"m/s", description = "O3 dry deposition velocity"]
        v_NO2(t), [unit = u"m/s", description = "NO2 dry deposition velocity"]
        v_NO(t), [unit = u"m/s", description = "NO dry deposition velocity"]
        v_HNO3(t), [unit = u"m/s", description = "HNO3 dry deposition velocity"]
        v_H2O2(t), [unit = u"m/s", description = "H2O2 dry deposition velocity"]
        v_Ald(t), [unit = u"m/s",
            description = "Acetaldehyde (aldehyde class) dry deposition velocity"]
        v_HCHO(t), [unit = u"m/s", description = "Formaldehyde dry deposition velocity"]
        v_OP(t), [unit = u"m/s",
            description = "Methyl hydroperoxide (organic peroxide class) dry deposition velocity"]
        v_PAA(t), [unit = u"m/s", description = "Peroxyacetyl nitrate dry deposition velocity"]
        v_ORA(t), [unit = u"m/s",
            description = "Formic acid (organic acid class) dry deposition velocity"]
        v_NH3(t), [unit = u"m/s", description = "NH3 dry deposition velocity"]
        v_PAN(t), [unit = u"m/s",
            description = "Peroxyacetyl nitrate dry deposition velocity"]
        v_HNO2(t), [unit = u"m/s", description = "Nitrous acid dry deposition velocity"]
    end
    deprate = @variables begin
        k_SO2(t), [unit = u"1/s", description = "SO2 dry deposition rate"]
        k_O3(t), [unit = u"1/s", description = "O3 dry deposition rate"]
        k_NO2(t), [unit = u"1/s", description = "NO2 dry deposition rate"]
        k_NO(t), [unit = u"1/s", description = "NO dry deposition rate"]
        k_HNO3(t), [unit = u"1/s", description = "HNO3 dry deposition rate"]
        k_H2O2(t), [unit = u"1/s", description = "H2O2 dry deposition rate"]
        k_Ald(t), [unit = u"1/s",
            description = "Acetaldehyde (aldehyde class) dry deposition rate"]
        k_HCHO(t), [unit = u"1/s", description = "Formaldehyde dry deposition rate"]
        k_OP(t), [unit = u"1/s",
            description = "Methyl hydroperoxide (organic peroxide class) dry deposition rate"]
        k_PAA(t), [unit = u"1/s", description = "Peroxyacetyl nitrate dry deposition rate"]
        k_ORA(t), [unit = u"1/s",
            description = "Formic acid (organic acid class) dry deposition rate"]
        k_NH3(t), [unit = u"1/s", description = "NH3 dry deposition rate"]
        k_PAN(t), [unit = u"1/s", description = "Peroxyacetyl nitrate dry deposition rate"]
        k_HNO2(t), [unit = u"1/s", description = "Nitrous acid dry deposition rate"]
    end

    datas = [So2Data, O3Data, No2Data, NoData, Hno3Data, H2o2Data, AldData, HchoData,
        OpData, PaaData, OraData, Nh3Data, PanData, Hno2Data]

    isSO2 = repeat([false], 14)
    isSO2[1] = true
    isO3 = repeat([false], 14)
    isO3[2] = true

    eqs = [
        depvel .~ DryDepGas.(lev, z, z₀, u_star, L, ρA, datas, G, Ts, θ,
            season, landuse, rain, dew, isSO2, isO3);
        deprate .~ depvel*g*ρA / del_P]

    ODESystem(
        eqs,
        t,
        [depvel; deprate],
        params;
        name = name,
        metadata = Dict(:coupletype => DryDepositionGasCoupler)
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
        SeinfeldSeason::Int=Int(seinfeldMidsummer),
        [description = "Index for season from Seinfeld and Pandis (2006)"]
        WesleySeason::Int=Int(wesleyMidsummer),
        [description = "Index for season from Wesley (1989)"]
        SeinfeldLandUse::Int=Int(seinfeldGrass),
        [description = "Index for land use from Seinfeld and Pandis (2006)"]
        WesleyLandUse::Int=Int(wesleyRangeAg),
        [description = "Index for land use from Wesley (1989)"]
        z=50, [unit = u"m", description = "Top of the surface layer"]
        z₀=0.04, [unit = u"m", description = "Roughness length"]
        u_star=0.44, [unit = u"m/s", description = "Friction velocity"]
        L=0, [unit = u"m", description = "Monin-Obukhov length"]
        ρA=1.2, [unit = u"kg*m^-3", description = "Air density"]
        Ts=298, [unit = u"K", description = "Surface air temperature"]
        lev=1, [description = "Level of the atmospheric layer"]
        Dp=0.8e-6, [unit = u"m", description = "Particle diameter"]
        P=101325, [unit = u"Pa", description = "Pressure"]
        ρParticle=1.0, [unit = u"kg*m^-3", description = "Particle density"]
    end

    @variables begin
        v(t), [unit = u"m/s", description = "Particle dry deposition velocity"]
        k(t), [unit = u"1/s", description = "Particle dry deposition rate"]
    end
    eqs = [
        v ~ DryDepParticle(lev, z, z₀, u_star, L, Dp, Ts, P, ρParticle, ρA,
            SeinfeldSeason, WesleySeason, SeinfeldLandUse, WesleyLandUse),
        k ~ v / z
    ]

    ODESystem(
        eqs,
        t,
        [v; k],
        params;
        name = name,
        metadata = Dict(:coupletype => DryDepositionAerosolCoupler)
    )
end
