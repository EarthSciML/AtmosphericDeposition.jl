export ra, mu, mfp, cc, vs, dParticle, dH2O, sc, stSmooth, stVeg, RbGas, z₀_table, A_table, α_table, γ_table, RbParticle, DryDepGas, DryDepParticle

g = 9.81u"m*s^-2" # gravitational acceleration [m/s2]
κ = 0.4 # von Karman constant
k = 1.3806488e-23u"m^2*kg*s^-2/K" # Boltzmann constant
M_air = 28.97e-3u"kg/mol" # molecular weight of air
R = 8.3144621u"J/K/mol" # Gas constant

"""
Function Ra calculates aerodynamic resistance to dry deposition 
where z is the top of the surface layer [m], z₀ is the roughness length [m], u_star is friction velocity [m/s], and L is Monin-Obukhov length [m]
Based on Seinfeld and Pandis (2006) [Seinfeld and Pandis (2006)] equation 19.13 & 19.14.
"""

function ra(z, z₀, u_star, L)
    if L == 0u"m"
        ζ = 0
        ζ₀= 0
    else
        ζ = z/L
        ζ₀= z₀/L
    end
    print(ζ)
    if 0 < ζ < 1
        rₐ = 1/(κ*u_star)*(log(z/z₀)+4.7*(ζ-ζ₀))
    elseif ζ == 0
        rₐ = 1/(κ*u_star)*log(z/z₀)
    elseif -1 < ζ < 0
        η₀=(1-15*ζ₀)^1/4
        η =(1-15*ζ)^1/4
        rₐ = 1/(κ*u_star)*[log(z/z₀)+log(((η₀^2+1)*(η₀+1)^2)/((η^2+1)*(η+1)^2))+2*(atan(η)-atan(η₀))]
    else
        print("wrong rₐ")
    end
    return rₐ
end

"""
Function mu calculates the dynamic viscosity of air [kg m-1 s-1] where T is temperature [K].
"""
function mu(T)
    return (1.8e-5*(T/298u"K")^0.85)u"kg/m/s"
end

"""
Function mfp calculates the mean free path of air [m] 
where T is temperature [K] P is pressure [Pa], and Mu is dynamic viscosity [kg/(m s)].
From Seinfeld and Pandis (2006) equation 9.6
"""
function mfp(T,P,μ)
    return 2*μ/(P*(8*M_air/(pi*R*T))^0.5)
end

"""
Function cc calculates the Cunnningham slip correction factor
where Dp is particle diameter [m], T is temperature [K], and P is pressure [Pa].
From Seinfeld and Pandis (2006) equation 9.34.
"""
function cc(Dₚ,T,P,μ)
    λ = mfp(T,P,μ)
    return 1+2*λ/Dₚ*(1.257+0.4*exp(-1.1*Dₚ/(2*λ)))
end

"""
Function vs calculates the terminal setting velocity of a
particle where Dp is particle diameter [m], ρₚ is particle density [kg/m3], Cc is the Cunningham slip correction factor, and μ is air dynamic viscosity [kg/(s m)]. 
From equation 9.42 in Seinfeld and Pandis (2006).
"""
function vs(Dₚ, ρₚ, Cc, μ)
    if Dₚ > 20.e-6u"m"
        print("Particle diameter ", Dₚ ," [m] is greater than 20um; Stokes settling no longer applies.")
    else
        return Dₚ^2*ρₚ*g*Cc/(18*μ)
    end
end

"""
Function dParticle calculates the brownian diffusivity of a particle [m2/s] using the Stokes-Einstein-Sutherland relation
(Seinfeld and Pandis eq. 9.73) 
where T is air temperature [K], P is pressure [Pa], Dp is particle diameter [m], and μ is air dynamic viscosity [kg/(s m)]
"""
function dParticle(T,P,Dₚ,Cc,μ)
    return k*T*Cc/(3*pi*μ*Dₚ)
end

"""
Function dH2O calculates molecular diffusivity of water vapor in air [m2/s] where T is temperature [K]
using a regression fit to data in Bolz and Tuve (1976) found here: http://www.cambridge.org/us/engineering/author/nellisandklein/downloads/examples/EXAMPLE_9.2-1.pdf
"""
function dH2O(T)
    T_unitless = T*1u"K^-1"
    return (-2.775e-6 + 4.479e-8*T_unitless + 1.656e-10*T_unitless^2)u"m^2/s"
end

"""
Function sc computes the dimensionless Schmidt number,
where μ is dynamic viscosity of air [kg/(s m)], ρ is air density [kg/m3], and D is the molecular diffusivity of the gas speciesof interest [m2/s]
"""
function sc(μ, ρ, D)
    return μ/(ρ*D)
end

"""
Function stSmooth computes the dimensionless Stokes number for dry deposition of particles on smooth surfaces or surfaces with bluff roughness elements, 
where vs is settling velocity [m/s], u_star is friction velocity [m/s], μ is dynamic viscosity of air [kg/(s m)], and ρ is air density [kg/m3], 
based on Seinfeld and Pandis (2006) equation 19.23.
"""
function stSmooth(vₛ, u_star, μ, ρ)
    return vₛ*u_star^2*ρ/(g*μ)
end

"""
Function stVeg computes the dimensionless Stokes number for dry deposition of particles on vegetated surfaces, 
where vs is settling velocity [m/s], u_star is friction velocity [m/s], and A is the characteristic collector radius [m],
based on Seinfeld and Pandis (2006) equation 19.24.
"""
function stVeg(vₛ, u_star, A)
    return vₛ*u_star/(g*A)
end

"""
Function RbGas calculates the quasi-laminar sublayer resistance to dry deposition for a gas species [s/m],
where Sc is the dimensionless Schmidt number and u_star is the friction velocity [m/s]. 
From Seinfeld and Pandis (2006) equation 19.17.
"""
function RbGas(Sc, u_star)
    return 5*Sc^(2/3)/u_star
end

"""
Values for the characteristic radii of collectors [m]
where the columns are land use categories and the rows are seasonal categories.
Land-use categories (LUCs)
1. Evergreen–needleleaf trees
2. Deciduous broadleaf trees
3. Grass
4. Desert
5. Shrubs and interrupted woodlands
Seasonal categories (SC)
1. Midsummer with lush vegetation
2. Autumn with cropland not harvested
3. Late autumn after frost, no snow
4. Winter, snow on ground
5. Transitional
given in Seinfeld and Pandis Table 19.2
"""
z₀_table = SA_F32[
    0.8 1.05 0.1 0.04 0.1
    0.9 1.05 0.1 0.04 0.1
    0.9 0.95 0.05 0.04 0.1
    0.9 0.55 0.02 0.04 0.1
    0.8 0.75 0.05 0.04 0.1
] # unit:[m]

A_table = SA_F32[
    2.0 5.0 2.0 Inf 10.0
    2.0 5.0 2.0 Inf 10.0
    2.0 10.0 5.0 Inf 10.0
    2.0 10.0 2.0 Inf 10.0
    2.0 5.0 2.0 Inf 10.0
] # unit:[mm]

α_table = SA_F32[
    1.0 0.8 1.2 50.0 1.3
]

γ_table = SA_F32[
    0.56 0.56 0.54 0.54 0.54
]

"""
Function RbParticle calculates the quasi-laminar sublayer resistance to dry deposition for a particles [s/m],
where Sc is the dimensionless Schmidt number, u_star is the friction velocity [m/s], St is the dimensionless Stokes number, 
Dp is particle diameter [m], and iSeason and iLandUse are season and land use indexes, respectively.
From Seinfeld and Pandis (2006) equation 19.27.
"""

function RbParticle(Sc, u_star, St, Dₚ, iSeason::Int, iLandUse::Int) 
    α = α_table[iLandUse]
    γ = γ_table[iLandUse]
    A = (A_table[iSeason,iLandUse]*10^(-3))u"m"
    R1 = exp(-St^0.5)
    term_1 = Sc^(-γ)
    term_2 = (St/(α+St))^2
    term_3 = 1/2*(Dₚ/A)^2
    return 1/(3*u_star*(term_1+term_2+term_3)*R1)
end

"""
Function DryDepGas calculates dry deposition velocity [m/s] for a gas species,
where z is the height of the surface layer [m], zo is roughness length [m], u_star is friction velocity [m/s], 
L is Monin-Obukhov length [m], T is surface air temperature [K], ρA is air density [kg/m3]
gasData is data about the gas species for surface resistance calculations, G is solar
irradiation [W m-2], Θ is the slope of the local terrain [radians], iSeason and iLandUse are indexes for the season and land use,
dew and rain indicate whether there is dew or rain on the ground, and isSO2 and isO3 indicate whether the gas species of interest is either SO2 or O3, respectively. 
Based on Seinfeld and Pandis (2006) equation 19.2.
"""
function DryDepGas(z, z₀, u_star, L, ρA, gasData::GasData, G, Ts, θ, iSeason::Int, iLandUse::Int, rain::Bool, dew::Bool, isSO2::Bool, isO3::Bool) 
    Ra = ra(z, z₀, u_star, L)
    μ = mu(Ts)
    Dg = dH2O(Ts)/gasData.Dh2oPerDx #Diffusivity of gas of interest [m2/s]
    Sc = sc(μ,ρA, Dg)
    Rb = RbGas(Sc, u_star)
    Rc = SurfaceResistance(gasData, G/u"W*m^-2", (Ts/u"K"-273), θ, iSeason::Int, iLandUse::Int, rain::Bool, dew::Bool, isSO2::Bool, isO3::Bool)u"s/m"
    return 1/(Ra+Rb+Rc)
end

"""
Function DryDepParticle calculates particle dry deposition velocity [m/s]
where z is the height of the surface layer [m], zo is roughness length [m], u_star is friction velocity [m/s], L is Monin-Obukhov length [m], 
Dp is particle diameter [m], Ts is surface air temperature [K], P is pressure [Pa], ρParticle is particle density [kg/m3], ρAir is air density [kg/m3],
and iSeason and iLandUse are indexes for the season and land use.
Based on Seinfeld and Pandis (2006) equation 19.7.
"""
function DryDepParticle(z, z₀, u_star, L, Dp, Ts, P, ρParticle, ρA, iSeason::Int, iLandUse::Int)
    Ra = ra(z, z₀, u_star, L)
    μ = mu(Ts)
    Cc = cc(Dp, Ts, P, μ)
    Vs = vs(Dp, ρParticle, Cc, μ)
    if iLandUse == 4 # dessert
        St = stSmooth(Vs, u_star, μ, ρA)
    else
        St = stVeg(Vs, u_star, (A_table[iSeason,iLandUse]*10^(-3))u"m")
    end
    D = dParticle(Ts, P, Dp, Cc, μ)
    Sc = sc(μ, ρA, D)
    Rb = RbParticle(Sc, u_star, St, Dp, iSeason, iLandUse)
    return 1/(Ra+Rb+Ra*Rb*Vs)+Vs
end

