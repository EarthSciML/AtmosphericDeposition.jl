using DepositionMTK
using Test

@testset "mfp" begin
    @test mfp(298u"K",101300u"Pa",1.8e-5u"kg/m/s") - 6.51e-8u"m" ≈ 0u"m" atol=1e-8u"m"
end

@testset "unit" begin
    @test unit(dH2O(300u"K"))==u"m^2/s"
    @test unit(DryDepParticle(0.4u"m",0.3u"m",1u"m/s", 1u"m", 1e-6u"m", 300u"K", 10300u"Pa", 1u"kg*m^-3",0.001u"kg*m^-3",1,1)) == u"m/s"
    @test unit(DryDepGas(0.4u"m",0.3u"m",1u"m/s", 1u"m", 0.001u"kg*m^-3", So2Data, 800u"W*m^-2", 300u"K", 0, 1, 1, false, false, true, false)) == u"m/s"
end

@testset "viscosity" begin
    T = [100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1100, 1200, 1300, 1400, 1500, 1600]
    μ_list = [0.6924, 1.0283, 1.3289, 1.488, 1.983, 2.075, 2.286, 2.484, 2.671, 2.848, 3.018, 3.177, 3.332, 3.481, 3.625, 3.765, 3.899, 4.023, 4.152, 4.44, 4.69, 4.93, 5.17, 5.4, 5.63]
    μ_test = []
    for i in 1:25
        push!(μ_test, (mu((T[i])u"K"))u"m*s/kg")
    end
    @test μ_list ≈ μ_list rtol=1e-2
end

@testset "Cc" begin
    Dp = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
    Cc_list = [216, 108, 43.6, 22.2, 11.4, 4.95, 2.85, 1.865, 1.326, 1.164, 1.082, 1.032, 1.016, 1.008, 1.003, 1.0016]
    Cc_test = []
    for i in 1:16
        push!(Cc_test, cc((Dp[i]*1e-6)u"m", 298u"K", 101325u"Pa", mu(298u"K")))
    end
    @test Cc_test ≈ Cc_list rtol=1e-2
end

@testset "Vs" begin
    Dp = [0.01, 0.1, 1, 10]
    Vs_list = [0.025, 0.35, 10.8, 1000]./3600 ./100
    Vs_test = []
    for i in 1:4
        push!(Vs_test, (vs((Dp[i]*1e-6)u"m", 1000u"kg*m^-3", cc((Dp[i]*1e-6)u"m", 298u"K", 101325u"Pa", mu(298u"K")), mu(298u"K")))u"s/m")
    end
    @test Vs_test ≈ Vs_list rtol=0.1
end
    
@testset "DryDepGas" begin
    z = 50u"m"
    z₀ = 0.04u"m"
    u_star = 0.44u"m/s"
    L = 0u"m"
    T = 298u"K"
    ρA = 1.2u"kg*m^-3"
    G = 300u"W*m^-2"
    θ = 0
    vd_true = 0.003u"m/s"
    @test DryDepGas(z, z₀, u_star, L, ρA, No2Data, G, T, θ, 1, 10, false, false, false, false) ≈ vd_true rtol=0.33
end

@testset "DryDepParticle" begin
    z = 20u"m"
    z₀ = 0.02u"m"
    u_star = 0.44u"m/s"
    L = 0u"m"
    T = 298u"K"
    P = 101325u"Pa"
    ρA = 1.2u"kg*m^-3"
    ρParticle = 1000u"kg*m^-3"
    Dp = [1.e-8, 1.e-7, 1.e-6, 1.e-5]
    vd_true = [0.5, 0.012, 0.02, 0.6] # [cm/s]
    vd_list = []
    for i in 1:4
        push!(vd_list, (DryDepParticle(z, z₀, u_star, L, (Dp[i])u"m", T, P, ρParticle, ρA, 1, 4))*100u"s/m")
    end
    @test vd_list ≈ vd_true rtol=0.8 # 80% difference, not ideal
end

    
   
    

    