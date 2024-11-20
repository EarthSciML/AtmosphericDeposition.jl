using AtmosphericDeposition
using Test, DynamicQuantities, ModelingToolkit

begin
    @parameters T [unit = u"K"]
    @parameters P [unit = u"kg*m^-1*s^-2"] # 1Pa = 1kg/m/s^-2
    @parameters μ [unit = u"kg/m/s"]
    @parameters z [unit = u"m"]
    @parameters z₀ [unit = u"m"]
    @parameters u_star [unit = u"m/s"]
    @parameters L [unit = u"m"]
    @parameters Dp [unit = u"m"]
    @parameters ρParticle [unit = u"kg*m^-3"]
    @parameters ρA [unit = u"kg*m^-3"]
    @parameters G [unit = u"W*m^-2"]
    @parameters θ, iLandUse, iSeason
end

@testset "mfp" begin
    @test substitute(mfp(T, P, μ), Dict(T => 298, P => 101300, μ => 1.8e-5, AtmosphericDeposition.defaults...)) ≈ 6.512893276888993e-8
    @test ModelingToolkit.get_unit(mfp(T, P, μ)) == u"m"
end


@testset "unit" begin
    @test ModelingToolkit.get_unit(dH2O(T)) == u"m^2/s"
    @test ModelingToolkit.get_unit(DryDepParticle(z, z₀, u_star, L, Dp, T, P, ρParticle, ρA, 1, 1)) == u"m/s"
    @test ModelingToolkit.get_unit(DryDepGas(z, z₀, u_star, L, ρA, AtmosphericDeposition.So2Data, G, T, θ, iSeason, iLandUse, false, false, true, false)) == u"m/s"
end

@testset "viscosity" begin
    T_ = [275, 300, 325, 350, 375, 400]
    μ_list = [1.725, 1.846, 1.962, 2.075, 2.181, 2.286] .* 10^-5
    μ_test = []
    for i in 1:6
        push!(μ_test, substitute(mu(T), Dict(T => T_[i], AtmosphericDeposition.defaults...)))
    end
    for i in 1:6
        @test (μ_test[i] - μ_list[i]) / μ_list[i] < 0.01
    end
end

@testset "Cc" begin
    Dp_ = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
    Cc_list = [216, 108, 43.6, 22.2, 11.4, 4.95, 2.85, 1.865, 1.326, 1.164, 1.082, 1.032, 1.016, 1.008, 1.003, 1.0016]
    Cc_test = []
    for i in 1:16
        push!(Cc_test, substitute(cc(Dp, T, P, μ), Dict(Dp => Dp_[i] * 1e-6, T => 298, P => 101325, μ => substitute(mu(T), Dict(T => 298, AtmosphericDeposition.defaults...)), AtmosphericDeposition.defaults...)))
    end
    for i in 1:16
        @test (Cc_test[i] - Cc_list[i]) / Cc_list[i] < 0.03
    end
end

@parameters Cc
@testset "Vs" begin
    Dp_ = [0.01, 0.1, 1, 10]
    Vs_list = [0.025, 0.35, 10.8, 1000] ./ 3600 ./ 100 # convert to m/s
    Vs_test = []
    Cc_list = []
    for i in 1:4
        push!(Cc_list, substitute(cc(Dp, T, P, μ), Dict(Dp => Dp_[i] * 1e-6, T => 298, P => 101325, μ => substitute(mu(T), Dict(T => 298, AtmosphericDeposition.defaults...)), AtmosphericDeposition.defaults...)))
        push!(Vs_test, substitute(vs(Dp, ρParticle, Cc, μ), Dict(Dp => Dp_[i] * 1e-6, ρParticle => 1000, Cc => Cc_list[i], μ => 1.836522217711828e-5, AtmosphericDeposition.defaults...)))
    end
    for i in 1:4
        @test (Vs_test[i] - Vs_list[i]) / Vs_list[i] < 1
    end
end

@testset "DryDepGas" begin
    vd_true = 0.03 # m/s
    @test (substitute(DryDepGas(z, z₀, u_star, L, ρA, AtmosphericDeposition.No2Data, G, T, 0, iSeason, iLandUse, false, false, false, false), Dict(z => 50, z₀ => 0.04, u_star => 0.44, L => 0, T => 298, ρA => 1.2, G => 300, iSeason => 1, iLandUse => 10, AtmosphericDeposition.defaults...)) - vd_true) / vd_true < 0.33
end

@testset "DryDepParticle" begin
    Dp_ = [1.e-8, 1.e-7, 1.e-6, 1.e-5]
    vd_true = [0.5, 0.012, 0.02, 0.6] ./ 100 # [m/s]
    vd_list = []
    for i in 1:4
        push!(vd_list, substitute(DryDepParticle(z, z₀, u_star, L, Dp, T, P, ρParticle, ρA, 1, 4), Dict(z => 20, z₀ => 0.02, u_star => 0.44, L => 0, T => 298, P => 101325, ρA => 1.2, ρParticle => 1000, Dp => Dp_[i], AtmosphericDeposition.defaults...)))
    end
    for i in 1:4
        @test vd_list[i] - vd_true[i] < 0.015
    end
end




