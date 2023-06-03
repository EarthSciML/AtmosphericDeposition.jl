using DepositionMTK
using Test,Unitful

@testset "unit" begin
    @test unit(WetDeposition(0.5,0.5,1.204u"kg*m^-3",1u"m")[1]) == u"s^-1"
    @test unit(WetDeposition(0.5,0.5,1.204u"kg*m^-3",1u"m")[2]) == u"s^-1"
    @test unit(WetDeposition(0.5,0.5,1.204u"kg*m^-3",1u"m")[3]) == u"s^-1"
end
