using DepositionMTK
using Test, SafeTestsets

@testset "DepositionMTK.jl" begin
    @safetestset "drydep_test.jl" begin
        include("drydep_test.jl")
    end
    
    @safetestset "wesely1989_test.jl" begin
        include("wesely1989_test.jl")
    end
end
