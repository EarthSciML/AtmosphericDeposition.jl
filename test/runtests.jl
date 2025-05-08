using AtmosphericDeposition
using Test, SafeTestsets

@testset "AtmosphericDeposition.jl" begin
    @safetestset "drydep_test.jl" begin
        include("drydep_test.jl")
    end

    @safetestset "wesely1989_test.jl" begin
        include("wesely1989_test.jl")
    end

    @safetestset "wetdep_test.jl" begin
        include("wetdep_test.jl")
    end

    @safetestset "connector_test.jl" begin
        include("connector_test.jl")
    end
end
