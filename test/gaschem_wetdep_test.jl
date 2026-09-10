# Tests for the GasChem GEOSChemGasPhase <-> wet-deposition coupling (GasChemExt):
# verify that O3 and NO2 are NOT wet-scavenged (GEOS-Chem does not wet-deposit
# these sparingly-soluble gases), while a soluble species (HNO3) still is.

@testitem "GEOSChem wet-dep excludes insoluble O3 and NO2" begin
    using AtmosphericDeposition, GasChem, EarthSciMLBase, ModelingToolkit

    sys = convert(System, couple(GEOSChemGasPhase(), WetDeposition()))

    # RHS of a species' tendency equation, identified by its LHS (D(species)).
    # Checking the RHS of the specific tendency avoids false positives from the
    # species appearing as a reactant in another (wet-deposited) species' equation.
    function tendency_rhs(sp)
        for eq in equations(sys)
            occursin("₊" * sp * "(t)", string(eq.lhs)) && return string(eq.rhs)
        end
        return nothing
    end

    # The wet-dep operator adds a `k_othergas` first-order sink. O3 and NO2 must not.
    @test tendency_rhs("O3") !== nothing
    @test !occursin("k_othergas", tendency_rhs("O3"))
    @test tendency_rhs("NO2") !== nothing
    @test !occursin("k_othergas", tendency_rhs("NO2"))
    # Sanity: a soluble species (HNO3) still IS wet-deposited via k_othergas.
    @test occursin("k_othergas", tendency_rhs("HNO3"))
end
