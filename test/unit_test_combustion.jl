@testset "combustion" begin
    CH4 = species_in_spdict("CH4")
    # Regression pin computed by this package; consistent with the literature
    # LHV of methane, ~50.0 MJ/kg.
    @test IdealGasThermo.LHV(CH4) == 5.002736488044851e7
    @test IdealGasThermo.LHV("CH4") == 5.002736488044851e7

    # Regression pins computed by this package. These are COMPLETE-combustion,
    # frozen-composition (no dissociation) flame temperatures, so they sit above
    # equilibrium literature values (CH4–air ≈ 2230 K, CH4–O2 ≈ 3080 K with
    # dissociation); the gap is large for pure O2 where dissociation dominates.
    @test IdealGasThermo.AFT(CH4) == 2376.6102617357988
    O2 = species_in_spdict("O2")
    @test IdealGasThermo.AFT(CH4, O2) ≈ 5280.53933225877 rtol = 1e-8
end
@testset "burnt gas funcs." begin
    CH4 = species_in_spdict("CH4")

    FAR = 0.05
    gas1 = IdealGasThermo.vitiated_species(CH4, DryAir, FAR)

    burntgas = IdealGasThermo.fixed_fuel_vitiated_species(CH4, DryAir)
    gas2 = burntgas(FAR)

    @test gas1.MW ≈ gas2.MW
    for (key, val) in gas1.composition
        @test gas1.composition[key] ≈ gas2.composition[key]
    end

    # Regression pins computed by this package (H2 in dry air, complete
    # combustion): outlet T at FAR = 0.01, and the FAR that reaches 1000 K.
    gas = Gas()
    gasburnt = IdealGasThermo.fuel_combustion(gas, "H2", 298.15, 0.01)
    @test gasburnt.T == 1293.4126150619875

    gas = Gas()
    FAR, _ = IdealGasThermo.gas_burn(gas, "H2", 298.15, 1000.0)
    @test FAR == 0.006636436988741555

end
