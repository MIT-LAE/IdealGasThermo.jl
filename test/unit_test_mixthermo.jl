# Legacy smoke (Gas{N} composition + thermo): kept until the mutable layer is
# removed at v2.0 (ADR-0002). The CEA "Air" temperature-range table that used
# to live here now anchors the PURE core (FrozenGas) in
# test/unit_test_cea_reference.jl — it should validate the pure core, not the
# legacy Gas/thermo_table path.

@testset "mix. comp." begin
    gas = Gas()
    Yair = Dict(
        "O2" => 0.231416,
        "Ar" => 0.012916,
        "Air" => 0.0,
        "H2O" => 0.0,
        "CO2" => 0.000484688,
        "N2" => 0.755184,
    )

    gas.X = IdealGasThermo.Xair

    for (key, val) in gas.Ydict
        if val != 0
            @test Yair[key] ≈ val atol = 1e-6
        end
    end
    Air = species_in_spdict("Air")
    @test gas.MW ≈ Air.MW

end

@testset "mix. thermo" begin
    gas = Gas()
    gas.X = IdealGasThermo.Xair

    Air = species_in_spdict("Air")

    # Low temp:
    gas.T = T = IdealGasThermo.Tstd
    gas.P = P = IdealGasThermo.Pstd

    @test gas.cp ≈ IdealGasThermo.Cp(T, Air) rtol = 1e-7
    @test gas.h ≈ IdealGasThermo.h(T, Air) rtol = 1e-7
    @test gas.s ≈ IdealGasThermo.s(T, P, Air) rtol = 1e-7

    # High temp:
    gas.T = T = 20 * T
    gas.P = P = 20 * P
    @test gas.cp ≈ IdealGasThermo.Cp(T, Air) rtol = 1e-7
    @test gas.h ≈ IdealGasThermo.h(T, Air) rtol = 1e-7
    @test gas.s ≈ IdealGasThermo.s(T, P, Air) rtol = 1e-7

end
