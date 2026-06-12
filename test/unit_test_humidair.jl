# NOTE: substantial parts of this file are MIGRATION TESTS — they pin the
# pure core to the legacy implementation (vitiated_species / Gas1D /
# set_Δh! era) to prove the refactor preserved behavior. They retire with
# the legacy layer in v2.0. The physics itself is guarded independently in
# unit_test_properties.jl.
@testset "humid_air" begin

    @testset "zero humidity reproduces dry air" begin
        air = FrozenGas(DryAir)
        for gas0 in (humid_air(SH = 0.0), humid_air(RH = 0.0))
            @test gas0 isa FrozenGas
            for T in [300.0, 1600.0]
                @test IdealGasThermo.cp(gas0, T) ≈ IdealGasThermo.cp(air, T) rtol = 1e-10
                @test IdealGasThermo.h(gas0, T) ≈ IdealGasThermo.h(air, T) rtol = 1e-10
            end
        end
    end

    @testset "agreement with legacy generate_humid_air" begin
        for RH in (0.1, 0.5, 1.0)
            legacy = FrozenGas(IdealGasThermo.generate_humid_air(RH))
            gas = humid_air(RH = RH) # same standard-day T, P defaults
            @test gas.MW ≈ legacy.MW rtol = 1e-12
            for T in [300.0, 1600.0]
                @test IdealGasThermo.cp(gas, T) ≈ IdealGasThermo.cp(legacy, T) rtol = 1e-10
                @test IdealGasThermo.h(gas, T) ≈ IdealGasThermo.h(legacy, T) rtol = 1e-10
                @test IdealGasThermo.s0(gas, T) ≈ IdealGasThermo.s0(legacy, T) rtol = 1e-10
            end
        end
        # round trip through the legacy specific-humidity definitions:
        # SH(composite of humid_air(RH)) == SH(RH, T, P)
        RH = 0.4
        SH = IdealGasThermo.specific_humidity(RH, IdealGasThermo.Tstd, IdealGasThermo.Pstd)
        gasRH = humid_air(RH = RH)
        gasSH = humid_air(SH = SH)
        for T in [300.0, 1600.0]
            @test IdealGasThermo.cp(gasSH, T) ≈ IdealGasThermo.cp(gasRH, T) rtol = 1e-12
            @test IdealGasThermo.h(gasSH, T) ≈ IdealGasThermo.h(gasRH, T) rtol = 1e-12
        end
    end

    @testset "agreement with the XwetAir constant" begin
        # XwetAir is Xair plus 0.018722 moles of H2O per mole of dry air
        # (entries are unnormalized mole ratios), i.e. specific humidity
        # ω = ε·0.018722.
        Xwet = IdealGasThermo.Xidict2Array(IdealGasThermo.XwetAir)
        Xwet = Xwet ./ sum(Xwet)
        reference = FrozenGas(generate_composite_species(Xwet))
        gas = humid_air(SH = IdealGasThermo.ε * 0.018722)
        @test gas.MW ≈ reference.MW rtol = 1e-12
        for T in [300.0, 1600.0]
            @test IdealGasThermo.cp(gas, T) ≈ IdealGasThermo.cp(reference, T) rtol = 1e-10
            @test IdealGasThermo.h(gas, T) ≈ IdealGasThermo.h(reference, T) rtol = 1e-10
            @test IdealGasThermo.s0(gas, T) ≈ IdealGasThermo.s0(reference, T) rtol = 1e-10
        end
    end

    # Water (MW ≈ 18) is lighter than dry air (MW ≈ 29), so adding vapor
    # lowers the mixture MW: R = 1000·Runiv/MW must increase strictly
    # monotonically with specific humidity.
    @testset "R increases monotonically with humidity" begin
        air = FrozenGas(DryAir)
        ωs = range(0.0, 0.05, length = 11)
        Rs = [IdealGasThermo.R(humid_air(SH = ω)) for ω in ωs]
        @test Rs[1] ≈ IdealGasThermo.R(air) rtol = 1e-12
        @test all(diff(Rs) .> 0)
        @test all(R -> R > IdealGasThermo.R(air), Rs[2:end])
    end

end
