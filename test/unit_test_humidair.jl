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

    @testset "RH ↔ SH round-trip consistency" begin
        # Specifying humidity by relative humidity or by the equivalent specific
        # humidity must build the same gas. RH is mapped to SH through the
        # saturation curve at the standard-day (T, P) — a distinct conversion
        # claim, not a passthrough — so converting RH → SH and back to a gas
        # must reproduce the RH-built gas.
        RH = 0.4
        SH = IdealGasThermo.specific_humidity(RH, IdealGasThermo.Tstd, IdealGasThermo.Pstd)
        gasRH = humid_air(RH = RH)
        gasSH = humid_air(SH = SH)
        @test gasSH.MW ≈ gasRH.MW rtol = 1e-12
        for T in [300.0, 1600.0]
            @test IdealGasThermo.cp(gasSH, T) ≈ IdealGasThermo.cp(gasRH, T) rtol = 1e-12
            @test IdealGasThermo.h(gasSH, T) ≈ IdealGasThermo.h(gasRH, T) rtol = 1e-12
        end
    end

    # Water (MW ≈ 18) is lighter than dry air (MW ≈ 29), so adding vapor
    # lowers the mixture MW: R = 1000·Runiv/MW must increase strictly
    # monotonically with specific humidity.
    @testset "R increases monotonically with humidity" begin
        air = FrozenGas(DryAir)
        ωs = range(0.0, 0.05, length = 3)
        Rs = [IdealGasThermo.R(humid_air(SH = ω)) for ω in ωs]
        @test Rs[1] ≈ IdealGasThermo.R(air) rtol = 1e-12
        @test all(diff(Rs) .> 0)
        @test all(R -> R > IdealGasThermo.R(air), Rs[2:end])
    end

end
