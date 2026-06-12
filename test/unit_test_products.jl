using ForwardDiff

@testset "Combustor products" begin

    @testset "FAR = 0 reproduces the oxidizer" begin
        sys = Combustor("CH4", DryAir)
        air = FrozenGas(DryAir)
        gas0 = products(sys, 0.0)
        @test gas0 isa FrozenGas
        for T in [300.0, 1600.0]
            @test IdealGasThermo.cp(gas0, T) ≈ IdealGasThermo.cp(air, T) rtol = 1e-10
            @test IdealGasThermo.h(gas0, T) ≈ IdealGasThermo.h(air, T) rtol = 1e-10
            @test IdealGasThermo.s0(gas0, T) ≈ IdealGasThermo.s0(air, T) rtol = 1e-10
        end
    end

    @testset "agreement with legacy vitiated_species (CH4 + Air)" begin
        sys = Combustor("CH4", "Air")
        legacy = FrozenGas(IdealGasThermo.vitiated_species("CH4", "Air", 0.03))
        gas = products(sys, 0.03)
        @test gas.MW ≈ legacy.MW rtol = 1e-12
        @test IdealGasThermo.R(gas) ≈ IdealGasThermo.R(legacy) rtol = 1e-12
        for T in [300.0, 999.9, 1000.0, 1600.0, 2200.0]
            @test IdealGasThermo.cp(gas, T) ≈ IdealGasThermo.cp(legacy, T) rtol = 1e-10
            @test IdealGasThermo.h(gas, T) ≈ IdealGasThermo.h(legacy, T) rtol = 1e-10
            @test IdealGasThermo.s0(gas, T) ≈ IdealGasThermo.s0(legacy, T) rtol = 1e-10
        end
    end

    @testset "agreement with legacy vitiated_species (Jet-A + Air)" begin
        FAR = 0.029681
        sys = Combustor("Jet-A(g)", "Air")
        legacy = FrozenGas(IdealGasThermo.vitiated_species("Jet-A(g)", "Air", FAR))
        gas = products(sys, FAR)
        @test gas.MW ≈ legacy.MW rtol = 1e-12
        @test IdealGasThermo.R(gas) ≈ IdealGasThermo.R(legacy) rtol = 1e-12
        for T in [300.0, 999.9, 1000.0, 1600.0, 2200.0]
            @test IdealGasThermo.cp(gas, T) ≈ IdealGasThermo.cp(legacy, T) rtol = 1e-10
            @test IdealGasThermo.h(gas, T) ≈ IdealGasThermo.h(legacy, T) rtol = 1e-10
            @test IdealGasThermo.s0(gas, T) ≈ IdealGasThermo.s0(legacy, T) rtol = 1e-10
        end
        # Requirements acceptance values: mean cp [J/kg/K] over two intervals
        h1600 = IdealGasThermo.h(gas, 1600.0)
        @test (h1600 - IdealGasThermo.h(gas, 298.15)) / 1301.85 ≈ 1172.7 rtol = 0.015
        @test (h1600 - IdealGasThermo.h(gas, 1373.78)) / 226.22 ≈ 1276.9 rtol = 0.015
    end

    @testset "zero allocations after warmup" begin
        sys = Combustor("CH4", "Air")
        measured(f::F, args...) where {F} = (f(args...); @allocated f(args...))
        @test measured(products, sys, 0.03) == 0
        # composed with a property read, still allocation-free
        h_at(sys, FAR, T) = IdealGasThermo.h(products(sys, FAR), T)
        @test measured(h_at, sys, 0.03, 1600.0) == 0
    end

    @testset "FAR-differentiability" begin
        sys = Combustor("CH4", "Air")
        D = ForwardDiff.derivative
        h_at(far) = IdealGasThermo.h(products(sys, far), 1600.0)
        for far in [0.03, 0.01]
            dh_ad = D(h_at, far)
            δ = 1e-6
            dh_fd = (h_at(far + δ) - h_at(far - δ)) / (2δ)
            @test dh_ad ≈ dh_fd rtol = 1e-6
        end
    end

    @testset "incomplete combustion (ηburn ≠ 1)" begin
        FAR, ηburn = 0.03, 0.9
        sys = Combustor("CH4", "Air"; ηburn = ηburn)
        legacy = FrozenGas(IdealGasThermo.vitiated_species("CH4", "Air", FAR; ηburn = ηburn))
        gas = products(sys, FAR)
        @test gas.MW ≈ legacy.MW rtol = 1e-10
        for T in [300.0, 1600.0]
            @test IdealGasThermo.cp(gas, T) ≈ IdealGasThermo.cp(legacy, T) rtol = 1e-10
            @test IdealGasThermo.h(gas, T) ≈ IdealGasThermo.h(legacy, T) rtol = 1e-10
            @test IdealGasThermo.s0(gas, T) ≈ IdealGasThermo.s0(legacy, T) rtol = 1e-10
        end
    end

end
