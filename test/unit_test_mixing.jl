using ForwardDiff

@testset "mix" begin

    @testset "mratio = 0 reproduces stream 1; mratio → ∞ approaches stream 2" begin
        vit = IdealGasThermo.vitiated_species("CH4", "Air", 0.03)
        gas1 = FrozenGas(DryAir)
        gas2 = FrozenGas(vit)

        gas0 = mix(gas1, gas2, 0.0)
        @test gas0 isa FrozenGas
        for T in [300.0, 1600.0]
            @test IdealGasThermo.cp(gas0, T) ≈ IdealGasThermo.cp(gas1, T) rtol = 1e-10
            @test IdealGasThermo.h(gas0, T) ≈ IdealGasThermo.h(gas1, T) rtol = 1e-10
            @test IdealGasThermo.s0(gas0, T) ≈ IdealGasThermo.s0(gas1, T) rtol = 1e-10
        end

        gasinf = mix(gas1, gas2, 1e9)
        for T in [300.0, 1600.0]
            @test IdealGasThermo.cp(gasinf, T) ≈ IdealGasThermo.cp(gas2, T) rtol = 1e-6
            @test IdealGasThermo.h(gasinf, T) ≈ IdealGasThermo.h(gas2, T) rtol = 1e-6
            @test IdealGasThermo.s0(gasinf, T) ≈ IdealGasThermo.s0(gas2, T) rtol = 1e-6
        end
    end

    @testset "energy balance: adiabatic mix conserves total enthalpy" begin
        P = 2.5e5
        core = GasState(products(Vitiator("CH4", DryAir), 0.03), 1500.0, P)
        bypass = GasState(FrozenGas(DryAir), 320.0, P)
        BPR = 5.0
        m = mix(core, bypass, BPR)
        @test m isa GasState
        @test m.P == P
        @test bypass.T < m.T < core.T
        Hin = IdealGasThermo.h(core.gas, core.T) + BPR * IdealGasThermo.h(bypass.gas, bypass.T)
        Hout = (1 + BPR) * IdealGasThermo.h(m.gas, m.T)
        @test Hout ≈ Hin rtol = 1e-12
        # unequal pressures need a momentum closure — rejected here
        @test_throws ArgumentError mix(core, GasState(FrozenGas(DryAir), 320.0, P / 2), BPR)
    end

    @testset "zero allocations after warmup" begin
        gas1 = FrozenGas(DryAir)
        gas2 = FrozenGas(IdealGasThermo.vitiated_species("CH4", "Air", 0.03))
        @test (@ballocated mix($gas1, $gas2, 0.25) samples = 1 evals = 1) == 0
        # composed with a property read, still allocation-free
        h_at(a, b, mratio, T) = IdealGasThermo.h(mix(a, b, mratio), T)
        @test (@ballocated $h_at($gas1, $gas2, 0.25, 1600.0) samples = 1 evals = 1) == 0
    end

    @testset "mratio-differentiability" begin
        gas1 = FrozenGas(DryAir)
        gas2 = FrozenGas(IdealGasThermo.vitiated_species("CH4", "Air", 0.03))
        D = ForwardDiff.derivative
        h_at(mr) = IdealGasThermo.h(mix(gas1, gas2, mr), 1600.0)
        mr = 0.25
        dh_ad = D(h_at, mr)
        δ = 1e-6
        dh_fd = (h_at(mr + δ) - h_at(mr - δ)) / (2δ)
        @test dh_ad ≈ dh_fd rtol = 1e-6
    end

end
