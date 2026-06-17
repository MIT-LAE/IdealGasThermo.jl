using ForwardDiff

@testset "Mixer mixed" begin

    @testset "mratio = 0 reproduces stream 1; mratio → ∞ approaches stream 2" begin
        vit = IdealGasThermo.vitiated_species("CH4", "Air", 0.03)
        sys = Mixer(DryAir, vit)
        gas1 = FrozenGas(DryAir)
        gas2 = FrozenGas(vit)

        gas0 = mixed(sys, 0.0)
        @test gas0 isa FrozenGas
        for T in [300.0, 1600.0]
            @test IdealGasThermo.cp(gas0, T) ≈ IdealGasThermo.cp(gas1, T) rtol = 1e-10
            @test IdealGasThermo.h(gas0, T) ≈ IdealGasThermo.h(gas1, T) rtol = 1e-10
            @test IdealGasThermo.s0(gas0, T) ≈ IdealGasThermo.s0(gas1, T) rtol = 1e-10
        end

        gasinf = mixed(sys, 1e9)
        for T in [300.0, 1600.0]
            @test IdealGasThermo.cp(gasinf, T) ≈ IdealGasThermo.cp(gas2, T) rtol = 1e-6
            @test IdealGasThermo.h(gasinf, T) ≈ IdealGasThermo.h(gas2, T) rtol = 1e-6
            @test IdealGasThermo.s0(gasinf, T) ≈ IdealGasThermo.s0(gas2, T) rtol = 1e-6
        end
    end

    @testset "zero allocations after warmup" begin
        vit = IdealGasThermo.vitiated_species("CH4", "Air", 0.03)
        sys = Mixer(DryAir, vit)
        measured(f::F, args...) where {F} = (f(args...); @allocated f(args...))
        @test measured(mixed, sys, 0.25) == 0
        # composed with a property read, still allocation-free
        h_at(sys, mratio, T) = IdealGasThermo.h(mixed(sys, mratio), T)
        @test measured(h_at, sys, 0.25, 1600.0) == 0
    end

    @testset "mratio-differentiability" begin
        vit = IdealGasThermo.vitiated_species("CH4", "Air", 0.03)
        sys = Mixer(DryAir, vit)
        D = ForwardDiff.derivative
        h_at(mr) = IdealGasThermo.h(mixed(sys, mr), 1600.0)
        mr = 0.25
        dh_ad = D(h_at, mr)
        δ = 1e-6
        dh_fd = (h_at(mr + δ) - h_at(mr - δ)) / (2δ)
        @test dh_ad ≈ dh_fd rtol = 1e-6
    end

end
