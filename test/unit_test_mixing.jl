# NOTE: substantial parts of this file are MIGRATION TESTS — they pin the
# pure core to the legacy implementation (vitiated_species / Gas1D /
# set_Δh! era) to prove the refactor preserved behavior. They retire with
# the legacy layer in v2.0. The physics itself is guarded independently in
# unit_test_properties.jl.
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

    # Cross-check against the full legacy gas_mixing on Gas{N} objects.
    # Compared: the merged COMPOSITION only — gas_mixing's resulting gas.X is
    # fed through generate_composite_species/FrozenGas and its cp(T) must
    # match mixed(sys, mratio). NOT compared: the legacy outlet temperature
    # (gas_mixing performs its own constant-enthalpy set_hP! solve, i.e. an
    # energy balance; mixed() returns a FrozenGas whose temperature is the
    # caller's argument) and the legacy cached gas.cp at that temperature.
    @testset "cross-check vs legacy gas_mixing composition" begin
        vit = IdealGasThermo.vitiated_species("CH4", "Air", 0.03)
        sys = Mixer(DryAir, vit)

        gas2 = Gas()
        gas2.X = vit.composition
        for mratio in (0.25, 1.0, 4.0)
            gas1 = Gas() # air at standard T, P ("Air" maps to Xair in gas_mixing)
            gas_prod = IdealGasThermo.gas_mixing(gas1, gas2, mratio)
            legacy = FrozenGas(generate_composite_species(gas_prod.X))
            gas = mixed(sys, mratio)
            for T in [300.0, 1600.0]
                @test IdealGasThermo.cp(gas, T) ≈ IdealGasThermo.cp(legacy, T) rtol = 1e-8
            end
        end
    end

    @testset "agreement with legacy law of mixtures (DryAir + vitiated CH4)" begin
        vit = IdealGasThermo.vitiated_species("CH4", "Air", 0.03)
        sys = Mixer(DryAir, vit)

        # Legacy merged composition via the existing utilities, exactly as
        # gas_mixing does it: X -> Y per stream, mass-fraction law of
        # mixtures, back to X, then a composite species.
        X1 = IdealGasThermo.Xidict2Array(DryAir.composition)
        X2 = IdealGasThermo.Xidict2Array(vit.composition)
        Y1 = X2Y(X1)
        Y2 = X2Y(X2)
        for mratio in (0.25, 1.0, 4.0)
            Yp = (Y1 + mratio * Y2) / (1 + mratio)
            X_merged = Y2X(Yp)
            legacy = FrozenGas(generate_composite_species(X_merged))
            gas = mixed(sys, mratio)
            @test gas.MW ≈ legacy.MW rtol = 1e-10
            @test IdealGasThermo.R(gas) ≈ IdealGasThermo.R(legacy) rtol = 1e-10
            for T in [300.0, 1600.0]
                @test IdealGasThermo.cp(gas, T) ≈ IdealGasThermo.cp(legacy, T) rtol = 1e-10
                @test IdealGasThermo.h(gas, T) ≈ IdealGasThermo.h(legacy, T) rtol = 1e-10
                @test IdealGasThermo.s0(gas, T) ≈ IdealGasThermo.s0(legacy, T) rtol = 1e-10
            end
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
        for mr in [0.25, 1.0, 4.0]
            dh_ad = D(h_at, mr)
            δ = 1e-6
            dh_fd = (h_at(mr + δ) - h_at(mr - δ)) / (2δ)
            @test dh_ad ≈ dh_fd rtol = 1e-6
        end
    end

end
