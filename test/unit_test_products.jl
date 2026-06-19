using ForwardDiff

@testset "Vitiator products" begin

    @testset "FAR = 0 reproduces the oxidizer" begin
        sys = Vitiator("CH4", DryAir)
        air = FrozenGas(DryAir)
        gas0 = products(sys, 0.0)
        @test gas0 isa FrozenGas
        for T in [300.0, 1600.0]
            @test IdealGasThermo.cp(gas0, T) ≈ IdealGasThermo.cp(air, T) rtol = 1e-10
            @test IdealGasThermo.h(gas0, T) ≈ IdealGasThermo.h(air, T) rtol = 1e-10
            @test IdealGasThermo.s0(gas0, T) ≈ IdealGasThermo.s0(air, T) rtol = 1e-10
        end
    end

    @testset "zero allocations after warmup" begin
        sys = Vitiator("CH4", "Air")
        @test (@ballocated products($sys, 0.03) samples = 1 evals = 1) == 0
        # composed with a property read, still allocation-free
        h_at(sys, FAR, T) = IdealGasThermo.h(products(sys, FAR), T)
        @test (@ballocated $h_at($sys, 0.03, 1600.0) samples = 1 evals = 1) == 0
    end

    @testset "FAR-differentiability" begin
        sys = Vitiator("CH4", "Air")
        D = ForwardDiff.derivative
        h_at(far) = IdealGasThermo.h(products(sys, far), 1600.0)
        far = 0.03
        dh_ad = D(h_at, far)
        δ = 1e-6
        dh_fd = (h_at(far + δ) - h_at(far - δ)) / (2δ)
        @test dh_ad ≈ dh_fd rtol = 1e-6
    end

    @testset "incomplete combustion (ηburn ≠ 1)" begin
        # ηburn is a real composition knob, not a passthrough: at the SAME FAR,
        # burning only a fraction of the fuel leaves a different product mixture
        # (more unburnt fuel / less CO2 + H2O) than complete combustion, so the
        # two gases must have measurably different properties. (Self-contained;
        # no comparison to the legacy vitiated_species path.)
        FAR = 0.03
        full = products(Vitiator("CH4", "Air"; ηburn = 1.0), FAR)
        partial = products(Vitiator("CH4", "Air"; ηburn = 0.9), FAR)
        # the compositions differ ⟹ cp differs by well over numerical noise
        @test !isapprox(IdealGasThermo.cp(partial, 1600.0),
                        IdealGasThermo.cp(full, 1600.0); rtol = 1e-3)
        # and ηburn has no effect when there is no fuel to burn: at FAR = 0 both
        # collapse to the pure oxidizer regardless of ηburn
        air = FrozenGas(DryAir)
        for ηburn in (0.9, 1.0)
            gas0 = products(Vitiator("CH4", DryAir; ηburn = ηburn), 0.0)
            @test IdealGasThermo.cp(gas0, 1600.0) ≈ IdealGasThermo.cp(air, 1600.0) rtol = 1e-10
        end
    end

end
