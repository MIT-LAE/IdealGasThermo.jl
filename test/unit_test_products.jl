using ForwardDiff

@testset "Vitiator products" begin

    @testset "FAR = 0 reproduces the oxidizer" begin
        sys = Vitiator("CH4")
        @test fieldnames(Vitiator) == (:name, :ηburn, :MWfuel, :sumΔX, :ΔX)
        @test_throws MethodError products(sys, 0.0)
        air = FrozenDryAir
        gas0 = products_in_air(sys, 0.0)
        @test gas0 isa FrozenGas
        for T in [300.0, 1600.0]
            @test IdealGasThermo.cp(gas0, T) ≈ IdealGasThermo.cp(air, T) rtol = 1e-10
            @test IdealGasThermo.h(gas0, T) ≈ IdealGasThermo.h(air, T) rtol = 1e-10
            @test IdealGasThermo.s0(gas0, T) ≈ IdealGasThermo.s0(air, T) rtol = 1e-10
        end

        # The reaction carries no oxidizer: an explicit non-air oxidizer is
        # reproduced at FAR = 0 as well.
        n2 = FrozenGas(species_in_spdict("N2"))
        @test products(sys, n2, 0.0).X == n2.X
    end

    @testset "zero allocations after warmup" begin
        sys = Vitiator("CH4")
        @test (@ballocated products_in_air($sys, 0.03) samples = 1 evals = 1) == 0
        # composed with a property read, still allocation-free
        h_at(sys, FAR, T) = IdealGasThermo.h(products_in_air(sys, FAR), T)
        @test (@ballocated $h_at($sys, 0.03, 1600.0) samples = 1 evals = 1) == 0
    end

    @testset "FAR-differentiability" begin
        sys = Vitiator("CH4")
        D = ForwardDiff.derivative
        h_at(far) = IdealGasThermo.h(products_in_air(sys, far), 1600.0)
        far = 0.03
        dh_ad = D(h_at, far)
        δ = 1e-6
        dh_fd = (h_at(far + δ) - h_at(far - δ)) / (2δ)
        @test dh_ad ≈ dh_fd rtol = 1e-6
    end

    @testset "oxidizer-composition differentiability" begin
        # A general re-burn may receive an oxidizer composition derived from an
        # upstream solve. Its tangent must reach the products' properties through
        # the live oxidizer argument; a same-fuel cumulative-FAR afterburner uses
        # one fixed system instead and does not need this path.
        air = FrozenGas(DryAir)
        iO2 = findfirst(==("O2"), IdealGasThermo.spdict.name)
        iN2 = findfirst(==("N2"), IdealGasThermo.spdict.name)
        o2_0 = air.X[iO2]

        function oxidizer(o2)
            X = map(x -> x + zero(o2), air.X) # promote the composition rail
            X = Base.setindex(X, o2, iO2)
            Base.setindex(X, air.X[iN2] + o2_0 - o2, iN2)
        end

        sys = Vitiator("CH4")
        h_product(o2) = IdealGasThermo.h(products(sys, oxidizer(o2), 0.02), 1600.0)
        d_ad = ForwardDiff.derivative(h_product, o2_0)
        δ = 1e-6
        d_fd_o2 = (h_product(o2_0 + δ) - h_product(o2_0 - δ)) / (2δ)
        @test d_ad ≈ d_fd_o2 rtol = 1e-6

        # The same derivative path occurs when a physical upstream mixer
        # supplies the oxidizer. `mix` returns a FrozenGas whose composition
        # follows `mratio`; products must preserve that tangent without
        # rebuilding the fuel reaction system.
        co2 = FrozenGas(species_in_spdict("CO2"))
        h_mixed_oxidizer(mratio) =
            IdealGasThermo.h(products(sys, mix(air, co2, mratio), 0.02), 1600.0)
        mratio = 0.05
        d_ad = ForwardDiff.derivative(h_mixed_oxidizer, mratio)
        d_fd_mixed =
            (h_mixed_oxidizer(mratio + δ) - h_mixed_oxidizer(mratio - δ)) / (2δ)
        @test d_ad ≈ d_fd_mixed rtol = 1e-6

        h_product_vector(o2) =
            IdealGasThermo.h(products(sys, collect(oxidizer(o2)), 0.02), 1600.0)
        @test ForwardDiff.derivative(h_product_vector, o2_0) ≈ d_fd_o2 rtol = 1e-6

        Xdual = oxidizer(ForwardDiff.Dual{:oxidizer}(o2_0, 1.0))
        @test eltype(products(sys, Xdual, 0.02).X) <: ForwardDiff.Dual
        @test eltype(products(sys, collect(Xdual), 0.02).X) <: ForwardDiff.Dual
        @test (@ballocated products($sys, $Xdual, 0.02) samples = 1 evals = 1) == 0
        @test products(sys, FrozenDryAir, 0.02).X ≈ products_in_air(sys, 0.02).X

        # The concrete re-burn shape: the oxidizer is itself a product gas
        # whose composition moves with an upstream FAR.
        upstream = Vitiator("CH4")
        h_reburn(far) = IdealGasThermo.h(
            products(sys, products_in_air(upstream, far), 0.005), 1600.0)
        d_ad = ForwardDiff.derivative(h_reburn, 0.01)
        d_fd = (h_reburn(0.01 + δ) - h_reburn(0.01 - δ)) / (2δ)
        @test d_ad ≈ d_fd rtol = 1e-6

        # A downstream burner may use a different fuel. Each Vitiator can be
        # constructed once, then consume the live product gas from upstream.
        downstream = Vitiator("H2")
        h_two_fuel_burners(far) = IdealGasThermo.h(
            products(downstream, products_in_air(upstream, far), 0.005), 1600.0)
        d_ad = ForwardDiff.derivative(h_two_fuel_burners, 0.01)
        d_fd = (h_two_fuel_burners(0.01 + δ) - h_two_fuel_burners(0.01 - δ)) / (2δ)
        @test d_ad ≈ d_fd rtol = 1e-6
    end

    @testset "incomplete combustion (ηburn ≠ 1)" begin
        # ηburn is a real composition knob, not a passthrough: at the SAME FAR,
        # burning only a fraction of the fuel leaves a different product mixture
        # (more unburnt fuel / less CO2 + H2O) than complete combustion, so the
        # two gases must have measurably different properties. (Self-contained;
        # no comparison to the legacy vitiated_species path.)
        FAR = 0.03
        full = products_in_air(Vitiator("CH4"; ηburn = 1.0), FAR)
        partial = products_in_air(Vitiator("CH4"; ηburn = 0.9), FAR)
        # the compositions differ ⟹ cp differs by well over numerical noise
        @test !isapprox(IdealGasThermo.cp(partial, 1600.0),
                        IdealGasThermo.cp(full, 1600.0); rtol = 1e-3)
        # and ηburn has no effect when there is no fuel to burn: at FAR = 0 both
        # collapse to the pure oxidizer regardless of ηburn
        air = FrozenDryAir
        for ηburn in (0.9, 1.0)
            gas0 = products_in_air(Vitiator("CH4"; ηburn = ηburn), 0.0)
            @test IdealGasThermo.cp(gas0, 1600.0) ≈ IdealGasThermo.cp(air, 1600.0) rtol = 1e-10
        end
    end

end
