# Physics property tests: invariants that must hold over the WHOLE space of
# gases this package can construct — not at hand-picked points, and not by
# agreement with the legacy implementation (see the migration tests for
# that). Randomized with a fixed seed and deliberately small N: each @test
# here is a distinct physical claim, not a coverage statistic.
using Random
using ForwardDiff

@testset "physics properties (randomized, seeded)" begin

    rng = Xoshiro(20260612)
    uni(a, b) = a + (b - a) * rand(rng) #uniform random between a and b

    air = FrozenGas(DryAir)
    sysCH4 = Vitiator("CH4", DryAir)
    sysJet = Vitiator("Jet-A(g)", DryAir)
    vit = IdealGasThermo.vitiated_species("CH4", "Air", 0.04)

    # a small pool of randomly-parameterized gases spanning what the
    # constructors can produce
    gas_pool() = [
        air,
        products(sysCH4, uni(0.005, 0.05)),
        products(sysJet, uni(0.005, 0.045)),
        mix(FrozenGas(DryAir), FrozenGas(vit), uni(0.1, 10.0)),
        humid_air(SH = uni(0.001, 0.05)),
    ]

    @testset "positivity, monotonicity, gamma bounds (random gases)" begin
        Tgrid = range(220.0, 2350.0; length = 40)
        for gas in gas_pool()
            hs = [IdealGasThermo.h(gas, T) for T in Tgrid]
            ss = [IdealGasThermo.s0(gas, T) for T in Tgrid]
            @test all(IdealGasThermo.cp(gas, T) > 0 for T in Tgrid)
            @test all(diff(hs) .> 0)  # dh/dT = cp > 0 ⟹ strictly increasing
            @test all(diff(ss) .> 0)  # ds0/dT = cp/T > 0
            @test all(1.0 < IdealGasThermo.gamma(gas, T) < 1.7 for T in Tgrid)
        end
    end

    @testset "second law: irreversible round trip generates entropy" begin
        for _ = 1:5
            st = GasState(air, uni(250.0, 800.0), 1.0e5)
            PR, ηp = uni(2.0, 30.0), uni(0.85, 0.95)
            rt = expand(compress(st, PR; ηp = ηp), PR; ηp = ηp)
            @test rt.P ≈ st.P rtol = 1e-12      # pressure round trip exact
            @test rt.T > st.T                   # lost work shows up as heat
            @test entropy(rt) > entropy(st)     # ΔS > 0, strictly
        end
        # the reversible limit recovers the identity
        st = GasState(air, 400.0, 1.0e5)
        rt = expand(compress(st, 12.0), 12.0)
        @test rt.T ≈ st.T rtol = 1e-9
    end

    @testset "process round trips (reversible limits)" begin
        for _ = 1:3
            st = GasState(air, uni(300.0, 900.0), uni(0.5e5, 5.0e5))
            q = uni(1.0e4, 5.0e5)
            back = add_heat(add_heat(st, q), -q)
            @test back.T ≈ st.T rtol = 1e-9
            @test back.P == st.P                # constant-P by construction
            w = uni(1.0e4, 3.0e5)
            rt = extract_work(add_work(st, w), w) # ηp = 1: thermodynamically reversible
            @test rt.T ≈ st.T rtol = 1e-9
            @test rt.P ≈ st.P rtol = 1e-9
        end
    end

    @testset "entropy pressure-scaling: s(T, kP) = s(T, P) - R ln k" begin
        for _ = 1:3
            st = GasState(air, uni(250.0, 1800.0), uni(0.2e5, 2.0e5))
            k = uni(0.1, 40.0)
            @test entropy(GasState(air, st.T, k * st.P)) ≈
                  entropy(st) - IdealGasThermo.R(air) * log(k) rtol = 1e-12
        end
    end

    @testset "combustion: mass conservation and heat release vs LHV" begin
        for (fuel, sys, FARmax) in (("CH4", sysCH4, 0.05), ("Jet-A(g)", sysJet, 0.045))
            fsp = species_in_spdict(fuel)
            ΣΔX = sum(IdealGasThermo.reaction_change_molar_fraction(fuel))
            for _ = 1:3
                f = uni(0.005, FARmax)
                gp = products(sys, f)
                # mass conservation through the mole-fraction algebra:
                # (1 + molFAR·ΣΔX)·MW_products = MW_ox + molFAR·MW_fuel
                molFAR = f * DryAir.MW / fsp.MW
                @test (1 + molFAR * ΣΔX) * gp.MW ≈
                      DryAir.MW + molFAR * fsp.MW rtol = 1e-12
            end
            # heat release at 298.15 K equals f·LHV — ties the products
            # composition algebra, the formation-inclusive datum (b₁
            # coefficients), and the independent Hf-header-based LHV()
            # computation together
            f = 0.03
            fgas = FrozenGas(fsp)
            release =
                -(
                    (1 + f) * IdealGasThermo.h(products(sys, f), 298.15) -
                    IdealGasThermo.h(air, 298.15) - f * IdealGasThermo.h(fgas, 298.15)
                )
            @test release ≈ f * IdealGasThermo.LHV(fsp) rtol = 1e-5
        end
    end

    @testset "mixing algebra: symmetry and self-identity" begin
        gA = FrozenGas(DryAir)
        gB = FrozenGas(vit)
        for _ = 1:4
            m = uni(0.1, 10.0)
            g1 = mix(gA, gB, m)      # m kg of vit per kg of air
            g2 = mix(gB, gA, 1 / m)  # the same mixture, described from the other side
            @test g1.MW ≈ g2.MW rtol = 1e-12
            for T in (300.0, 1600.0)
                @test IdealGasThermo.cp(g1, T) ≈ IdealGasThermo.cp(g2, T) rtol = 1e-12
                @test IdealGasThermo.h(g1, T) ≈ IdealGasThermo.h(g2, T) rtol = 1e-12
                @test IdealGasThermo.s0(g1, T) ≈ IdealGasThermo.s0(g2, T) rtol = 1e-12
            end
        end
        # mixing a gas with itself is the identity, at any ratio
        for m in (0.0, uni(0.1, 10.0), 1.0e6)
            g = mix(gA, gA, m)
            @test IdealGasThermo.cp(g, 600.0) ≈ IdealGasThermo.cp(air, 600.0) rtol = 1e-12
            @test IdealGasThermo.s0(g, 600.0) ≈ IdealGasThermo.s0(air, 600.0) rtol = 1e-12
        end
    end

    @testset "entropy of mixing: blend s0 vs a hand-written −R·Σ X·ln X" begin
        # A sign flip in the entropy-of-mixing fold (species.jl: alow[end] -=
        # Δs_mix) is invisible to every pure-species check (X·ln X = 0 at X = 1)
        # and to the legacy-agreement tests (they share the same fold). Pin it
        # directly: the molar s0 of a blend must equal Σ Xᵢ·s0ᵢ(T) − R·Σ Xᵢ·ln Xᵢ,
        # with the mixing term written out HERE, not taken from
        # generate_composite_species. A flipped sign would shift this by
        # 2R·|Σ X ln X| (~11.5 J/mol/K for a 50/50 blend). (test audit 2026-06-17)
        R = IdealGasThermo.Runiv
        for X in ([("N2", 0.5), ("O2", 0.5)],
                  [("N2", 0.78), ("O2", 0.21), ("CO2", 0.01)])
            xs = last.(X)
            sps = [species_in_spdict(first(p)) for p in X]
            comps = [FrozenGas(sp) for sp in sps]
            blend = FrozenGas(IdealGasThermo.Xidict2Array(Dict(X)))
            Δs_mix = sum(x * log(x) for x in xs)
            for T in (300.0, 1000.0, 2200.0)
                s0_blend = IdealGasThermo.s0(blend, T) * blend.MW / 1000   # molar
                expected = sum(x * IdealGasThermo.s0(c, T) * sp.MW / 1000
                               for (x, c, sp) in zip(xs, comps, sps)) - R * Δs_mix
                @test s0_blend ≈ expected rtol = 1e-10
            end
        end
    end

    @testset "coefficient self-consistency: s0 and h against ∫cp" begin
        # metamorphic check independent of any reference implementation:
        # the h and s0 polynomials must be the integrals of the cp
        # polynomial (composite Simpson, fine enough that quadrature error
        # is far below the tolerance; one interval straddles the 1000 K
        # coefficient breakpoint by construction)
        simpson(f, a, b, n) =
            (b - a) / (3n) * sum(
                (i == 0 || i == n ? 1 : iseven(i) ? 2 : 4) * f(a + (b - a) * i / n)
                for i = 0:n
            )
        cases = [(air, 300.0, 800.0), (products(sysJet, 0.03), 700.0, 1500.0)]
        push!(cases, (gas_pool()[end], uni(250.0, 900.0), uni(1100.0, 2200.0)))
        for (gas, T1, T2) in cases
            ∫cp = simpson(T -> IdealGasThermo.cp(gas, T), T1, T2, 2048)
            ∫cpT = simpson(T -> IdealGasThermo.cp(gas, T) / T, T1, T2, 2048)
            @test IdealGasThermo.h(gas, T2) - IdealGasThermo.h(gas, T1) ≈ ∫cp rtol = 1e-8
            @test IdealGasThermo.s0(gas, T2) - IdealGasThermo.s0(gas, T1) ≈ ∫cpT rtol = 1e-8
        end
    end

    @testset "error paths: pathological inputs fail loudly, not silently" begin
        # inversion driven into negative-temperature territory: the Newton
        # iterate hits log(T ≤ 0)
        @test_throws DomainError T_from_h(air, -1.0e10)
        # inversion target beyond any representable temperature: the bounded
        # Newton loop exhausts NEWTON_MAXITER and raises (the non-convergence
        # branch is live code, not decoration)
        @test_throws ErrorException T_from_h(air, 1.0e12)
        # beyond-stoichiometric FAR would need negative O2 — errors loudly
        # rather than returning an unphysical composition
        @test_throws DomainError products(sysCH4, 0.5)
    end

    @testset "Dual-carrying gas: AD inversions match finite differences" begin
        # products(sys, FAR::Dual) yields a FrozenGas{<:Dual} whose coefficients
        # carry the FAR-tangent, so inverting through it exercises the full
        # three-term IFT rule (the "composition moves" term). The invariant that
        # actually defines the rule's correctness — independent of fuel, oxidizer
        # composition, and the underlying NASA-9 data — is AD == central FD. (A
        # single pinned physics value like dT4/dFAR would instead track
        # thermo.inp and the air composition, not the rule.) We also assert the
        # result is a real number, not the nested-Dual tree the constant-
        # substance rule produced for a Dual-carrying gas.
        D = ForwardDiff.derivative
        δ = 1e-6
        fdcheck(f, x) = (f(x + δ) - f(x - δ)) / (2δ)
        for fuel in ["CH4", "H2", "Jet-A(g)"]   # distinct C/H/O ⟹ distinct ∂X/∂FAR
            sys = Vitiator(fuel, DryAir)
            fsp = species_in_spdict(fuel)
            hF = 1000 * fsp.Hf / fsp.MW
            # Stay sub-stoichiometric for the leanest-ceiling fuel in the pool
            # (H2 reaches stoichiometric at FAR ≈ 0.029; beyond it `products`
            # correctly throws on negative O₂).
            FAR0 = uni(0.005, 0.02)

            # (1) burner energy balance: Dual gas + Dual target
            hA = IdealGasThermo.h(air, uni(500.0, 900.0))
            T4of(far) = T_from_h(products(sys, far), (hA + far * hF) / (1 + far))
            d1 = D(T4of, FAR0)
            @test d1 isa Real                       # not a nested Dual
            @test d1 ≈ fdcheck(T4of, FAR0) rtol = 1e-6

            # (2) Dual gas + plain-Float target: invert a moving gas to fixed h
            hfix = IdealGasThermo.h(products(sys, FAR0), uni(1200.0, 1900.0))
            Tfix(far) = T_from_h(products(sys, far), hfix)
            @test D(Tfix, FAR0) ≈ fdcheck(Tfix, FAR0) rtol = 1e-6

            # (3) Dual gas through the isentropic inversion
            T1, PR = uni(1100.0, 1500.0), uni(2.0, 10.0)
            Tisen(far) = IdealGasThermo._T_polytropic(products(sys, far), T1, PR)
            @test D(Tisen, FAR0) ≈ fdcheck(Tisen, FAR0) rtol = 1e-6
        end
    end

    @testset "Dual-carrying gas: AD forward properties match finite differences" begin
        # products(sys, FAR::Dual) yields a FrozenGas{<:Dual} whose coefficients
        # carry the FAR-tangent. When T is also a same-tag Dual, both composition
        # and temperature move with one seed and the total derivative is the sum
        # of the composition tangent and the temperature tangent. As for the
        # inversions, the fuel/oxidizer-independent invariant is AD == central FD;
        # seeding a single ξ through both products(sys, FAR0+ξ) and T0+ξ is exactly
        # the same-tag Dual-gas-at-Dual-T path PowerCycles exercises.
        D = ForwardDiff.derivative
        δ = 1e-6
        fdcheck(f) = (f(δ) - f(-δ)) / (2δ)
        for fuel in ["CH4", "H2", "Jet-A(g)"]   # distinct C/H/O ⟹ distinct ∂X/∂FAR
            sys  = Vitiator(fuel, DryAir)
            FAR0 = uni(0.005, 0.02)             # sub-stoichiometric (H2 ceiling ≈ 0.029)
            T0   = uni(900.0, 1800.0)
            gasf(ξ) = products(sys, FAR0 + ξ)   # composition moves with the seed…
            # …and T with it. h/cp/s0/props are the new rules; gamma/speed_of_sound/
            # pressure_ratio carry none of their own and must inherit the total
            # derivative through generic dispatch.
            funcs = (
                ξ -> IdealGasThermo.h(gasf(ξ), T0 + ξ),
                ξ -> IdealGasThermo.cp(gasf(ξ), T0 + ξ),
                ξ -> IdealGasThermo.s0(gasf(ξ), T0 + ξ),
                ξ -> props(gasf(ξ), T0 + ξ).h,
                ξ -> props(gasf(ξ), T0 + ξ).cp,
                ξ -> props(gasf(ξ), T0 + ξ).s0,
                ξ -> IdealGasThermo.gamma(gasf(ξ), T0 + ξ),
                ξ -> IdealGasThermo.speed_of_sound(gasf(ξ), T0 + ξ),
                ξ -> IdealGasThermo.pressure_ratio(gasf(ξ), T0 + ξ, T0 + 200.0 + 2ξ),
            )
            for f in funcs
                @test D(f, 0.0) ≈ fdcheck(f) rtol = 1e-6
            end
        end

        # Single-layer return type. The nested same-tag Dual the constant-
        # substance rule produces for a Dual-carrying gas is what breaks downstream
        # Jacobian assembly — and the value-only ≈ check above can pass on it — so
        # pin the concrete type directly. It is fuel-invariant, so assert it once.
        tag    = :fwdtype
        gasd   = products(Vitiator("CH4", DryAir), ForwardDiff.Dual{tag}(0.012, 1.0))
        Td     = ForwardDiff.Dual{tag}(1500.0, 1.0)
        prd    = props(gasd, Td)
        single = ForwardDiff.Dual{tag,Float64,1}
        @test IdealGasThermo.h(gasd, Td)              isa single
        @test IdealGasThermo.cp(gasd, Td)             isa single
        @test IdealGasThermo.s0(gasd, Td)             isa single
        @test prd.h  isa single
        @test prd.cp isa single
        @test prd.s0 isa single
        @test IdealGasThermo.gamma(gasd, Td)          isa single   # derived, self-heals
        @test IdealGasThermo.speed_of_sound(gasd, Td) isa single   # derived, self-heals
        @test IdealGasThermo.pressure_ratio(gasd, Td, ForwardDiff.Dual{tag}(1700.0, 1.0)) isa single

        # A different-tag temperature is a legitimate nested (higher-order) AD,
        # not the same-tag bug — the new rules must leave it nested: the gas has
        # tag :inner and h(gas, value(T)) is itself a Dual{:inner}, so the value
        # rail of the result is correctly a Dual{:inner}.
        gas_inner = products(Vitiator("CH4", DryAir), ForwardDiff.Dual{:inner}(0.03, 1.0))
        T_outer   = ForwardDiff.Dual{:outer}(1600.0, 1.0)
        @test ForwardDiff.value(IdealGasThermo.h(gas_inner, T_outer)) isa ForwardDiff.Dual{:inner}
    end

end
