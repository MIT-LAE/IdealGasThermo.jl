using ForwardDiff

@testset "FastFrozenGas inversions" begin

    # the public inversion verb, used throughout: the keyword names what is
    # known; the isentrope is a *process* and lives in compress/expand
    Th(g, x) = temperature(g, h = x)
    # single-method positional helpers: local closures with kwarg defaults
    # have non-specializing kwsorters, and multi-method local functions get
    # boxed when captured — both would charge allocations to the test helper
    # rather than the package verbs (the verbs themselves are verified 0-alloc)
    Cmp(g, T1, PR) = compress(g, T1, PR)
    Cmp_eta(g, T1, PR, etap) = compress(g, T1, PR; ηp = etap)
    Exp(g, T1, PR) = expand(g, T1, PR)
    Exp_eta(g, T1, PR, etap) = expand(g, T1, PR; ηp = etap)

    @testset "construction and :seeded temperature(h = ...)" begin
        air = FrozenGas(DryAir)
        fg = FastFrozenGas(air)
        @test fg isa FastFrozenGas{:seeded} # exact tier is the default
        # default-mode inversion is exact: table-seeded Newton against the
        # same polynomials, same convergence criterion as FrozenGas
        for T in [250.0, 500.0, 999.0, 1000.5, 1600.0, 2200.0]
            @test Th(fg, IdealGasThermo.h(air, T)) ≈ T rtol = 1e-10
        end
        hspec = IdealGasThermo.h(air, 1234.5)
        # the same verb works on the plain gas
        @test Th(air, hspec) ≈ Th(fg, hspec) rtol = 1e-12
        # mode is validated
        @test_throws ArgumentError FastFrozenGas(air, mode = :warp)
    end

    @testset ":seeded compress/expand (the isentrope process verbs)" begin
        air = FrozenGas(DryAir)
        fg = FastFrozenGas(air)
        # same answers as the FrozenGas Newton solve; both verbs take ratios
        # ≥ 1 — the old PR = 0.25 isentrope case is expand with ratio 4
        for (T1, PR) in [(288.15, 12.0), (500.0, 30.0), (900.0, 1.05)]
            @test Cmp(fg, T1, PR) ≈ Cmp(air, T1, PR) rtol = 1e-10
            @test Cmp_eta(fg, T1, PR, 0.9) ≈ Cmp_eta(air, T1, PR, 0.9) rtol = 1e-10
        end
        for (T1, PR) in [(1600.0, 4.0), (900.0, 1.05)]
            @test Exp(fg, T1, PR) ≈ Exp(air, T1, PR) rtol = 1e-10
            @test Exp_eta(fg, T1, PR, 0.9) ≈ Exp_eta(air, T1, PR, 0.9) rtol = 1e-10
        end
        # integer arguments are accepted, like the FrozenGas methods
        @test Cmp(fg, 288.15, 12) === Cmp(fg, 288.15, 12.0)
        @test Th(fg, 500000) === Th(fg, 500000.0)
    end

    @testset ":seeded out-of-range targets fall back to exact solve" begin
        air = FrozenGas(DryAir)
        fg = FastFrozenGas(air) # default range T ∈ [200, 2400]
        # h above the table ceiling (T = 2600 K > Tmax = 2400 K): no error,
        # no extrapolation — the exact cold-start Newton answer
        hhot = IdealGasThermo.h(air, 2600.0)
        @test hhot > fg.hmax # the target really is outside the table
        @test Th(fg, hhot) ≈ Th(air, hhot) rtol = 1e-10
        # h below the table floor (T = 150 K < Tmin = 200 K)
        hcold = IdealGasThermo.h(air, 150.0)
        @test hcold < fg.hmin
        @test Th(fg, hcold) ≈ Th(air, hcold) rtol = 1e-10
        # the prototype blowout case: T1 = 1500 K, PR = 40 lands above s0max
        T1, PR = 1500.0, 40.0
        @test IdealGasThermo.s0(air, T1) + air.R * log(PR) > fg.s0max
        @test Cmp(fg, T1, PR) ≈ Cmp(air, T1, PR) rtol = 1e-10
    end

    @testset ":fast mode — pure table lookup, same verb" begin
        air = FrozenGas(DryAir)
        fg = FastFrozenGas(air, mode = :fast)
        @test fg isa FastFrozenGas{:fast}
        # documented accuracy |ΔT/T| ≲ 2e-9 at N = 256 (measured 5.8e-10 for
        # h, 1.4e-9 for s0); test at 2e-9 across a dense in-range sweep
        for T in range(250.0, 2200.0; length = 1000)
            @test Th(fg, IdealGasThermo.h(air, T)) ≈ T rtol = 2e-9
        end
        # the isentrope verbs agree with the exact solve to the same tier
        # accuracy (the old PR = 0.25 case is expand with ratio 4)
        for (T1, PR) in [(288.15, 12.0), (500.0, 30.0), (900.0, 1.05)]
            @test Cmp(fg, T1, PR) ≈ Cmp(air, T1, PR) rtol = 2e-9
            @test Cmp_eta(fg, T1, PR, 0.9) ≈ Cmp_eta(air, T1, PR, 0.9) rtol = 2e-9
        end
        for (T1, PR) in [(1600.0, 4.0), (900.0, 1.05)]
            @test Exp(fg, T1, PR) ≈ Exp(air, T1, PR) rtol = 2e-9
            @test Exp_eta(fg, T1, PR, 0.9) ≈ Exp_eta(air, T1, PR, 0.9) rtol = 2e-9
        end
        # out of range: the approximate mode is loud — DomainError, no fallback
        @test_throws DomainError Th(fg, IdealGasThermo.h(air, 2600.0))
        @test_throws DomainError Th(fg, IdealGasThermo.h(air, 150.0))
        @test_throws DomainError Cmp(fg, 1500.0, 40.0)
    end

    @testset "forward properties forward exactly to the wrapped gas" begin
        air = FrozenGas(DryAir)
        for fg in (FastFrozenGas(air), FastFrozenGas(air, mode = :fast))
            @test IdealGasThermo.R(fg) === IdealGasThermo.R(air)
            for T in [300.0, 1200.0] # two temperatures straddling the 1000 K seam
                @test IdealGasThermo.cp(fg, T) === IdealGasThermo.cp(air, T)
                @test IdealGasThermo.h(fg, T) === IdealGasThermo.h(air, T)
                @test IdealGasThermo.s0(fg, T) === IdealGasThermo.s0(air, T)
                @test IdealGasThermo.gamma(fg, T) === IdealGasThermo.gamma(air, T)
                @test props(fg, T) === props(air, T)
            end
            @test IdealGasThermo.pressure_ratio(fg, 288.15, 600.0) ===
                  IdealGasThermo.pressure_ratio(air, 288.15, 600.0)
        end
    end

    @testset "zero allocations after warmup" begin
        air = FrozenGas(DryAir)
        seeded = FastFrozenGas(air) # construction allocates (tables); calls must not
        fast = FastFrozenGas(air, mode = :fast)
        hin = IdealGasThermo.h(air, 700.0) # in-range target
        @test (@ballocated $Th($seeded, $hin) samples = 1 evals = 1) == 0
        @test (@ballocated $Cmp($seeded, 288.15, 12.0) samples = 1 evals = 1) == 0
        @test (@ballocated $Exp($seeded, 1600.0, 4.0) samples = 1 evals = 1) == 0
        @test (@ballocated $Th($fast, $hin) samples = 1 evals = 1) == 0
        @test (@ballocated $Cmp($fast, 288.15, 12.0) samples = 1 evals = 1) == 0
        # the out-of-range fallback path is allocation-free too
        @test (@ballocated $Th($seeded, IdealGasThermo.h($air, 2600.0)) samples = 1 evals = 1) == 0
        @test (@ballocated $Cmp($seeded, 1500.0, 40.0) samples = 1 evals = 1) == 0
    end

    @testset "ForwardDiff through the verb (IFT rules)" begin
        @test Base.get_extension(IdealGasThermo, :IdealGasThermoForwardDiffExt) !== nothing
        air = FrozenGas(DryAir)
        fg = FastFrozenGas(air)
        D = ForwardDiff.derivative
        # dT/dh = 1/cp at the solution (implicit function theorem)
        for T in [400.0, 1500.0]
            hT = IdealGasThermo.h(air, T)
            @test D(hh -> Th(fg, hh), hT) ≈ 1 / IdealGasThermo.cp(air, T) rtol = 1e-10
        end
        # ∂T2/∂PR = R·T2 / (PR·cp(T2)) from s0(T2) = s0(T1) + R ln(PR)
        T1, PR = 288.15, 12.0
        T2 = Cmp(fg, T1, PR)
        @test D(pr -> Cmp(fg, T1, pr), PR) ≈
              air.R * T2 / (PR * IdealGasThermo.cp(air, T2)) rtol = 1e-10
        # ∂T2/∂T1 = cp(T1)·T2 / (cp(T2)·T1)
        @test D(t1 -> Cmp(fg, t1, PR), T1) ≈
              IdealGasThermo.cp(air, T1) * T2 /
              (IdealGasThermo.cp(air, T2) * T1) rtol = 1e-10
        # both arguments Dual (gradient): sum of the two partials
        g = ForwardDiff.gradient(x -> Cmp(fg, x[1], x[2]), [T1, PR])
        @test g[1] ≈ D(t1 -> Cmp(air, t1, PR), T1) rtol = 1e-10
        @test g[2] ≈ D(pr -> Cmp(air, T1, pr), PR) rtol = 1e-10
        # the :fast mode carries the same IFT tangents (primal ≲ 2e-9)
        fast = FastFrozenGas(air, mode = :fast)
        hT = IdealGasThermo.h(air, 700.0)
        @test D(hh -> Th(fast, hh), hT) ≈
              1 / IdealGasThermo.cp(air, Th(air, hT)) rtol = 1e-8
        # Dual-typed evaluation stays allocation-free through the rules
        @test (@ballocated $D(hh -> $Th($fg, hh), IdealGasThermo.h($air, 700.0)) samples = 1 evals = 1) == 0
        @test (@ballocated $D(pr -> $Cmp($fg, 288.15, pr), 12.0) samples = 1 evals = 1) == 0
    end

end
