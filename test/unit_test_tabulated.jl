using ForwardDiff

@testset "TabulatedGas fast inversions" begin

    @testset "construction and seeded T_of_h" begin
        air = FrozenGas(DryAir)
        tg = tabulate(air)
        @test tg isa TabulatedGas
        # default-tier inversion is exact: table-seeded Newton against the
        # same polynomials, same convergence criterion as FrozenGas
        for T in [250.0, 500.0, 999.0, 1000.5, 1600.0, 2200.0]
            @test T_of_h(tg, IdealGasThermo.h(air, T)) ≈ T rtol = 1e-10
        end
        # deterministic: same input, same output, bitwise
        hspec = IdealGasThermo.h(air, 1234.5)
        @test T_of_h(tg, hspec) === T_of_h(tg, hspec)
    end

    @testset "seeded T_isentropic" begin
        air = FrozenGas(DryAir)
        tg = tabulate(air)
        # same answers as the FrozenGas Newton solve, compression and expansion
        for (T1, PR) in [(288.15, 12.0), (500.0, 30.0), (1600.0, 0.25), (900.0, 1.05)]
            @test T_isentropic(tg, T1, PR) ≈ T_isentropic(air, T1, PR) rtol = 1e-10
            @test T_isentropic(tg, T1, PR; ηp = 0.9) ≈
                  T_isentropic(air, T1, PR; ηp = 0.9) rtol = 1e-10
        end
        # deterministic
        @test T_isentropic(tg, 288.15, 12.0) === T_isentropic(tg, 288.15, 12.0)
        # integer arguments are accepted, like the FrozenGas methods
        @test T_isentropic(tg, 288.15, 12) === T_isentropic(tg, 288.15, 12.0)
        @test T_of_h(tg, 500000) === T_of_h(tg, 500000.0)
    end

    @testset "out-of-range targets fall back to exact FrozenGas solve" begin
        air = FrozenGas(DryAir)
        tg = tabulate(air) # default range T ∈ [200, 2400]
        # h above the table ceiling (T = 2600 K > Tmax = 2400 K): no error,
        # no extrapolation — the exact cold-start Newton answer
        hhot = IdealGasThermo.h(air, 2600.0)
        @test hhot > tg.hmax # the target really is outside the table
        @test T_of_h(tg, hhot) ≈ T_of_h(air, hhot) rtol = 1e-10
        # h below the table floor (T = 150 K < Tmin = 200 K)
        hcold = IdealGasThermo.h(air, 150.0)
        @test hcold < tg.hmin
        @test T_of_h(tg, hcold) ≈ T_of_h(air, hcold) rtol = 1e-10
        # the prototype blowout case: T1 = 1500 K, PR = 40 lands above s0max
        T1, PR = 1500.0, 40.0
        @test IdealGasThermo.s0(air, T1) + air.R * log(PR) > tg.s0max
        @test T_isentropic(tg, T1, PR) ≈ T_isentropic(air, T1, PR) rtol = 1e-10
    end

    @testset "interp tier: pure table lookup" begin
        air = FrozenGas(DryAir)
        tg = tabulate(air)
        # documented accuracy |ΔT/T| ≤ 1e-9 at N = 256 (measured 5.8e-10 for
        # h, 1.4e-9 for s0); test at 2e-9 across a dense in-range sweep
        for T in range(250.0, 2200.0; length = 1000)
            @test T_of_h_interp(tg, IdealGasThermo.h(air, T)) ≈ T rtol = 2e-9
        end
        # isentropic interp agrees with the exact solve to the same tier accuracy
        for (T1, PR) in [(288.15, 12.0), (500.0, 30.0), (1600.0, 0.25), (900.0, 1.05)]
            @test T_isentropic_interp(tg, T1, PR) ≈
                  T_isentropic(air, T1, PR) rtol = 2e-9
            @test T_isentropic_interp(tg, T1, PR; ηp = 0.9) ≈
                  T_isentropic(air, T1, PR; ηp = 0.9) rtol = 2e-9
        end
        # out of range: the approximate tier is loud — DomainError, no fallback
        @test_throws DomainError T_of_h_interp(tg, IdealGasThermo.h(air, 2600.0))
        @test_throws DomainError T_of_h_interp(tg, IdealGasThermo.h(air, 150.0))
        @test_throws DomainError T_isentropic_interp(tg, 1500.0, 40.0)
    end

    @testset "forward properties forward exactly to the wrapped gas" begin
        air = FrozenGas(DryAir)
        tg = tabulate(air)
        @test IdealGasThermo.R(tg) === IdealGasThermo.R(air)
        for T in [250.0, 999.9, 1000.0, 1600.0, 2600.0] # incl. outside table range
            @test IdealGasThermo.cp(tg, T) === IdealGasThermo.cp(air, T)
            @test IdealGasThermo.h(tg, T) === IdealGasThermo.h(air, T)
            @test IdealGasThermo.s0(tg, T) === IdealGasThermo.s0(air, T)
            @test IdealGasThermo.gamma(tg, T) === IdealGasThermo.gamma(air, T)
            @test props(tg, T) === props(air, T)
        end
        @test IdealGasThermo.pressure_ratio(tg, 288.15, 600.0) ===
              IdealGasThermo.pressure_ratio(air, 288.15, 600.0)
    end

    @testset "zero allocations after warmup" begin
        air = FrozenGas(DryAir)
        tg = tabulate(air) # construction allocates (the tables); calls must not
        measured(f::F, args...) where {F} = (f(args...); @allocated f(args...))
        hin = IdealGasThermo.h(air, 700.0) # in-range target
        @test measured(T_of_h, tg, hin) == 0
        @test measured(T_isentropic, tg, 288.15, 12.0) == 0
        @test measured(T_of_h_interp, tg, hin) == 0
        @test measured(T_isentropic_interp, tg, 288.15, 12.0) == 0
        # the out-of-range fallback path is allocation-free too
        @test measured(T_of_h, tg, IdealGasThermo.h(air, 2600.0)) == 0
        @test measured(T_isentropic, tg, 1500.0, 40.0) == 0
    end

    @testset "ForwardDiff through seeded inversions (IFT rules)" begin
        @test Base.get_extension(IdealGasThermo, :IdealGasThermoForwardDiffExt) !== nothing
        air = FrozenGas(DryAir)
        tg = tabulate(air)
        D = ForwardDiff.derivative
        # dT/dh = 1/cp at the solution (implicit function theorem)
        for T in [400.0, 1500.0]
            hT = IdealGasThermo.h(air, T)
            @test D(hh -> T_of_h(tg, hh), hT) ≈
                  1 / IdealGasThermo.cp(air, T) rtol = 1e-10
        end
        # ∂T2/∂PR = R·T2 / (PR·cp(T2)) from s0(T2) = s0(T1) + R ln(PR)
        T1, PR = 288.15, 12.0
        T2 = T_isentropic(tg, T1, PR)
        @test D(pr -> T_isentropic(tg, T1, pr), PR) ≈
              air.R * T2 / (PR * IdealGasThermo.cp(air, T2)) rtol = 1e-10
        # ∂T2/∂T1 = cp(T1)·T2 / (cp(T2)·T1)
        @test D(t1 -> T_isentropic(tg, t1, PR), T1) ≈
              IdealGasThermo.cp(air, T1) * T2 /
              (IdealGasThermo.cp(air, T2) * T1) rtol = 1e-10
        # both arguments Dual (gradient): sum of the two partials
        g = ForwardDiff.gradient(x -> T_isentropic(tg, x[1], x[2]), [T1, PR])
        @test g[1] ≈ D(t1 -> T_isentropic(air, t1, PR), T1) rtol = 1e-10
        @test g[2] ≈ D(pr -> T_isentropic(air, T1, pr), PR) rtol = 1e-10
        # Dual-typed evaluation stays allocation-free through the rules
        measured(f::F, args...) where {F} = (f(args...); @allocated f(args...))
        @test measured(D, hh -> T_of_h(tg, hh), IdealGasThermo.h(air, 700.0)) == 0
        @test measured(D, pr -> T_isentropic(tg, 288.15, pr), 12.0) == 0
    end

end
