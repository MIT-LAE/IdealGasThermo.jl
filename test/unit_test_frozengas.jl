using ForwardDiff

@testset "FrozenGas pure core" begin

    @testset "construction from composite species" begin
        air = FrozenGas(DryAir)
        # isbits ⟹ stack-allocated, usable in immutable params NamedTuples,
        # thread-safe by construction
        @test isbitstype(typeof(air))
        # Dry air specific gas constant [J/kg/K]
        @test IdealGasThermo.R(air) ≈ 287.05 atol = 0.1
    end

    @testset "absolute values, datum, composite construction" begin
        air = FrozenGas(DryAir)
        # cp/h/s0 ABSOLUTE values are anchored to CEA (an external authority) in
        # unit_test_cea_reference.jl — per-species and the dry-air pseudo-species.
        # Here we keep only the claims that file does not make; the old
        # `cp ≈ 1005 atol=5` textbook ballpark and the `≈ Gas1D` agreement loop
        # (a same-NASA-9-kernel twin) are gone — they added a loose/circular
        # echo of the now-tight CEA anchor.

        # γ ≈ 1.4 for diatomic-dominated air near room T (sanity; γ = cp/(cp−R),
        # and cp/R is CEA-anchored)
        @test IdealGasThermo.gamma(air, 300.0) ≈ 1.400 atol = 0.003

        # formation-inclusive datum: h(298.15) equals the mixture mass-specific
        # formation enthalpy, so sensible enthalpy is h(T) − h(298.15)
        @test IdealGasThermo.h(air, 298.15) ≈ air.Hf / air.MW * 1000.0 rtol = 1e-4

        # the DryAir composite reproduces the fitted "Air" pseudo-species — the
        # gas anchored against the CEA Air table — tying DryAir to that external
        # anchor without re-listing the values (composite path vs single fit)
        airsp = FrozenGas(species_in_spdict("Air"))
        for T in (300.0, 1000.0, 1600.0)
            @test IdealGasThermo.cp(air, T) ≈ IdealGasThermo.cp(airsp, T) rtol = 1e-4
            @test IdealGasThermo.s0(air, T) ≈ IdealGasThermo.s0(airsp, T) rtol = 1e-4
        end
    end

    @testset "props" begin
        air = FrozenGas(DryAir)
        for T in [250.0, 999.9, 1000.0, 1600.0]
            p = props(air, T)
            @test p.cp == IdealGasThermo.cp(air, T)
            @test p.h ≈ IdealGasThermo.h(air, T) rtol = 1e-14
            @test p.s0 ≈ IdealGasThermo.s0(air, T) rtol = 1e-14
        end
    end

    @testset "exported property API (cₚ/c_p aliases + unqualified access)" begin
        # The pure-core property accessors are exported, so a consumer that did
        # `using IdealGasThermo` can call them WITHOUT qualifying. These bare names
        # (no `IdealGasThermo.` prefix) only resolve if the export is in place — a
        # regression in the export list makes this testset fail to even run.
        air = FrozenGas(DryAir)
        T = 800.0
        # cₚ and c_p are aliases of the same (unexported, Base.cp-conflicting) `cp`
        @test cₚ === IdealGasThermo.cp
        @test c_p === IdealGasThermo.cp
        @test cₚ(air, T) == IdealGasThermo.cp(air, T)
        @test c_p(air, T) === cₚ(air, T)
        # the rest of the set resolves unqualified and agrees with the qualified call
        @test h(air, T) == IdealGasThermo.h(air, T)
        @test s0(air, T) == IdealGasThermo.s0(air, T)
        @test gamma(air, T) == IdealGasThermo.gamma(air, T)
        @test γ === gamma                       # Unicode alias of gamma
        @test γ(air, T) == gamma(air, T)
        @test R(air) == IdealGasThermo.R(air)
        @test T_of_h(air, h(air, T)) ≈ T rtol = 1e-10
        T2 = T_isentropic(air, 300.0, 8.0)
        @test pressure_ratio(air, 300.0, T2) ≈ 8.0 rtol = 1e-10
        # the props bundle deliberately keeps the field name `cp` (a field accessor
        # never collides with Base.cp), and it matches the cₚ alias
        @test props(air, T).cp == cₚ(air, T)
    end

    @testset "T_of_h inversion" begin
        air = FrozenGas(DryAir)
        for T in [250.0, 500.0, 999.5, 1000.5, 1600.0, 2200.0]
            @test IdealGasThermo.T_of_h(air, IdealGasThermo.h(air, T)) ≈ T rtol = 1e-10
        end
        # deterministic: same input, same output, bitwise
        hspec = IdealGasThermo.h(air, 1234.5)
        @test IdealGasThermo.T_of_h(air, hspec) === IdealGasThermo.T_of_h(air, hspec)
    end

    @testset "isentropic relations" begin
        air = FrozenGas(DryAir)
        # identity: no pressure change, no temperature change
        @test IdealGasThermo.T_isentropic(air, 500.0, 1.0) ≈ 500.0 rtol = 1e-12
        # round trip: s0(T2) = s0(T1) + R ln(PR)  ⟺  pressure_ratio inverts it
        for (T1, PR) in [(288.15, 12.0), (500.0, 30.0), (1600.0, 0.25), (900.0, 1.05)]
            T2 = IdealGasThermo.T_isentropic(air, T1, PR)
            @test IdealGasThermo.pressure_ratio(air, T1, T2) ≈ PR rtol = 1e-10
        end
    end

    @testset "zero allocations after warmup" begin
        air = FrozenGas(DryAir)
        measured(f::F, args...) where {F} = (f(args...); @allocated f(args...))
        @test measured(IdealGasThermo.cp, air, 600.0) == 0
        @test measured(IdealGasThermo.h, air, 600.0) == 0
        @test measured(IdealGasThermo.s0, air, 600.0) == 0
        @test measured(IdealGasThermo.gamma, air, 600.0) == 0
        @test measured(props, air, 600.0) == 0
        @test measured(IdealGasThermo.T_of_h, air, 5e5) == 0
        @test measured(IdealGasThermo.T_isentropic, air, 288.15, 12.0) == 0
        @test measured(IdealGasThermo.pressure_ratio, air, 288.15, 600.0) == 0
    end

    @testset "derivatives: ForwardDiff vs analytic" begin
        air = FrozenGas(DryAir)
        D = ForwardDiff.derivative
        for T in [250.0, 500.0, 999.0, 1001.0, 1600.0, 2200.0]
            # closed forms: dh/dT = cp, ds0/dT = cp/T
            @test D(t -> IdealGasThermo.h(air, t), T) ≈
                  IdealGasThermo.cp(air, T) rtol = 1e-10
            @test D(t -> IdealGasThermo.s0(air, t), T) ≈
                  IdealGasThermo.cp(air, T) / T rtol = 1e-10
            # props must carry the same derivatives as the scalar functions
            @test D(t -> props(air, t).h, T) ≈ IdealGasThermo.cp(air, T) rtol = 1e-10
            # dcp/dT against central finite difference
            dcp_fd = (IdealGasThermo.cp(air, T + 0.5) - IdealGasThermo.cp(air, T - 0.5))
            @test D(t -> IdealGasThermo.cp(air, t), T) ≈ dcp_fd rtol = 1e-4
        end
        # inversion derivatives via implicit function theorem: dT/dh = 1/cp
        for T in [400.0, 1500.0]
            hT = IdealGasThermo.h(air, T)
            @test D(hh -> IdealGasThermo.T_of_h(air, hh), hT) ≈
                  1 / IdealGasThermo.cp(air, T) rtol = 1e-10
        end
        # T_isentropic: ∂T2/∂PR = R·T2 / (ηp·PR·cp(T2)) from s0(T2) = s0(T1) + R ln(PR)/ηp
        T1, PR = 288.15, 12.0
        T2 = IdealGasThermo.T_isentropic(air, T1, PR)
        @test D(pr -> IdealGasThermo.T_isentropic(air, T1, pr), PR) ≈
              air.R * T2 / (PR * IdealGasThermo.cp(air, T2)) rtol = 1e-10
        # ∂T2/∂T1 = cp(T1)·T2 / (cp(T2)·T1)
        @test D(t1 -> IdealGasThermo.T_isentropic(air, t1, PR), T1) ≈
              IdealGasThermo.cp(air, T1) * T2 / (IdealGasThermo.cp(air, T2) * T1) rtol = 1e-10
    end

    @testset "ForwardDiff extension (analytic fast paths)" begin
        # the analytic-rule extension is active, not just generic fallback
        @test Base.get_extension(IdealGasThermo, :IdealGasThermoForwardDiffExt) !== nothing
        air = FrozenGas(DryAir)
        D = ForwardDiff.derivative
        measured(f::F, args...) where {F} = (f(args...); @allocated f(args...))
        # Dual-typed evaluation stays allocation-free through the rules
        @test measured(D, t -> IdealGasThermo.h(air, t), 600.0) == 0
        @test measured(D, t -> props(air, t).h, 600.0) == 0
        @test measured(D, hh -> IdealGasThermo.T_of_h(air, hh), 5e5) == 0
        @test measured(D, pr -> IdealGasThermo.T_isentropic(air, 288.15, pr), 12.0) == 0
        # nested duals: d²h/dT² == dcp/dT
        d2h = D(t -> D(s -> IdealGasThermo.h(air, s), t), 1600.0)
        @test d2h ≈ D(t -> IdealGasThermo.cp(air, t), 1600.0) rtol = 1e-10
    end

    @testset "construction from mole fractions" begin
        X = IdealGasThermo.Xidict2Array(IdealGasThermo.Xair)
        X = X ./ sum(X)
        air_from_X = FrozenGas(X)
        air = FrozenGas(DryAir)
        @test IdealGasThermo.R(air_from_X) ≈ IdealGasThermo.R(air) rtol = 1e-12
        for T in [300.0, 1600.0]
            @test IdealGasThermo.cp(air_from_X, T) ≈ IdealGasThermo.cp(air, T) rtol = 1e-12
            @test IdealGasThermo.h(air_from_X, T) ≈ IdealGasThermo.h(air, T) rtol = 1e-12
        end
    end

end
