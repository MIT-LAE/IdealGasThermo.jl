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

    @testset "cp" begin
        air = FrozenGas(DryAir)
        # CEA-derived acceptance value for dry air
        @test IdealGasThermo.cp(air, 300.0) ≈ 1005.0 atol = 5.0
        # Agreement with the established Gas1D mixture thermodynamics
        g1d = Gas1D()
        for T in [250.0, 500.0, 999.9, 1000.0, 1600.0, 2200.0]
            g1d.T = T
            @test IdealGasThermo.cp(air, T) ≈ g1d.cp rtol = 1e-10
        end
    end

    @testset "h, s0, gamma" begin
        air = FrozenGas(DryAir)
        @test IdealGasThermo.gamma(air, 300.0) ≈ 1.400 atol = 0.003
        g1d = Gas1D()
        for T in [250.0, 500.0, 999.9, 1000.0, 1600.0, 2200.0]
            g1d.T = T
            @test IdealGasThermo.h(air, T) ≈ g1d.h rtol = 1e-10
            @test IdealGasThermo.s0(air, T) ≈ g1d.ϕ rtol = 1e-10
        end
        # Enthalpy datum is formation-inclusive: h at 298.15 K equals the
        # mixture mass-specific formation enthalpy, so sensible enthalpy is
        # h(T) - h(298.15)
        @test IdealGasThermo.h(air, 298.15) ≈ air.Hf / air.MW * 1000.0 rtol = 1e-4
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

    @testset "T_of_h inversion" begin
        air = FrozenGas(DryAir)
        for T in [250.0, 500.0, 999.5, 1000.5, 1600.0, 2200.0]
            @test T_of_h(air, IdealGasThermo.h(air, T)) ≈ T rtol = 1e-10
        end
        # deterministic: same input, same output, bitwise
        hspec = IdealGasThermo.h(air, 1234.5)
        @test T_of_h(air, hspec) === T_of_h(air, hspec)
    end

    @testset "isentropic relations" begin
        air = FrozenGas(DryAir)
        # identity: no pressure change, no temperature change
        @test T_isentropic(air, 500.0, 1.0) ≈ 500.0 rtol = 1e-12
        # round trip: s0(T2) = s0(T1) + R ln(PR)  ⟺  pressure_ratio inverts it
        for (T1, PR) in [(288.15, 12.0), (500.0, 30.0), (1600.0, 0.25), (900.0, 1.05)]
            T2 = T_isentropic(air, T1, PR)
            @test IdealGasThermo.pressure_ratio(air, T1, T2) ≈ PR rtol = 1e-10
        end
        # cross-validation against the established compress (ideal, ηp = 1)
        g1d = Gas1D()
        set_TP!(g1d, 288.15, 101325.0)
        IdealGasThermo.compress(g1d, 12.0)
        @test T_isentropic(air, 288.15, 12.0) ≈ g1d.T rtol = 1e-7
    end

    @testset "zero allocations after warmup" begin
        air = FrozenGas(DryAir)
        measured(f::F, args...) where {F} = (f(args...); @allocated f(args...))
        @test measured(IdealGasThermo.cp, air, 600.0) == 0
        @test measured(IdealGasThermo.h, air, 600.0) == 0
        @test measured(IdealGasThermo.s0, air, 600.0) == 0
        @test measured(IdealGasThermo.gamma, air, 600.0) == 0
        @test measured(props, air, 600.0) == 0
        @test measured(T_of_h, air, 5e5) == 0
        @test measured(T_isentropic, air, 288.15, 12.0) == 0
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
            @test D(hh -> T_of_h(air, hh), hT) ≈
                  1 / IdealGasThermo.cp(air, T) rtol = 1e-10
        end
        # T_isentropic: ∂T2/∂PR = R·T2 / (ηp·PR·cp(T2)) from s0(T2) = s0(T1) + R ln(PR)/ηp
        T1, PR = 288.15, 12.0
        T2 = T_isentropic(air, T1, PR)
        @test D(pr -> T_isentropic(air, T1, pr), PR) ≈
              air.R * T2 / (PR * IdealGasThermo.cp(air, T2)) rtol = 1e-10
        # ∂T2/∂T1 = cp(T1)·T2 / (cp(T2)·T1)
        @test D(t1 -> T_isentropic(air, t1, PR), T1) ≈
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
        @test measured(D, hh -> T_of_h(air, hh), 5e5) == 0
        @test measured(D, pr -> T_isentropic(air, 288.15, pr), 12.0) == 0
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
