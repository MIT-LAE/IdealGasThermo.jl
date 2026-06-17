using ForwardDiff

@testset "gas dynamics: speed of sound, Mach, stagnation/static" begin

    air = FrozenGas(DryAir)

    @testset "speed_of_sound is √(γRT), pure in (gas, T)" begin
        # NB: not `a ≈ sqrt(gamma*R*T)` — that re-derives the value with
        # speed_of_sound's own expression (a tautology). γ, R, cp are anchored
        # against CEA in unit_test_cea_reference.jl; here we pin the one thing
        # that catches a units/scale error in speed_of_sound itself: the
        # absolute sea-level dry-air value (~340.3 m/s).
        @test speed_of_sound(air, 288.15) ≈ 340.3 atol = 0.1
        @test speed_of_sound(air, 1600.0) > speed_of_sound(air, 288.15)  # rises with T
        # FastFrozenGas forwards the property unchanged
        fg = FastFrozenGas(air)
        @test speed_of_sound(fg, 700.0) === speed_of_sound(air, 700.0)
        # GasState accessor reads at st.T, needs no pressure
        for (T, P) in ((288.15, 101325.0), (1600.0, 12.0e5))
            @test speed_of_sound(GasState(air, T, P)) === speed_of_sound(air, T)
        end
    end

    @testset "mach = V / a" begin
        st = GasState(air, 288.15, 101325.0)
        a = speed_of_sound(air, 288.15)
        @test mach(air, 288.15, a) ≈ 1.0 rtol = 1e-14   # V = a ⇒ M = 1
        @test mach(st, 0.5 * a) ≈ 0.5 rtol = 1e-14
        @test mach(st, 100.0) === 100.0 / speed_of_sound(st)
    end

    @testset "stagnation: energy balance and isentropic pressure" begin
        st = GasState(air, 250.0, 5.0e4)
        for M in (0.0, 0.3, 0.8, 1.0, 2.0)
            stt = stagnation_state(st, M)
            V = M * speed_of_sound(st)
            # total enthalpy = static enthalpy + ½V² (energy conservation)
            @test IdealGasThermo.h(stt) ≈ IdealGasThermo.h(st) + V^2 / 2 rtol = 1e-12
            # the bring-to-rest is isentropic: full entropy unchanged
            @test entropy(stt) ≈ entropy(st) rtol = 1e-12
            # stagnation never lowers T or P (equality at M = 0, to within the
            # enthalpy-inversion tolerance)
            @test stt.T ≥ st.T - 1e-9
            @test stt.P ≥ st.P - 1e-6
        end
        @test stagnation_state(st, 0.0).T ≈ st.T rtol = 1e-12
        @test stagnation_state(st, 0.0).P ≈ st.P rtol = 1e-12
        @test_throws ArgumentError stagnation_state(st, -0.1)
    end

    @testset "static inverts stagnation across the Mach range" begin
        st = GasState(air, 288.15, 101325.0)
        for M in (0.0, 0.2, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0)
            back = static_state(stagnation_state(st, M), M)
            @test back.T ≈ st.T rtol = 1e-10
            @test back.P ≈ st.P rtol = 1e-9
        end
        @test static_state(st, 0.0).T ≈ st.T rtol = 1e-12   # M = 0 ⇒ static == total
        @test static_state(st, 0.0).P ≈ st.P rtol = 1e-12
        @test_throws ArgumentError static_state(st, -0.1)
        # static lowers T and P below the total state for M > 0
        sts = static_state(st, 0.9)
        @test sts.T < st.T
        @test sts.P < st.P
    end

    @testset "low-Mach limit matches the constant-γ ratio" begin
        # the exact enthalpy/entropy result reduces to 1 + ½(γ−1)M² + O(M⁴);
        # at M = 0.1 the constant-γ relation is good to ~1e-3 relative
        st = GasState(air, 300.0, 1.0e5)
        γ = IdealGasThermo.gamma(air, 300.0)
        for M in (0.05, 0.1)
            stt = stagnation_state(st, M)
            TR_cg = 1 + 0.5 * (γ - 1) * M^2
            PR_cg = TR_cg^(γ / (γ - 1))
            @test stt.T / st.T ≈ TR_cg rtol = 2e-3
            @test stt.P / st.P ≈ PR_cg rtol = 2e-3
        end
    end

    @testset "works for FastFrozenGas too" begin
        fg = FastFrozenGas(air)
        stf = GasState(fg, 288.15, 101325.0)
        st = GasState(air, 288.15, 101325.0)
        @test stagnation_state(stf, 0.8).T ≈ stagnation_state(st, 0.8).T rtol = 1e-9
        @test static_state(stf, 0.8).T ≈ static_state(st, 0.8).T rtol = 1e-9
    end

    @testset "zero allocations after warmup" begin
        measured(f::F, args...) where {F} = (f(args...); @allocated f(args...))
        st = GasState(air, 288.15, 101325.0)
        sos(g, T) = speed_of_sound(g, T)
        sosst(x) = speed_of_sound(x)
        mch(x, V) = mach(x, V)
        stag(x, M) = stagnation_state(x, M)
        stat(x, M) = static_state(x, M)
        @test measured(sos, air, 288.15) == 0
        @test measured(sosst, st) == 0
        @test measured(mch, st, 100.0) == 0
        @test measured(stag, st, 0.8) == 0
        @test measured(stat, st, 0.8) == 0
    end

    @testset "ForwardDiff: analytic tangents through the flow verbs" begin
        D = ForwardDiff.derivative
        st = GasState(air, 288.15, 101325.0)
        # d a/dT for √(γRT): finite-difference check (γ varies with T)
        T0 = 600.0
        @test D(T -> speed_of_sound(air, T), T0) ≈
              (speed_of_sound(air, T0 + 1e-3) - speed_of_sound(air, T0 - 1e-3)) / 2e-3 rtol = 1e-6
        # ∂Tt/∂M of stagnation at fixed static state:
        #   Tt = T_of_h(h(Ts) + ½M²a²),  a² = γ(Ts)·R·Ts
        #   ⇒ dTt/dM = M·a² / cp(Tt)   (IFT for T_of_h)
        M = 0.8
        a2 = IdealGasThermo.gamma(air, st.T) * IdealGasThermo.R(air) * st.T
        Tt = stagnation_state(st, M).T
        @test D(m -> stagnation_state(st, m).T, M) ≈
              M * a2 / IdealGasThermo.cp(air, Tt) rtol = 1e-9
        # the stagnation pressure rail differentiates too (vs central FD)
        @test D(m -> stagnation_state(st, m).P, M) ≈
              (stagnation_state(st, M + 1e-5).P - stagnation_state(st, M - 1e-5).P) / 2e-5 rtol = 1e-6
    end

end
