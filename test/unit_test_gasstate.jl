using ForwardDiff

@testset "GasState and process verbs" begin

    # Gas1D is the deprecated legacy reference (ADR-0002); silence its
    # depwarn where it is constructed for agreement checks
    quiet_Gas1D() =
        Base.CoreLogging.with_logger(() -> Gas1D(), Base.CoreLogging.NullLogger())

    @testset "GasState construction and accessors" begin
        air = FrozenGas(DryAir)
        st = GasState(air, 288.15, 101325.0)
        # an isbits value: stack-allocated, immutable, thread-safe
        @test isbitstype(typeof(st))
        @test st.gas === air
        @test st.T === 288.15
        @test st.P === 101325.0
        # immutable: states are replaced, never mutated
        @test_throws ErrorException st.T = 300.0
        # mixed Real arguments promote (the F type parameter is shared)
        @test GasState(air, 288, 101325.0) === GasState(air, 288.0, 101325.0)
        # accessors forward to the gas property functions at st.T
        for (T, P) in [(288.15, 101325.0), (1600.0, 12.0e5), (999.9, 5.0e4)]
            s = GasState(air, T, P)
            @test IdealGasThermo.cp(s) === IdealGasThermo.cp(air, T)
            @test IdealGasThermo.h(s) === IdealGasThermo.h(air, T)
            @test IdealGasThermo.s0(s) === IdealGasThermo.s0(air, T)
            @test IdealGasThermo.gamma(s) === IdealGasThermo.gamma(air, T)
            @test IdealGasThermo.R(s) === IdealGasThermo.R(air)
        end
        # works for FastFrozenGas too (not isbits — holds table Vectors)
        fg = FastFrozenGas(air)
        stf = GasState(fg, 600.0, 2.0e5)
        @test IdealGasThermo.h(stf) === IdealGasThermo.h(air, 600.0)
    end

    @testset "entropy and density agree with legacy Gas1D" begin
        air = FrozenGas(DryAir)
        g1d = quiet_Gas1D() # legacy stateful layer, same DryAir composition
        for (T, P) in [(288.15, 101325.0), (626.0, 12.159e5), (1600.0, 12.159e5),
                       (882.6, 101325.0), (250.0, 2.0e4)]
            set_TP!(g1d, T, P)
            st = GasState(air, T, P)
            @test entropy(st) ≈ g1d.s rtol = 1e-10   # s = s0(T) − R·ln(P/Pstd)
            @test density(st) ≈ g1d.ρ rtol = 1e-10   # ρ = P/(R·T) [kg/m³]
        end
        # unexported short aliases for qualified use
        st = GasState(air, 288.15, 101325.0)
        @test IdealGasThermo.s(st) === entropy(st)
        @test IdealGasThermo.rho(st) === density(st)
        # entropy at standard pressure is the entropy complement itself
        @test entropy(GasState(air, 500.0, IdealGasThermo.Pstd)) ≈
              IdealGasThermo.s0(air, 500.0) rtol = 1e-14
    end

    @testset "scalar compress/expand vs legacy turbo verbs" begin
        air = FrozenGas(DryAir)
        # legacy compress(gas::AbstractGas, PR, ηp): tolerance-limited Newton
        for (T1, PR, etap) in [(288.15, 12.0, 1.0), (288.15, 12.0, 0.9),
                               (500.0, 30.0, 0.85), (900.0, 1.05, 1.0)]
            g1d = quiet_Gas1D()
            set_TP!(g1d, T1, 101325.0)
            IdealGasThermo.compress(g1d, PR, etap)
            @test compress(air, T1, PR; ηp = etap) ≈ g1d.T rtol = 1e-7
        end
        # legacy expand(gas::AbstractGas, PR, ηp) takes PR ≤ 1; the new verb
        # takes the SAME physical process as ratio ≥ 1 — direction is in the verb
        for (T1, PR, etap) in [(1600.0, 4.0, 1.0), (1600.0, 4.0, 0.9),
                               (1323.8, 5.3, 0.92), (900.0, 1.05, 1.0)]
            g1d = quiet_Gas1D()
            set_TP!(g1d, T1, 101325.0)
            IdealGasThermo.expand(g1d, 1 / PR, etap)
            @test expand(air, T1, PR; ηp = etap) ≈ g1d.T rtol = 1e-7
        end
        # both verbs demand ratio ≥ 1: direction lives in the verb, not the number
        @test_throws ArgumentError compress(air, 288.15, 0.5)
        @test_throws ArgumentError expand(air, 1600.0, 0.25)
        # round trip at ηp = 1: expand undoes compress exactly
        for (T1, PR) in [(288.15, 12.0), (500.0, 30.0)]
            @test expand(air, compress(air, T1, PR), PR) ≈ T1 rtol = 1e-10
        end
        # ηp = 1 reduces both to the pure isentrope
        @test compress(air, 288.15, 12.0) ≈
              IdealGasThermo.T_isentropic(air, 288.15, 12.0) rtol = 1e-14
        @test expand(air, 1600.0, 4.0) ≈
              IdealGasThermo.T_isentropic(air, 1600.0, 0.25) rtol = 1e-14
        # the verbs work for every gas flavor through the same engine
        fg = FastFrozenGas(air)
        @test compress(fg, 288.15, 12.0) ≈ compress(air, 288.15, 12.0) rtol = 1e-10
        @test expand(fg, 1600.0, 4.0; ηp = 0.9) ≈
              expand(air, 1600.0, 4.0; ηp = 0.9) rtol = 1e-10
    end

    @testset "state compress/expand/expand_to" begin
        air = FrozenGas(DryAir)
        st2 = GasState(air, 288.15, 101325.0)
        # compress: T from the scalar verb, P multiplied by PR
        st3 = compress(st2, 12.0; ηp = 0.9)
        @test st3 isa GasState
        @test st3.T ≈ compress(air, 288.15, 12.0; ηp = 0.9) rtol = 1e-14
        @test st3.P ≈ 12.0 * 101325.0 rtol = 1e-14
        @test st3.gas === air
        # expand: T from the scalar verb, P divided by PR
        st5 = GasState(air, 1600.0, 12.0 * 101325.0)
        st9 = expand(st5, 4.0; ηp = 0.92)
        @test st9.T ≈ expand(air, 1600.0, 4.0; ηp = 0.92) rtol = 1e-14
        @test st9.P ≈ 3.0 * 101325.0 rtol = 1e-14
        # ratio validation carries through the state layer
        @test_throws ArgumentError compress(st2, 0.5)
        @test_throws ArgumentError expand(st5, 0.25)
        # expand_to: the nozzle convenience — hits the target pressure exactly
        stn = expand_to(st5, 101325.0)
        @test stn.P == 101325.0
        @test stn.T ≈ expand(air, 1600.0, 12.0) rtol = 1e-14
        stn9 = expand_to(st5, 101325.0; ηp = 0.92)
        @test stn9.T ≈ expand(air, 1600.0, 12.0; ηp = 0.92) rtol = 1e-14
        # expand_to cannot raise the pressure
        @test_throws ArgumentError expand_to(st2, 2.0 * 101325.0)
        # round trip: expand undoes compress, state-for-state, at ηp = 1
        rt = expand(compress(st2, 12.0), 12.0)
        @test rt.T ≈ st2.T rtol = 1e-10
        @test rt.P ≈ st2.P rtol = 1e-12
    end

    @testset "isentropic-efficiency ηs form of compress/expand" begin
        air = FrozenGas(DryAir)
        h = IdealGasThermo.h
        # ηs = 1 is the loss-free isentrope, identical to the default and ηp = 1
        @test compress(air, 288.15, 12.0; ηs = 1.0) ≈ compress(air, 288.15, 12.0) rtol = 1e-12
        @test expand(air, 1600.0, 4.0; ηs = 1.0) ≈ expand(air, 1600.0, 4.0) rtol = 1e-12

        # compressor: a real machine (ηs < 1) overshoots the ideal outlet T;
        # the enthalpy rise is the ideal rise divided by ηs
        T1, PR, ηs = 288.15, 12.0, 0.85
        T2s = compress(air, T1, PR)                    # loss-free outlet
        T2 = compress(air, T1, PR; ηs = ηs)
        @test T2 > T2s
        @test h(air, T2) - h(air, T1) ≈ (h(air, T2s) - h(air, T1)) / ηs rtol = 1e-10

        # turbine: a real machine (ηs < 1) leaves the gas hotter than ideal;
        # the enthalpy drop is ηs times the ideal drop
        T1t, PRt, ηst = 1600.0, 4.0, 0.9
        T2st = expand(air, T1t, PRt)                   # loss-free outlet
        T2t = expand(air, T1t, PRt; ηs = ηst)
        @test T2t > T2st
        @test h(air, T1t) - h(air, T2t) ≈ ηst * (h(air, T1t) - h(air, T2st)) rtol = 1e-10

        # the two conventions are consistent: take a polytropic outlet, measure
        # its isentropic efficiency, feed it back as ηs → the same outlet T
        for (T0, P0, ηp) in ((288.15, 12.0, 0.9), (500.0, 30.0, 0.85))
            T2p = compress(air, T0, P0; ηp = ηp)
            ηs_eff = (h(air, compress(air, T0, P0)) - h(air, T0)) / (h(air, T2p) - h(air, T0))
            @test compress(air, T0, P0; ηs = ηs_eff) ≈ T2p rtol = 1e-10
        end
        for (T0, P0, ηp) in ((1600.0, 4.0, 0.92), (1323.8, 5.3, 0.88))
            T2p = expand(air, T0, P0; ηp = ηp)
            ηs_eff = (h(air, T0) - h(air, T2p)) / (h(air, T0) - h(air, expand(air, T0, P0)))
            @test expand(air, T0, P0; ηs = ηs_eff) ≈ T2p rtol = 1e-10
        end

        # at most one efficiency convention — both is an error, scalar and state
        st = GasState(air, 288.15, 101325.0)
        sthot = GasState(air, 1600.0, 12.159e5)
        @test_throws ArgumentError compress(air, 288.15, 12.0; ηp = 0.9, ηs = 0.9)
        @test_throws ArgumentError expand(air, 1600.0, 4.0; ηp = 0.9, ηs = 0.9)
        @test_throws ArgumentError compress(st, 12.0; ηp = 0.9, ηs = 0.9)
        @test_throws ArgumentError expand(sthot, 4.0; ηp = 0.9, ηs = 0.9)
        @test_throws ArgumentError expand_to(sthot, 101325.0; ηp = 0.9, ηs = 0.9)

        # state layer: ηs degrades T but lands on the SAME pressure as ηp
        c_ηp = compress(st, 12.0; ηp = 0.85)
        c_ηs = compress(st, 12.0; ηs = 0.85)
        @test c_ηp.P ≈ c_ηs.P rtol = 1e-14           # outlet pressure is P·PR either way
        @test c_ηs.T ≈ compress(air, st.T, 12.0; ηs = 0.85) rtol = 1e-14
        e_ηs = expand(sthot, 4.0; ηs = 0.9)
        @test e_ηs.P ≈ sthot.P / 4.0 rtol = 1e-14
        @test e_ηs.T ≈ expand(air, sthot.T, 4.0; ηs = 0.9) rtol = 1e-14
        # expand_to with ηs hits the named pressure exactly
        en = expand_to(sthot, 101325.0; ηs = 0.9)
        @test en.P == 101325.0
        @test en.T ≈ expand(air, sthot.T, sthot.P / 101325.0; ηs = 0.9) rtol = 1e-14

        # works through the accelerated flavor
        fg = FastFrozenGas(air)
        @test compress(fg, 288.15, 12.0; ηs = 0.85) ≈ compress(air, 288.15, 12.0; ηs = 0.85) rtol = 1e-10

        # zero allocations on the ηs path
        measured(f::F, args...) where {F} = (f(args...); @allocated f(args...))
        cmp_ηs(x, PR, e) = compress(x, PR; ηs = e)
        expd_ηs(x, PR, e) = expand(x, PR; ηs = e)
        @test measured(cmp_ηs, st, 12.0, 0.85) == 0
        @test measured(expd_ηs, sthot, 4.0, 0.9) == 0

        # ForwardDiff flows through the ηs path (built on the IFT-ruled engines)
        T2 = compress(air, 288.15, 12.0; ηs = 0.85)
        @test ForwardDiff.derivative(pr -> compress(air, 288.15, pr; ηs = 0.85), 12.0) ≈
              (compress(air, 288.15, 12.0 + 1e-5; ηs = 0.85) -
               compress(air, 288.15, 12.0 - 1e-5; ηs = 0.85)) / 2e-5 rtol = 1e-6
    end

    @testset "add_heat: constant-P enthalpy change, signed q" begin
        air = FrozenGas(DryAir)
        st = GasState(air, 626.0, 12.159e5)
        q = 8.0e5 # J/kg of heating
        sth = add_heat(st, q)
        @test sth.P === st.P                                # constant pressure
        @test IdealGasThermo.h(sth) ≈ IdealGasThermo.h(st) + q rtol = 1e-12
        @test sth.T > st.T
        # negative q cools
        stc = add_heat(st, -3.0e5)
        @test stc.P === st.P
        @test IdealGasThermo.h(stc) ≈ IdealGasThermo.h(st) - 3.0e5 rtol = 1e-12
        @test stc.T < st.T
        # q = 0 is the identity to the inversion tolerance
        @test add_heat(st, 0.0).T ≈ st.T rtol = 1e-12
    end

    @testset "add_work/extract_work vs legacy set_Δh!" begin
        air = FrozenGas(DryAir)
        w = 3.5e5 # J/kg, always nonnegative — direction lives in the verb
        for etap in (1.0, 0.9)
            # add_work ↔ legacy set_Δh!(gas, +w, ηp): P2 = P1·exp(ηp/R·Δs0)
            st1 = GasState(air, 288.15, 101325.0)
            g1d = quiet_Gas1D()
            set_TP!(g1d, st1.T, st1.P)
            set_Δh!(g1d, w, etap)
            st2 = add_work(st1, w; ηp = etap)
            @test st2.T ≈ g1d.T rtol = 1e-9
            @test st2.P ≈ g1d.P rtol = 1e-8
            @test IdealGasThermo.h(st2) ≈ IdealGasThermo.h(st1) + w rtol = 1e-12
            # extract_work ↔ legacy set_Δh!(gas, −w, 1/ηp): P2 = P1·exp(Δs0/(ηp·R))
            st5 = GasState(air, 1600.0, 12.159e5)
            g1d = quiet_Gas1D()
            set_TP!(g1d, st5.T, st5.P)
            set_Δh!(g1d, -w, 1 / etap)
            st6 = extract_work(st5, w; ηp = etap)
            @test st6.T ≈ g1d.T rtol = 1e-9
            @test st6.P ≈ g1d.P rtol = 1e-8
            @test IdealGasThermo.h(st6) ≈ IdealGasThermo.h(st5) - w rtol = 1e-12
        end
        # at ηp = 1 both reduce to the isentrope: P follows pressure_ratio
        st1 = GasState(air, 288.15, 101325.0)
        st2 = add_work(st1, w)
        @test st2.P ≈ st1.P * IdealGasThermo.pressure_ratio(air, st1.T, st2.T) rtol = 1e-12
        st6 = extract_work(st2, w)
        @test st6.T ≈ st1.T rtol = 1e-10
        @test st6.P ≈ st1.P rtol = 1e-9
        # w must be nonnegative — the sign convention is the verb's job
        @test_throws ArgumentError add_work(st1, -1.0)
        @test_throws ArgumentError extract_work(st1, -1.0)
    end

    @testset "temperature facade: h-inversion only, isentrope form removed" begin
        air = FrozenGas(DryAir)
        # the h form is unchanged
        for T in [250.0, 500.0, 1600.0]
            @test temperature(air, h = IdealGasThermo.h(air, T)) ≈ T rtol = 1e-10
        end
        # the isentrope form is a process, not an inversion — it now lives in
        # the compress/expand verbs; the old kwargs throw with a pointer there
        @test_throws ArgumentError temperature(air, T1 = 288.15, PR = 12.0)
        @test_throws ArgumentError temperature(air, T1 = 288.15, PR = 12.0, ηp = 0.9)
        @test_throws ArgumentError temperature(air, T1 = 288.15)
        @test_throws ArgumentError temperature(air, PR = 12.0)
        @test_throws ArgumentError temperature(air, h = 5.0e5, T1 = 288.15, PR = 12.0)
        @test_throws ArgumentError temperature(air)
        # same facade behavior for the accelerated flavors
        fg = FastFrozenGas(air)
        @test temperature(fg, h = IdealGasThermo.h(air, 700.0)) ≈ 700.0 rtol = 1e-10
        @test_throws ArgumentError temperature(fg, T1 = 288.15, PR = 12.0)
    end

    @testset "zero allocations after warmup" begin
        air = FrozenGas(DryAir)
        st = GasState(air, 288.15, 101325.0)
        sthot = GasState(air, 1600.0, 12.159e5)
        # single-method local helpers (see unit_test_fastgas.jl): closures
        # with kwarg defaults have non-specializing kwsorters, and
        # multi-method local functions get boxed when captured — both would
        # charge allocations to the helper, not the package verbs
        measured(f::F, args...) where {F} = (f(args...); @allocated f(args...))
        mk(g, T, P) = GasState(g, T, P)
        ent(x) = entropy(x)
        den(x) = density(x)
        cpst(x) = IdealGasThermo.cp(x)
        hst(x) = IdealGasThermo.h(x)
        s0st(x) = IdealGasThermo.s0(x)
        gammast(x) = IdealGasThermo.gamma(x)
        Rst(x) = IdealGasThermo.R(x)
        cmp(x, PR) = compress(x, PR)
        cmp_eta(x, PR, e) = compress(x, PR; ηp = e)
        expd(x, PR) = expand(x, PR)
        expd_eta(x, PR, e) = expand(x, PR; ηp = e)
        expto(x, P2) = expand_to(x, P2)
        expto_eta(x, P2, e) = expand_to(x, P2; ηp = e)
        addq(x, q) = add_heat(x, q)
        addw(x, w) = add_work(x, w)
        addw_eta(x, w, e) = add_work(x, w; ηp = e)
        extw(x, w) = extract_work(x, w)
        extw_eta(x, w, e) = extract_work(x, w; ηp = e)
        @test measured(mk, air, 288.15, 101325.0) == 0
        @test measured(ent, st) == 0
        @test measured(den, st) == 0
        @test measured(cpst, st) == 0
        @test measured(hst, st) == 0
        @test measured(s0st, st) == 0
        @test measured(gammast, st) == 0
        @test measured(Rst, st) == 0
        @test measured(cmp, st, 12.0) == 0
        @test measured(cmp_eta, st, 12.0, 0.9) == 0
        @test measured(expd, sthot, 4.0) == 0
        @test measured(expd_eta, sthot, 4.0, 0.9) == 0
        @test measured(expto, sthot, 101325.0) == 0
        @test measured(expto_eta, sthot, 101325.0, 0.9) == 0
        @test measured(addq, st, 5.0e5) == 0
        @test measured(addw, st, 3.5e5) == 0
        @test measured(addw_eta, st, 3.5e5, 0.9) == 0
        @test measured(extw, sthot, 3.5e5) == 0
        @test measured(extw_eta, sthot, 3.5e5, 0.9) == 0
        # the scalar kernel verbs too
        cmps(g, T1, PR) = compress(g, T1, PR)
        expds(g, T1, PR) = expand(g, T1, PR)
        @test measured(cmps, air, 288.15, 12.0) == 0
        @test measured(expds, air, 1600.0, 4.0) == 0
    end

    @testset "ForwardDiff through the state verbs" begin
        air = FrozenGas(DryAir)
        st = GasState(air, 288.15, 101325.0)
        D = ForwardDiff.derivative
        # GasState carries Duals through its parametric F
        # ∂(compress(st, PR).T)/∂PR = R·T2/(PR·cp(T2)) — the existing IFT
        # closed form (s0(T2) = s0(T1) + R·ln(PR), ηp = 1)
        PR = 12.0
        T2 = compress(air, st.T, PR)
        @test D(pr -> compress(st, pr).T, PR) ≈
              air.R * T2 / (PR * IdealGasThermo.cp(air, T2)) rtol = 1e-10
        # the pressure rail differentiates trivially: ∂(P·PR)/∂PR = P
        @test D(pr -> compress(st, pr).P, PR) ≈ st.P rtol = 1e-14
        # Dual T and P in the state: derivative w.r.t. the inlet temperature
        @test D(t -> compress(GasState(air, t, 101325.0), PR).T, 288.15) ≈
              IdealGasThermo.cp(air, 288.15) * T2 /
              (IdealGasThermo.cp(air, T2) * 288.15) rtol = 1e-10
        # Duals through q: ∂(add_heat(st, q).T)/∂q = 1/cp(T2) (IFT for T_of_h)
        q = 5.0e5
        Tq = add_heat(st, q).T
        @test D(qq -> add_heat(st, qq).T, q) ≈
              1 / IdealGasThermo.cp(air, Tq) rtol = 1e-10
        # Duals through w in both work verbs (T and P rails)
        w = 3.5e5
        @test D(ww -> add_work(st, ww; ηp = 0.9).T, w) ≈
              1 / IdealGasThermo.cp(air, add_work(st, w; ηp = 0.9).T) rtol = 1e-10
        sthot = GasState(air, 1600.0, 12.159e5)
        @test D(ww -> extract_work(sthot, ww; ηp = 0.9).P, w) ≈
              (extract_work(sthot, w + 0.5; ηp = 0.9).P -
               extract_work(sthot, w - 0.5; ηp = 0.9).P) rtol = 1e-6
        # gradient of a 3-step mini-chain (compress → add_heat → extract_work)
        # w.r.t. [PR, q] against central finite differences (fixed composition)
        chain(x) = extract_work(add_heat(compress(st, x[1]; ηp = 0.9), x[2]),
                                4.0e5; ηp = 0.92)
        chainT(x) = chain(x).T
        chainP(x) = chain(x).P
        x0 = [12.0, 5.0e5]
        for f in (chainT, chainP)
            g = ForwardDiff.gradient(f, x0)
            for i in 1:2
                δ = x0[i] * 1e-6
                xp = copy(x0); xp[i] += δ
                xm = copy(x0); xm[i] -= δ
                @test g[i] ≈ (f(xp) - f(xm)) / (2δ) rtol = 1e-6
            end
        end
        # Dual-typed state evaluation stays allocation-free through the rules
        measured(f::F, args...) where {F} = (f(args...); @allocated f(args...))
        cmpT(x, pr) = compress(x, pr).T
        @test measured(D, pr -> cmpT(st, pr), 12.0) == 0
    end

end
