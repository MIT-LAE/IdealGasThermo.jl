# Gas dynamics on the pure core: speed of sound (in frozengas.jl), Mach
# number, and the `static_state` ↔ `stagnation_state` pair. These supersede
# the legacy constant-γ stagnation helpers (isenTR/isenPR/get_static) and the
# mutable gas_Mach!. See docs/adr/0005.
#
# Stagnation is treated as an *isentropic reference state* (textbook
# definition): total enthalpy is the static enthalpy plus ½V², and the total
# pressure is the loss-free pressure that entropy-conservation implies. The
# relations are built on the enthalpy/entropy curves and the inversion, NOT
# the constant-γ shortcuts `1 + ½(γ−1)M²` and `(…)^(γ/(γ−1))` — for a
# thermally-perfect gas those ratios drift with temperature (ADR-0004's
# temp_ratio lesson). A lossy ram/diffuser is the caller's job: compose with a
# pressure-recovery factor or an `expand`.

"""
    mach(gas, T, V)

Mach number `M = V / a` of a flow of speed `V` [m/s] through `gas` at
temperature `T` [K], with the local speed of sound
`a = `[`speed_of_sound`](@ref)`(gas, T)`. Pure function of `(gas, T, V)`.
"""
@inline mach(gas::PureGas, T, V) = V / speed_of_sound(gas, T)

"""
    mach(st::GasState, V)

Mach number of a flow of speed `V` [m/s] at the state `st`, using the speed
of sound at `st.T`.
"""
@inline mach(st::GasState, V) = V / speed_of_sound(st)

"""
    stagnation_state(st::GasState, M) -> GasState

The stagnation (total) state reached by **isentropically** bringing the flow
at static state `st`, Mach number `M ≥ 0`, to rest. The flow speed is
`V = M·a` with `a = `[`speed_of_sound`](@ref)`(st)`; the total enthalpy is
`h_t = h(st) + ½V²`, giving `T_t` by the enthalpy inversion, and the total
pressure is the loss-free `P_t = P·exp((s0(T_t) − s0(T))/R) ≥ P`.

For a thermally-perfect gas this is the exact enthalpy/entropy result, *not*
the constant-γ relation `1 + ½(γ−1)M²` (ADR-0004/0005). The inverse is
[`static_state`](@ref): `static_state(stagnation_state(st, M), M)` returns
`st`. Pure, zero-allocation; `M` may be a ForwardDiff `Dual`.
"""
function stagnation_state(st::GasState, M)
    M ≥ 0 ||
        throw(ArgumentError("stagnation_state: Mach number must be ≥ 0 (got M = $M)"))
    gas = st.gas
    V = M * speed_of_sound(gas, st.T)
    Tt = T_of_h(gas, h(gas, st.T) + V^2 / 2)
    Pt = st.P * exp((s0(gas, Tt) - s0(gas, st.T)) / R(gas))
    GasState(gas, Tt, Pt)
end

"""
    static_state(st::GasState, M) -> GasState

The static state at Mach number `M ≥ 0` whose isentropic stagnation state is
`st` — the inverse of [`stagnation_state`](@ref). Solves
`h(T) + ½(M·a(T))² = h(st)` for the static temperature `T` (the Mach number
is defined with the *local* sound speed `a(T)`, so `a` is evaluated at the
static `T`), then `P = P_t·exp((s0(T) − s0(T_t))/R) ≤ P_t`.

The residual `h(T) + ½M²·γ(T)·R·T − h_t` is strictly increasing in `T`, so a
bounded Newton solve from the constant-γ estimate converges (relative
tolerance 1e-12, ≤ 30 iterations; errors otherwise). Pure, zero-allocation.
"""
function static_state(st::GasState, M)
    M ≥ 0 || throw(ArgumentError("static_state: Mach number must be ≥ 0 (got M = $M)"))
    gas = st.gas
    Tt, Pt = st.T, st.P
    ht = h(gas, Tt)
    Rg = R(gas)
    halfM2 = M^2 / 2
    # constant-γ estimate seeds the Newton solve
    γt = gamma(gas, Tt)
    T = Tt / (1 + (γt - 1) * halfM2)
    for _ = 1:NEWTON_MAXITER
        γ = gamma(gas, T)
        # f = h(T) + ½M²·a²(T) − h_t with a² = γRT; f′ ≈ cp + ½M²γR (the
        # slowly-varying dγ/dT term is dropped — it only affects the rate,
        # not the converged root, which satisfies f = 0 to tolerance)
        f = h(gas, T) + halfM2 * γ * Rg * T - ht
        df = cp(gas, T) + 2 * halfM2 * γ * Rg
        dT = -f / df
        T += dT
        if abs(dT) ≤ NEWTON_RTOL * abs(T)
            Ps = Pt * exp((s0(gas, T) - s0(gas, Tt)) / Rg)
            return GasState(gas, T, Ps)
        end
    end
    error("static_state: Mach inversion did not converge for M = $M (last T = $T)")
end
