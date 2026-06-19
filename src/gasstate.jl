# GasState: an immutable (gas, T, P) record, plus the process verbs that
# transform one state into the next. See docs/adr/0004.

"""
    GasState(gas, T, P)

Immutable thermodynamic state: a substance (a [`FrozenGas`](@ref) or
[`FastFrozenGas`](@ref)) pinned at temperature `T` [K] and pressure `P`
[Pa]. A `GasState` is a *value* — it is never mutated; process verbs
([`compress`](@ref), [`expand`](@ref), [`add_heat`](@ref),
[`add_work`](@ref), [`extract_work`](@ref), [`expand_to`](@ref)) return a
NEW state. `GasState{FrozenGas{Float64},Float64}` is `isbits`
(stack-allocated, thread-safe by construction).

This is ergonomics, not architecture (ADR-0004): the substance remains a
pure set of property curves; the record only carries the caller's (T, P)
pair so the two travel together through a process chain.

Property access is by read-only accessor functions, not field magic:
`cp(st)`, `h(st)`, `s0(st)`, `gamma(st)`, `R(st)` forward to the gas at
`st.T`; [`entropy`](@ref) and [`density`](@ref) are the two quantities the
(T, P) pair uniquely enables. The numeric type `F` widens automatically
(e.g. to `ForwardDiff.Dual`), and mixed `Real` arguments to the
constructor promote.

```julia-repl
julia> st = GasState(FrozenGas(DryAir), 288.15, 101325.0);

julia> density(st)
1.2250291442077996
```
"""
struct GasState{G,F<:Real}
    gas::G
    T::F
    P::F
end

# promoting outer constructor, so GasState(gas, 1600.0, P::Dual) works
GasState(gas, T::Real, P::Real) = GasState(gas, promote(T, P)...)

# accessors: the gas property functions at st.T, arity-1
"""
    cp(st::GasState)

Specific heat at constant pressure [J/kg/K] of the state's gas at `st.T`.
"""
@inline cp(st::GasState) = cp(st.gas, st.T)

"""
    h(st::GasState)

Specific enthalpy [J/kg] of the state's gas at `st.T` (formation-inclusive
datum, as [`h(gas, T)`](@ref h)).
"""
@inline h(st::GasState) = h(st.gas, st.T)

"""
    s0(st::GasState)

Entropy complement φ(T) = ∫cp/T dT [J/kg/K] of the state's gas at `st.T`
(the temperature-only part of entropy; see [`entropy`](@ref) for the full
pressure-dependent entropy).
"""
@inline s0(st::GasState) = s0(st.gas, st.T)

"""
    gamma(st::GasState)

Ratio of specific heats cp/(cp − R) of the state's gas at `st.T`.
"""
@inline gamma(st::GasState) = gamma(st.gas, st.T)

"""
    R(st::GasState)

Specific gas constant [J/kg/K] of the state's gas.
"""
@inline R(st::GasState) = R(st.gas)

"""
    speed_of_sound(st::GasState)

Speed of sound [m/s] of the state's gas at `st.T` — a property of the gas
and temperature only (see [`speed_of_sound(gas, T)`](@ref speed_of_sound)),
needing no pressure.
"""
@inline speed_of_sound(st::GasState) = speed_of_sound(st.gas, st.T)

"""
    entropy(st::GasState)

Full specific entropy `s = s0(T) − R·ln(P/Pstd)` [J/kg/K] — the entropy
complement corrected to the state's pressure (`Pstd` = 101325 Pa).
Also available as the unexported alias `IdealGasThermo.s(st)`.
"""
@inline entropy(st::GasState) = s0(st.gas, st.T) - R(st.gas) * log(st.P / Pstd)

"""
    density(st::GasState)

Mass density `ρ = P/(R·T)` [kg/m³] from the ideal-gas law (`P` [Pa], `R`
[J/kg/K], `T` [K]). Also available as the unexported alias
`IdealGasThermo.rho(st)`.
"""
@inline density(st::GasState) = st.P / (R(st.gas) * st.T)

# unexported short aliases (`s` and `rho` are too collision-prone to export)
@inline s(st::GasState) = entropy(st)
@inline rho(st::GasState) = density(st)

# ---------------------------------------------------------------------------
# Process verbs, scalar kernel layer: T1 -> T2 for any gas flavor.
# Both verbs take pressure ratios ≥ 1 — the direction of the process lives
# in the verb, never in the number (ADR-0004).
# ---------------------------------------------------------------------------

# the gas flavors the process verbs operate on (the pure core); the legacy
# compress/expand methods on the mutable AbstractGas layer live in turbo.jl
const PureGas = Union{FrozenGas,FastFrozenGas}

# Efficiency resolution shared by compress/expand. The ratio verbs accept
# *either* a polytropic efficiency `ηp` (entropy distributed along the path)
# *or* an isentropic efficiency `ηs` (ideal enthalpy change degraded at the
# same pressure ratio) — never both (ADR-0005). `nothing` defaults stand in
# for "unspecified" so the mutual exclusion is detectable; an unspecified
# pair means the loss-free isentrope (ηp = 1). These compile away per
# specialization (the `=== nothing` tests fold to constants), so the verbs
# stay zero-allocation.
@inline _resolve_ηp(ηp) = ηp === nothing ? 1.0 : ηp
@inline function _check_one_efficiency(ηp, ηs, verb)
    ηp === nothing ||
        ηs === nothing ||
        throw(
            ArgumentError(
                "$verb: specify at most one efficiency — ηp (polytropic) or " *
                "ηs (isentropic), not both",
            ),
        )
end

"""
    compress(gas, T1, PR; ηp = 1.0) -> T2
    compress(gas, T1, PR; ηs)       -> T2

Temperature [K] after compression of a [`FrozenGas`](@ref) or
[`FastFrozenGas`](@ref) from `T1` [K] by pressure ratio `PR ≥ 1`
(`ArgumentError` otherwise — to lower the pressure, use [`expand`](@ref):
the direction lives in the verb, never in the number).

The loss model is chosen by *which* efficiency keyword is given (at most one;
`ArgumentError` if both — ADR-0005):

- **polytropic `ηp`** (default `ηp = 1`, the isentrope): solves
  `s0(T2) = s0(T1) + R·ln(PR)/ηp`.
- **isentropic `ηs`**: the ideal enthalpy rise to the same `PR` is degraded
  by `1/ηs`, `h(T2) = h(T1) + (h(T2s) − h(T1))/ηs` with `T2s` the loss-free
  outlet. Both conventions reach the same outlet pressure `P1·PR` (only the
  temperature differs); `ηs = 1` reproduces the isentrope.

Pure, zero-allocation; arguments may be ForwardDiff `Dual`s
(implicit-function-theorem rules from the package extension).
"""
function compress(gas::PureGas, T1, PR; ηp = nothing, ηs = nothing)
    PR ≥ 1 || throw(
        ArgumentError(
            "compress: pressure ratio must be ≥ 1 (got PR = $PR); " *
            "both process verbs take ratios ≥ 1 — use expand to lower the pressure",
        ),
    )
    _check_one_efficiency(ηp, ηs, "compress")
    if ηs === nothing
        _T_polytropic(gas, T1, PR; ηp = _resolve_ηp(ηp))
    else
        h1 = h(gas, T1)
        T2s = _T_polytropic(gas, T1, PR) # loss-free outlet at this PR
        T_from_h(gas, h1 + (h(gas, T2s) - h1) / ηs)
    end
end

"""
    expand(gas, T1, PR; ηp = 1.0) -> T2
    expand(gas, T1, PR; ηs)       -> T2

Temperature [K] after expansion of a [`FrozenGas`](@ref) or
[`FastFrozenGas`](@ref) from `T1` [K] by pressure ratio `PR ≥ 1` — the
ratio of inlet to outlet pressure, so `PR = 4` quarters the pressure
(`ArgumentError` if `PR < 1`; to raise the pressure use [`compress`](@ref)).

The loss model is chosen by *which* efficiency keyword is given (at most one;
`ArgumentError` if both — ADR-0005):

- **polytropic `ηp`** (default `ηp = 1`, the isentrope): solves
  `s0(T2) = s0(T1) + R·ηp·ln(1/PR)` — the expansion convention, matching the
  legacy `expand(gas::AbstractGas, 1/PR, ηp)`.
- **isentropic `ηs`**: the actual enthalpy drop is `ηs` times the loss-free
  drop to the same `PR`, `h(T2) = h(T1) − ηs·(h(T1) − h(T2s))` with `T2s` the
  loss-free outlet. Both conventions reach the same outlet pressure `P1/PR`
  (only the temperature differs); `ηs = 1` reproduces the isentrope.

Pure, zero-allocation; arguments may be ForwardDiff `Dual`s.
"""
function expand(gas::PureGas, T1, PR; ηp = nothing, ηs = nothing)
    PR ≥ 1 || throw(
        ArgumentError(
            "expand: pressure ratio must be ≥ 1 (got PR = $PR); " *
            "both process verbs take ratios ≥ 1 — use compress to raise the pressure",
        ),
    )
    _check_one_efficiency(ηp, ηs, "expand")
    if ηs === nothing
        _T_polytropic(gas, T1, inv(PR); ηp = inv(_resolve_ηp(ηp)))
    else
        h1 = h(gas, T1)
        T2s = _T_polytropic(gas, T1, inv(PR)) # loss-free outlet at this PR
        T_from_h(gas, h1 - ηs * (h1 - h(gas, T2s)))
    end
end

# ---------------------------------------------------------------------------
# Process verbs, state layer: GasState -> NEW GasState (replace, never mutate)
# ---------------------------------------------------------------------------

"""
    compress(st::GasState, PR; ηp = 1.0) -> GasState
    compress(st::GasState, PR; ηs)       -> GasState

Compression of a state by pressure ratio `PR ≥ 1`: a new state at
`(T2, P·PR)` with `T2` from the scalar [`compress`](@ref) verb. The loss
model follows the efficiency keyword — polytropic `ηp` (default, the
isentrope) or isentropic `ηs`, at most one (ADR-0005). The outlet pressure
`P·PR` is the same either way. The input state is untouched.
"""
compress(st::GasState, PR; ηp = nothing, ηs = nothing) =
    GasState(st.gas, compress(st.gas, st.T, PR; ηp = ηp, ηs = ηs), st.P * PR)

"""
    expand(st::GasState, PR; ηp = 1.0) -> GasState
    expand(st::GasState, PR; ηs)       -> GasState

Expansion of a state by pressure ratio `PR ≥ 1` (inlet over outlet
pressure): a new state at `(T2, P/PR)` with `T2` from the scalar
[`expand`](@ref) verb. The loss model follows the efficiency keyword —
polytropic `ηp` (default, the isentrope) or isentropic `ηs`, at most one
(ADR-0005). The input state is untouched. See [`expand_to`](@ref) to name
the outlet pressure instead of the ratio.
"""
expand(st::GasState, PR; ηp = nothing, ηs = nothing) =
    GasState(st.gas, expand(st.gas, st.T, PR; ηp = ηp, ηs = ηs), st.P / PR)

"""
    expand_to(st::GasState, P2; ηp = 1.0) -> GasState
    expand_to(st::GasState, P2; ηs)       -> GasState

Expansion of a state to outlet pressure `P2 ≤ st.P` (`ArgumentError`
otherwise) — the nozzle convenience, equivalent to `expand(st, st.P / P2;
…)`. The loss model follows the efficiency keyword (polytropic `ηp` default,
or isentropic `ηs`, at most one). The returned state's pressure is exactly
the `P2` requested.
"""
function expand_to(st::GasState, P2; ηp = nothing, ηs = nothing)
    P2 ≤ st.P || throw(
        ArgumentError(
            "expand_to: target pressure P2 = $P2 exceeds the state pressure " *
            "$(st.P); use compress to raise the pressure",
        ),
    )
    GasState(st.gas, expand(st.gas, st.T, st.P / P2; ηp = ηp, ηs = ηs), P2)
end

"""
    add_heat(st::GasState, q) -> GasState

Constant-pressure heat addition of `q` `J/kg` (signed: negative `q`
cools): a new state at the same pressure with `h(out) = h(st) + q`, the
temperature from the enthalpy inversion. The input state is untouched.
"""
add_heat(st::GasState, q) = GasState(st.gas, T_from_h(st.gas, h(st.gas, st.T) + q), st.P)

"""
    add_work(st::GasState, w; ηp = 1.0) -> GasState

Polytropic work *addition* of `w ≥ 0` `J/kg` (`ArgumentError` if `w < 0` —
the sign convention lives in the verb; see [`extract_work`](@ref)): a new
state with `h(out) = h(st) + w` (temperature from the enthalpy inversion)
and the pressure walked along the polytrope,
`P2 = P1·exp(ηp/R·(s0(T2) − s0(T1)))` — the compression ηp convention,
matching the legacy `set_Δh!(gas, +w, ηp)`. At `ηp = 1` the pressure
follows the exact isentrope ([`pressure_ratio`](@ref)).
"""
add_work(st::GasState, w; ηp = 1.0) = _work(st, w, w, ηp)

"""
    extract_work(st::GasState, w; ηp = 1.0) -> GasState

Polytropic work *extraction* of `w ≥ 0` `J/kg` (`ArgumentError` if
`w < 0`): a new state with `h(out) = h(st) − w` (temperature from the
enthalpy inversion) and the pressure walked along the polytrope,
`P2 = P1·exp((s0(T2) − s0(T1))/(ηp·R))` — the expansion ηp convention,
matching the legacy `set_Δh!(gas, −w, 1/ηp)` (the legacy verb applies
`exp(ηp/R·Δs0)` with a caller-chosen ηp convention; here the convention is
the verb's job). The turbine work-balance verb.
"""
extract_work(st::GasState, w; ηp = 1.0) = _work(st, w, -w, ηp, true)

# shared kernel: enthalpy change Δh = ±w, pressure along the polytrope
# P2 = P1·exp(K/R·Δs0) with K = ηp (work added) or 1/ηp (work extracted)
function _work(st::GasState, w, Δh, ηp, extracting::Bool = false)
    w ≥ 0 || throw(
        ArgumentError(
            (extracting ? "extract_work" : "add_work") *
            ": w must be ≥ 0 (got w = $w); the sign of the enthalpy change " *
            "lives in the verb — use the opposite verb to reverse it",
        ),
    )
    gas = st.gas
    K = extracting ? inv(ηp) : ηp
    T2 = T_from_h(gas, h(gas, st.T) + Δh)
    P2 = st.P * exp(K / R(gas) * (s0(gas, T2) - s0(gas, st.T)))
    GasState(gas, T2, P2)
end

Base.show(io::IO, st::GasState) = print(
    io,
    "GasState(T = ",
    round(Float64(st.T), digits = 3),
    " K, P = ",
    round(Float64(st.P), digits = 1),
    " Pa)",
)
