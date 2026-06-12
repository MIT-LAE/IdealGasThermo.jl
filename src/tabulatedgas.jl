"""
    TabulatedGas

A [`FrozenGas`](@ref) bundled with two precomputed cubic-Hermite *inverse*
tables on uniform grids — h → T and s0 → T — that accelerate the inversions
[`T_of_h`](@ref) and [`T_isentropic`](@ref). Construct with
[`tabulate`](@ref).

Only the inverses are accelerated. The forward property functions
(`cp`, `h`, `s0`, `props`, `gamma`, `R`, `pressure_ratio`) forward to the
wrapped `FrozenGas` unchanged: the forward functions never change; both
tables sample the same NASA-9 polynomials.

Two tiers:

- **Default tier (exact)** — `T_of_h(tg, hspec)` and
  `T_isentropic(tg, T1, PR; ηp)` use the table only as a *seed* for the same
  bounded Newton iteration and convergence criterion as `FrozenGas`
  (relative tolerance 1e-12 on the temperature step, ≤ 30 iterations).
  The answer satisfies the package's documented inversion contract; the seed
  just cuts the iteration count to typically 1–2 polish steps. If the target
  falls outside the tabulated range the table is bypassed entirely and the
  plain `FrozenGas` Newton solve (cold start) is used — never extrapolation.
- **Interp tier (opt-in, approximate)** — [`T_of_h_interp`](@ref) and
  [`T_isentropic_interp`](@ref) evaluate the Hermite table only, no polish.
  Accuracy |ΔT/T| ≲ 2e-9 over the table range for the default N = 256
  (measured for dry air over T ∈ [250, 2200] K: 5.8e-10 for the h table,
  1.4e-9 for the s0 table). Out of range they throw a `DomainError` — the
  approximate tier has no silent fallback.
"""
struct TabulatedGas
    gas::FrozenGas{Float64}
    # h → T table: T at uniform h nodes, node slopes dT/dh = 1/cp
    hmin::Float64
    hmax::Float64
    invdh::Float64        # 1/spacing for O(1) cell lookup
    dh::Float64
    Th::Vector{Float64}
    Mh::Vector{Float64}
    # s0 → T table: T at uniform s0 nodes, node slopes dT/ds0 = T/cp
    s0min::Float64
    s0max::Float64
    invds0::Float64
    ds0::Float64
    Ts::Vector{Float64}
    Ms::Vector{Float64}
end

# Forward property functions forward to the wrapped gas exactly — the
# forward functions never change; both tables sample the same NASA-9
# polynomials. Only the inversions are accelerated.
"""
    R(tg::TabulatedGas)

Specific gas constant [J/kg/K] of the wrapped [`FrozenGas`](@ref).
"""
@inline R(tg::TabulatedGas) = R(tg.gas)
@inline cp(tg::TabulatedGas, T) = cp(tg.gas, T)
@inline h(tg::TabulatedGas, T) = h(tg.gas, T)
@inline s0(tg::TabulatedGas, T) = s0(tg.gas, T)
@inline gamma(tg::TabulatedGas, T) = gamma(tg.gas, T)
@inline props(tg::TabulatedGas, T) = props(tg.gas, T)
@inline pressure_ratio(tg::TabulatedGas, T1, T2) = pressure_ratio(tg.gas, T1, T2)

# Exact Newton inverse of s0 (table construction only; same tolerance as
# the T_of_h / T_isentropic Newton solves). s0 is strictly monotonic in T,
# so for targets inside [s0(Tlo), s0(Thi)] the root lies in [Tlo, Thi];
# clamping each iterate to that bracket keeps the log(T) evaluations valid.
function _T_of_s0(gas::FrozenGas, starget, Tlo, Thi)
    T = 0.5 * (Tlo + Thi)
    for _ = 1:NEWTON_MAXITER
        Tnew = clamp(T + (starget - s0(gas, T)) * T / cp(gas, T), Tlo, Thi) # ds0/dT = cp/T
        dT = Tnew - T
        T = Tnew
        if abs(dT) ≤ NEWTON_RTOL * abs(T)
            return T
        end
    end
    error("_T_of_s0 did not converge for starget = $starget (last T = $T)")
end

"""
    tabulate(gas::FrozenGas; N::Int=256, Tmin=200.0, Tmax=2400.0)

Build a [`TabulatedGas`](@ref): precompute cubic-Hermite inverse tables
h → T and s0 → T on uniform N-node grids spanning `T ∈ [Tmin, Tmax]` [K].
Node temperatures come from the exact `FrozenGas` Newton inversions; node
slopes are the exact analytic inverses dT/dh = 1/cp and dT/ds0 = T/cp, so
the interpolant is C¹ and matches the polynomials to |ΔT/T| ≲ 1e-9 at the
default N = 256. Construction allocates (~4 KB per table); all subsequent
inversion calls are pure and zero-allocation.
"""
function tabulate(gas::FrozenGas{Float64}; N::Int = 256, Tmin = 200.0, Tmax = 2400.0)
    N ≥ 2 || throw(ArgumentError("tabulate needs N ≥ 2 nodes, got N = $N"))
    Tmin < Tmax || throw(ArgumentError("tabulate needs Tmin < Tmax"))
    # h → T
    hmin, hmax = h(gas, Tmin), h(gas, Tmax)
    dh = (hmax - hmin) / (N - 1)
    Th = [T_of_h(gas, hmin + (i - 1) * dh; Tguess = 0.5 * (Tmin + Tmax)) for i = 1:N]
    Mh = [1 / cp(gas, T) for T in Th]
    # s0 → T
    s0min, s0max = s0(gas, Tmin), s0(gas, Tmax)
    ds0 = (s0max - s0min) / (N - 1)
    Ts = [_T_of_s0(gas, s0min + (i - 1) * ds0, Tmin, Tmax) for i = 1:N]
    Ms = [T / cp(gas, T) for T in Ts]
    TabulatedGas(gas, hmin, hmax, 1 / dh, dh, Th, Mh, s0min, s0max, 1 / ds0, ds0, Ts, Ms)
end

# Cubic Hermite evaluation on a uniform table (x0, invdx, dx, nodes Tn,
# node slopes Mn in dT/dx). Caller guarantees x ∈ [x0, x0 + (N-1)·dx];
# the clamp only guards the floor() at the exact top edge.
@inline function _hermite(x0, invdx, dx, Tn::Vector{Float64}, Mn::Vector{Float64}, x)
    u = (x - x0) * invdx
    i = clamp(floor(Int, u), 0, length(Tn) - 2) # 0-based cell index
    t = u - i
    @inbounds begin
        T0 = Tn[i+1]
        T1 = Tn[i+2]
        m0 = Mn[i+1] * dx
        m1 = Mn[i+2] * dx
    end
    t2 = t * t
    t3 = t2 * t
    (2t3 - 3t2 + 1) * T0 + (t3 - 2t2 + t) * m0 + (-2t3 + 3t2) * T1 + (t3 - t2) * m1
end

@inline _hermite_h(tg::TabulatedGas, hspec) =
    _hermite(tg.hmin, tg.invdh, tg.dh, tg.Th, tg.Mh, hspec)
@inline _hermite_s0(tg::TabulatedGas, starget) =
    _hermite(tg.s0min, tg.invds0, tg.ds0, tg.Ts, tg.Ms, starget)

"""
    T_of_h(tg::TabulatedGas, hspec)

Temperature [K] at which the gas has specific enthalpy `hspec` [J/kg] —
exact, table-seeded: a Hermite-interpolated seed followed by the same Newton
iteration and convergence criterion as the `FrozenGas` method (relative
tolerance 1e-12 on the temperature step, ≤ 30 iterations, typically 1–2
polish steps from the seed). If `hspec` lies outside the tabulated range the
call falls back to the plain `FrozenGas` Newton solve and still returns the
exact answer — no extrapolation. Pure and zero-allocation.
"""
function T_of_h(tg::TabulatedGas, hspec::AbstractFloat)
    if !(tg.hmin ≤ hspec ≤ tg.hmax)
        return T_of_h(tg.gas, hspec) # out of table range: exact cold-start solve
    end
    gas = tg.gas
    T = _hermite_h(tg, hspec)
    for _ = 1:NEWTON_MAXITER
        dT = (hspec - h(gas, T)) / cp(gas, T)
        T += dT
        if abs(dT) ≤ NEWTON_RTOL * abs(T)
            return T
        end
    end
    error("T_of_h did not converge for hspec = $hspec (last T = $T)")
end

# Convenience: integers/rationals promote to float, as the untyped FrozenGas
# methods accept. ForwardDiff Duals dispatch to the (more specific) rules in
# the package extension — they never reach the float Newton loop here.
T_of_h(tg::TabulatedGas, hspec::Real) = T_of_h(tg, float(hspec))

"""
    T_isentropic(tg::TabulatedGas, T1, PR; ηp=1.0)

Temperature [K] after an ideal compression/expansion from `T1` [K] by
pressure ratio `PR`, solving `s0(T2) = s0(T1) + R·ln(PR)/ηp` — exact,
table-seeded: the target s0 is computed from the wrapped gas's exact
polynomials, the s0 → T Hermite table provides the seed, and the same
Newton iteration and convergence criterion as the `FrozenGas` method
finishes the solve (relative tolerance 1e-12, ≤ 30 iterations, typically
1–2 polish steps). `ηp` is the polytropic efficiency, as in the
`FrozenGas` method. If the target s0 lies outside the tabulated range the
call falls back to the plain `FrozenGas` Newton solve and still returns
the exact answer — no extrapolation. Pure and zero-allocation.
"""
function T_isentropic(tg::TabulatedGas, T1::AbstractFloat, PR::AbstractFloat; ηp = 1.0)
    gas = tg.gas
    target = s0(gas, T1) + gas.R * log(PR) / ηp
    if !(tg.s0min ≤ target ≤ tg.s0max)
        return T_isentropic(gas, T1, PR; ηp = ηp) # out of table range: exact cold start
    end
    T = _hermite_s0(tg, target)
    for _ = 1:NEWTON_MAXITER
        dT = (target - s0(gas, T)) * T / cp(gas, T) # ds0/dT = cp/T
        T += dT
        if abs(dT) ≤ NEWTON_RTOL * abs(T)
            return T
        end
    end
    error("T_isentropic did not converge for T1 = $T1, PR = $PR (last T = $T)")
end

T_isentropic(tg::TabulatedGas, T1::Real, PR::Real; ηp = 1.0) =
    T_isentropic(tg, float(T1), float(PR); ηp = ηp)

"""
    T_of_h_interp(tg::TabulatedGas, hspec)

Temperature [K] at specific enthalpy `hspec` [J/kg] from the cubic-Hermite
h → T table alone — **the opt-in approximate tier**: no Newton polish.
Accuracy over the tabulated range at the default N = 256: |ΔT/T| ≲ 1e-9
(measured 5.8e-10 for dry air over T ∈ [250, 2200] K). Throws a
`DomainError` if `hspec` is outside the tabulated range — this tier never
falls back and never extrapolates. For the exact answer use
`T_of_h(tg, hspec)`. Pure and zero-allocation.
"""
function T_of_h_interp(tg::TabulatedGas, hspec)
    if !(tg.hmin ≤ hspec ≤ tg.hmax)
        throw(DomainError(hspec,
            "T_of_h_interp: hspec outside the tabulated range " *
            "[$(tg.hmin), $(tg.hmax)] J/kg; use T_of_h(tg, hspec) for the " *
            "exact out-of-range solve or rebuild the table with a wider [Tmin, Tmax]"))
    end
    _hermite_h(tg, hspec)
end

"""
    T_isentropic_interp(tg::TabulatedGas, T1, PR; ηp=1.0)

Temperature [K] after an ideal compression/expansion from `T1` [K] by
pressure ratio `PR` (target `s0(T2) = s0(T1) + R·ln(PR)/ηp`, computed from
the exact polynomials) from the cubic-Hermite s0 → T table alone — **the
opt-in approximate tier**: no Newton polish. Accuracy over the tabulated
range at the default N = 256: |ΔT/T| ≲ 2e-9 (measured 1.4e-9 for dry air
over T ∈ [250, 2200] K). Throws a `DomainError` if the target s0 is
outside the tabulated range — this tier never falls back and never
extrapolates. For the exact answer use `T_isentropic(tg, T1, PR; ηp)`.
Pure and zero-allocation.
"""
function T_isentropic_interp(tg::TabulatedGas, T1, PR; ηp = 1.0)
    gas = tg.gas
    target = s0(gas, T1) + gas.R * log(PR) / ηp
    if !(tg.s0min ≤ target ≤ tg.s0max)
        throw(DomainError(target,
            "T_isentropic_interp: target s0 for T1 = $T1, PR = $PR, ηp = $ηp " *
            "is outside the tabulated range [$(tg.s0min), $(tg.s0max)] J/kg/K; " *
            "use T_isentropic(tg, T1, PR; ηp) for the exact out-of-range solve " *
            "or rebuild the table with a wider [Tmin, Tmax]"))
    end
    _hermite_s0(tg, target)
end
