"""
    FastFrozenGas{mode}

A [`FrozenGas`](@ref) bundled with two precomputed cubic-Hermite *inverse*
tables on uniform grids — h → T and s0 → T — that accelerate the
[`temperature`](@ref) inversions. The `mode` type parameter selects the
accuracy tier, so a `FastFrozenGas` drops into any code written against the
gas interface with no call-site changes:

```julia-repl
julia> air  = FrozenGas(DryAir);
julia> fast = FastFrozenGas(air);              # mode = :seeded (default)
julia> fast = FastFrozenGas(air, mode = :fast) # pure-lookup tier

julia> temperature(fast, h = 5.0e5)            # same verb for every gas
```

Only the inverses are accelerated. The forward property functions
(`cp`, `h`, `s0`, `props`, `gamma`, `R`, `pressure_ratio`) forward to the
wrapped `FrozenGas` unchanged: the forward functions never change; both
tables sample the same NASA-9 polynomials.

The two modes:

- **`:seeded` (default — exact)**: the table provides only the *seed* for
  the same bounded Newton iteration and convergence criterion as
  `FrozenGas` (relative tolerance 1e-12 on the temperature step, ≤ 30
  iterations, typically 1–2 polish steps). The answer satisfies the
  package's documented inversion contract. Targets outside the tabulated
  range bypass the table entirely and use the plain cold-start `FrozenGas`
  solve — never extrapolation.
- **`:fast` (opt-in — approximate)**: the Hermite table is evaluated
  directly, no polish. Accuracy |ΔT/T| ≲ 2e-9 over the table range at the
  default N = 256 (measured for dry air over T ∈ [250, 2200] K: 5.8e-10
  for the h table, 1.4e-9 for the s0 table). Out-of-range targets throw a
  `DomainError` — the approximate tier has no silent fallback and never
  extrapolates.
"""
struct FastFrozenGas{M}
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
    R(fg::FastFrozenGas)

Specific gas constant [J/kg/K] of the wrapped [`FrozenGas`](@ref).
"""
@inline R(fg::FastFrozenGas) = R(fg.gas)
@inline cp(fg::FastFrozenGas, T) = cp(fg.gas, T)
@inline h(fg::FastFrozenGas, T::Real) = h(fg.gas, T)
@inline s0(fg::FastFrozenGas, T) = s0(fg.gas, T)
@inline gamma(fg::FastFrozenGas, T) = gamma(fg.gas, T)
@inline speed_of_sound(fg::FastFrozenGas, T) = speed_of_sound(fg.gas, T)
@inline props(fg::FastFrozenGas, T) = props(fg.gas, T)
@inline pressure_ratio(fg::FastFrozenGas, T1, T2) = pressure_ratio(fg.gas, T1, T2)

# Exact Newton inverse of s0 (table construction only; same tolerance as
# the temperature inversions). s0 is strictly monotonic in T, so for
# targets inside [s0(Tlo), s0(Thi)] the root lies in [Tlo, Thi]; clamping
# each iterate to that bracket keeps the log(T) evaluations valid.
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
    FastFrozenGas(gas::FrozenGas; mode::Symbol=:seeded, N::Int=256,
                  Tmin=200.0, Tmax=2400.0)

Build a [`FastFrozenGas{mode}`](@ref FastFrozenGas): precompute
cubic-Hermite inverse tables h → T and s0 → T on uniform N-node grids
spanning `T ∈ [Tmin, Tmax]` [K]. `mode` is `:seeded` (exact, default) or
`:fast` (pure lookup). Node temperatures come from the exact `FrozenGas`
Newton inversions; node slopes are the exact analytic inverses
dT/dh = 1/cp and dT/ds0 = T/cp, so the interpolant is C¹ and matches the
polynomials to |ΔT/T| ≲ 2e-9 at the default N = 256. Construction
allocates (~4 KB per table); all subsequent inversion calls are pure and
zero-allocation.
"""
function FastFrozenGas(
    gas::FrozenGas{Float64};
    mode::Symbol = :seeded,
    N::Int = 256,
    Tmin = 200.0,
    Tmax = 2400.0,
)
    mode === :seeded || mode === :fast ||
        throw(ArgumentError("FastFrozenGas mode must be :seeded or :fast, got :$mode"))
    N ≥ 2 || throw(ArgumentError("FastFrozenGas needs N ≥ 2 nodes, got N = $N"))
    Tmin < Tmax || throw(ArgumentError("FastFrozenGas needs Tmin < Tmax"))
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
    FastFrozenGas{mode}(
        gas, hmin, hmax, 1 / dh, dh, Th, Mh, s0min, s0max, 1 / ds0, ds0, Ts, Ms,
    )
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

@inline _hermite_h(fg::FastFrozenGas, hspec) =
    _hermite(fg.hmin, fg.invdh, fg.dh, fg.Th, fg.Mh, hspec)
@inline _hermite_s0(fg::FastFrozenGas, starget) =
    _hermite(fg.s0min, fg.invds0, fg.ds0, fg.Ts, fg.Ms, starget)

# ---- :seeded mode — exact, table-seeded Newton -----------------------------

function T_of_h(fg::FastFrozenGas{:seeded}, hspec::AbstractFloat)
    if !(fg.hmin ≤ hspec ≤ fg.hmax)
        return T_of_h(fg.gas, hspec) # out of table range: exact cold-start solve
    end
    gas = fg.gas
    T = _hermite_h(fg, hspec)
    for _ = 1:NEWTON_MAXITER
        dT = (hspec - h(gas, T)) / cp(gas, T)
        T += dT
        if abs(dT) ≤ NEWTON_RTOL * abs(T)
            return T
        end
    end
    error("temperature(fg, h = $hspec) did not converge (last T = $T)")
end

function T_isentropic(fg::FastFrozenGas{:seeded}, T1::AbstractFloat, PR::AbstractFloat; ηp = 1.0)
    gas = fg.gas
    target = s0(gas, T1) + gas.R * log(PR) / ηp
    if !(fg.s0min ≤ target ≤ fg.s0max)
        return T_isentropic(gas, T1, PR; ηp = ηp) # out of table range: exact cold start
    end
    T = _hermite_s0(fg, target)
    for _ = 1:NEWTON_MAXITER
        dT = (target - s0(gas, T)) * T / cp(gas, T) # ds0/dT = cp/T
        T += dT
        if abs(dT) ≤ NEWTON_RTOL * abs(T)
            return T
        end
    end
    error("T_isentropic(fg, T1 = $T1, PR = $PR) did not converge (last T = $T)")
end

# ---- :fast mode — pure Hermite lookup, loud out of range -------------------

function T_of_h(fg::FastFrozenGas{:fast}, hspec::AbstractFloat)
    if !(fg.hmin ≤ hspec ≤ fg.hmax)
        throw(DomainError(hspec,
            "temperature(fg, h = ...): h outside the tabulated range " *
            "[$(fg.hmin), $(fg.hmax)] J/kg; use the :seeded mode for the " *
            "exact out-of-range solve or rebuild with a wider [Tmin, Tmax]"))
    end
    _hermite_h(fg, hspec)
end

function T_isentropic(fg::FastFrozenGas{:fast}, T1::AbstractFloat, PR::AbstractFloat; ηp = 1.0)
    gas = fg.gas
    target = s0(gas, T1) + gas.R * log(PR) / ηp
    if !(fg.s0min ≤ target ≤ fg.s0max)
        throw(DomainError(target,
            "T_isentropic(fg, T1 = $T1, PR = $PR, ηp = $ηp): target s0 outside " *
            "the tabulated range [$(fg.s0min), $(fg.s0max)] J/kg/K; use the " *
            ":seeded mode for the exact out-of-range solve or rebuild with a " *
            "wider [Tmin, Tmax]"))
    end
    _hermite_s0(fg, target)
end

# Convenience: integers/rationals promote to float, as the untyped FrozenGas
# methods accept. ForwardDiff Duals dispatch to the (more specific) rules in
# the package extension — they never reach the float Newton loop here.
T_of_h(fg::FastFrozenGas, hspec::Real) = T_of_h(fg, float(hspec))
T_isentropic(fg::FastFrozenGas, T1::Real, PR::Real; ηp = 1.0) =
    T_isentropic(fg, float(T1), float(PR); ηp = ηp)
