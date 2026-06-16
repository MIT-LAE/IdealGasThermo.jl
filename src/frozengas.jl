"""
    FrozenGas

Immutable, `isbits` representation of a fixed ("frozen") composition gas:
the pure property core of IdealGasThermo.

Holds a single equivalent NASA-9 coefficient set for the mixture,
**mass-scaled at construction** (coefficients premultiplied by `1000/MW` so
property polynomials evaluate directly in J/kg-based SI units, with no
per-call division or species summation). Construct from any
[`AbstractSpecies`](@ref) (a database [`species`](@ref) or a
[`composite_species`](@ref)):

```julia-repl
julia> air = FrozenGas(DryAir);
```

All property functions of a `FrozenGas` are pure functions of `(gas, T)`:
no state, no global lookups, zero allocations, generic over `Real`.

Enthalpy datum: formation-inclusive (CEA-style) — `h(gas, 298.15)` is the
mixture's mass-specific formation enthalpy. Sensible enthalpy from 298.15 K
is `h(gas, T) - h(gas, 298.15)`.

The field eltype `TF` is `Float64` for ordinary use (keeping the struct
`isbits`); it widens to e.g. `ForwardDiff.Dual` when the *composition*
itself carries derivative information, as in [`products`](@ref)
differentiated with respect to FAR.
"""
struct FrozenGas{TF<:Real}
    alow::SVector{9,TF}  # mass-scaled coefficients, T < Tmid (1000 K)
    ahigh::SVector{9,TF} # mass-scaled coefficients, T ≥ Tmid
    MW::TF               # molecular weight [g/mol]
    R::TF                # specific gas constant [J/kg/K]
    Hf::TF               # formation enthalpy at 298.15 K [J/mol]
end

function FrozenGas(sp::AbstractSpecies)
    scale = 1000.0 / sp.MW # molar (J/mol) → mass-specific (J/kg)
    FrozenGas(
        SVector{9,Float64}(sp.alow) * scale,
        SVector{9,Float64}(sp.ahigh) * scale,
        sp.MW,
        Runiv / sp.MW * 1000.0,
        sp.Hf,
    )
end

"""
    FrozenGas(X::AbstractVector, name::AbstractString="frozen gas")

Construct from mole fractions `X` ordered as the species database
(`spdict`); `X` must sum to 1. Consults the database once, here — property
calls never do.
"""
FrozenGas(X::AbstractVector, name::AbstractString = "frozen gas") =
    FrozenGas(generate_composite_species(X, name))

"""
    R(gas::FrozenGas)

Specific gas constant [J/kg/K].
"""
R(gas::FrozenGas) = gas.R

# Coefficient set for the NASA-9 interval containing T.
# Tmid == 1000 K always (validated in readThermo.jl); same branch convention
# as Gas1D: alow strictly below 1000 K.
@inline coeffs(gas::FrozenGas, T) = T < 1000.0 ? gas.alow : gas.ahigh

# NASA-9 polynomial kernels in dimensionless (per-R) form, shared between
# the scalar property functions and `props` so the two paths are
# bit-identical.
@inline poly_cp_R(a, T) =
    (a[1] / T + a[2]) / T + a[3] + T * (a[4] + T * (a[5] + T * (a[6] + T * a[7])))
@inline poly_h_R(a, T, lnT) =
    -a[1] / T + a[2] * lnT + a[8] +
    T * (a[3] + T * (a[4] / 2 + T * (a[5] / 3 + T * (a[6] / 4 + T * (a[7] / 5)))))
@inline poly_s0_R(a, T, lnT) =
    (-a[1] / T / 2 - a[2]) / T + a[3] * lnT + a[9] +
    T * (a[4] + T * (a[5] / 2 + T * (a[6] / 3 + T * (a[7] / 4))))

"""
    cp(gas::FrozenGas, T)

Specific heat at constant pressure [J/kg/K] at temperature `T` [K].
Pure, zero-allocation, generic over `Real`.
"""
@inline cp(gas::FrozenGas, T) = Runiv * poly_cp_R(coeffs(gas, T), T)

"""
    h(gas::FrozenGas, T)

Specific enthalpy [J/kg] at temperature `T` [K]. Formation-inclusive
(CEA-style) datum: `h(gas, 298.15)` is the mixture formation enthalpy;
sensible enthalpy is `h(gas, T) - h(gas, 298.15)`.
Pure, zero-allocation, generic over `Real`.
"""
@inline h(gas::FrozenGas, T::Real) = Runiv * poly_h_R(coeffs(gas, T), T, log(T))

"""
    s0(gas::FrozenGas, T)

Standard-state entropy function φ(T) = ∫cp/T dT [J/kg/K] (entropy
complement). Entropy at pressure `P` is `s0(gas, T) - R(gas)*log(P/Pstd)`.
Pure, zero-allocation, generic over `Real`.
"""
@inline s0(gas::FrozenGas, T) = Runiv * poly_s0_R(coeffs(gas, T), T, log(T))

"""
    props(gas::FrozenGas, T)

All temperature-dependent properties in one call, sharing the temperature
powers and the single `log(T)`: returns `(cp = ..., h = ..., s0 = ...)`
([J/kg/K], [J/kg], [J/kg/K]). Equivalent to calling [`cp`](@ref),
[`h`](@ref), [`s0`](@ref) individually, ~2x faster when more than one
property is needed. Pure, zero-allocation, generic over `Real`.
"""
@inline function props(gas::FrozenGas, T)
    a = coeffs(gas, T)
    lnT = log(T)
    (
        cp = Runiv * poly_cp_R(a, T),
        h = Runiv * poly_h_R(a, T, lnT),
        s0 = Runiv * poly_s0_R(a, T, lnT),
    )
end

"""
    gamma(gas::FrozenGas, T)

Ratio of specific heats cp/(cp - R) at temperature `T` [K].
"""
@inline function gamma(gas::FrozenGas, T)
    c = cp(gas, T)
    c / (c - gas.R)
end

"""
    speed_of_sound(gas::FrozenGas, T)

Speed of sound `a = √(γ·R·T)` [m/s] at temperature `T` [K], with
`γ = `[`gamma`](@ref)`(gas, T)` and `R = `[`R`](@ref)`(gas)`. A pure
function of `(gas, T)`: composition and temperature are all it needs — no
pressure, no state. Zero-allocation, generic over `Real`.
"""
@inline speed_of_sound(gas::FrozenGas, T) = sqrt(gamma(gas, T) * gas.R * T)

# Inversion contract (T_of_h, T_isentropic): Newton iteration, relative
# tolerance 1e-12 on the temperature step (well inside the documented
# ≤ 1e-10), at most 30 iterations, deterministic fixed algorithm, errors if
# not converged. dh/dT = cp > 0 makes h strictly monotonic in T, so the
# solve is well-posed over the data's validity range.
const NEWTON_RTOL = 1e-12
const NEWTON_MAXITER = 30

"""
    T_of_h(gas::FrozenGas, hspec; Tguess=500.0)

Temperature [K] at which the gas has specific enthalpy `hspec` [J/kg]
(same formation-inclusive datum as [`h`](@ref)). Deterministic bounded
Newton solve: relative tolerance 1e-12, ≤ 30 iterations; errors if not
converged. Pure and zero-allocation.
"""
function T_of_h(gas::FrozenGas, hspec; Tguess = 500.0)
    T = one(hspec / oneunit(hspec)) * Tguess # promote to eltype of hspec
    for _ = 1:NEWTON_MAXITER
        dT = (hspec - h(gas, T)) / cp(gas, T)
        T += dT
        if abs(dT) ≤ NEWTON_RTOL * abs(T)
            return T
        end
    end
    error("T_of_h did not converge for hspec = $hspec (last T = $T)")
end

"""
    pressure_ratio(gas::FrozenGas, T1, T2)

The pressure ratio P2/P1 across an ideal (isentropic) process taking the gas
from `T1` to `T2` [K]: `exp((s0(T2) - s0(T1))/R)`. Inverse of
[`T_isentropic`](@ref). Pure, zero-allocation, generic over `Real`.
"""
@inline pressure_ratio(gas::FrozenGas, T1, T2) =
    exp((s0(gas, T2) - s0(gas, T1)) / gas.R)

"""
    T_isentropic(gas::FrozenGas, T1, PR; ηp=1.0)

Temperature [K] after an ideal compression/expansion from `T1` [K] by
pressure ratio `PR`, solving `s0(T2) = s0(T1) + R·ln(PR)/ηp` (Newton with a
constant-γ initial guess; relative tolerance 1e-12, ≤ 30 iterations; errors
if not converged). `ηp` is the polytropic efficiency: pass `ηp` for
compression (PR > 1) and `1/ηp` for expansion (PR < 1) conventions, or leave
at 1 for the isentropic relation. Pure and zero-allocation.
"""
function T_isentropic(gas::FrozenGas, T1, PR; ηp = 1.0)
    target = s0(gas, T1) + gas.R * log(PR) / ηp
    T = T1 * PR^(gas.R / cp(gas, T1) / ηp) # constant-γ initial guess
    for _ = 1:NEWTON_MAXITER
        dT = (target - s0(gas, T)) * T / cp(gas, T) # ds0/dT = cp/T
        T += dT
        if abs(dT) ≤ NEWTON_RTOL * abs(T)
            return T
        end
    end
    error("T_isentropic did not converge for T1 = $T1, PR = $PR (last T = $T)")
end

"""
    temperature(gas; h = ...)

The public inversion verb: the temperature [K] at which the gas has
specific enthalpy `h` [J/kg] (same formation-inclusive datum as
[`h`](@ref)).

The verb is identical for every gas flavor — `FrozenGas` (plain Newton),
`FastFrozenGas{:seeded}` (table-seeded Newton, same exact contract), and
`FastFrozenGas{:fast}` (pure table lookup, ≲ 2e-9) — so accelerated gases
drop into existing call sites unchanged. The keyword value may be a
ForwardDiff `Dual`; derivatives use the implicit-function-theorem rules
from the package extension, never differentiation of the Newton loop.

`temperature` inverts a property relation; it does not run a process. The
former isentrope form `temperature(gas, T1 = ..., PR = ...; ηp)` is now an
`ArgumentError` — a polytropic change of state is a *process*, expressed by
the process verbs [`compress`](@ref) and [`expand`](@ref) (both with
pressure ratio ≥ 1; the direction lives in the verb). See ADR-0004.
"""
function temperature(gas; h = nothing, T1 = nothing, PR = nothing, ηp = nothing)
    if T1 !== nothing || PR !== nothing || ηp !== nothing
        throw(
            ArgumentError(
                "temperature no longer takes the isentrope form (T1/PR/ηp): " *
                "that is a process, not an inversion — use " *
                "compress(gas, T1, PR; ηp) or expand(gas, T1, PR; ηp), " *
                "both with pressure ratio ≥ 1",
            ),
        )
    end
    h !== nothing ||
        throw(ArgumentError("temperature: specify h = ... (enthalpy inversion)"))
    return T_of_h(gas, h)
end
