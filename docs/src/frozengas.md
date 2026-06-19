# [The immutable pure core](@id frozengas-page)

The headline API of `IdealGasThermo.jl` is an **immutable, allocation-free**
representation of a gas. A substance is a set of property curves
([`FrozenGas`](@ref), or [`FastFrozenGas`](@ref) with tabulated inverses); a
thermodynamic point pairs a substance with ``T`` and ``P`` ([`GasState`](@ref));
process and flow *verbs* map states to states. Combustion ([`Vitiator`](@ref) /
[`products`](@ref)), mixing ([`mix`](@ref)) and humidity ([`humid_air`](@ref))
all produce immutable gases.

New here? Start with the [Getting started tutorial](@ref tutorial-getting-started).

!!! note "Specific heat is exported as `cₚ` / `c_p`"
    The specific-heat accessor is exported under the interchangeable aliases
    `cₚ` and `c_p` (the bare name `cp` would shadow `Base.cp`); the ratio of
    specific heats is exported as `gamma` and the Unicode alias `γ`.

## [Substance vs. state](@id frozengas-substance-state)

A `FrozenGas` is a **substance**, not a parcel of gas. It stores the
fixed-composition NASA-9 property curves — an *equivalent* set of polynomial
coefficients for the mixture (see [Representing mixtures with fixed
composition](@ref gas1dthermo)) — and nothing about the current conditions.
**Temperature is an argument** to every property function, never a field:

```math
c_p = c_p(\text{gas}, T), \qquad h = h(\text{gas}, T), \qquad s^0 = s^0(\text{gas}, T).
```

The legacy mutable `Gas`/`Gas1D` types did the opposite: they stored ``T`` on the
struct and cached ``c_p, h, \phi`` into fields, so assigning `gas.T = …` silently
recomputed the cache. That coupling — a value that is both substance *and* state,
mutated in place — is exactly what the pure core removes.

When you *do* need the conditions, a [`GasState`](@ref) is the immutable
``(\text{substance}, T, P)`` value that carries them. Process verbs
(`compress`, `expand`, `add_heat`, …) take a state and return a **new** state;
the substance is never touched.

## Why immutability is better

### 1. Speed — there is no cache to maintain

Property reads are pure functions of ``(\text{gas}, T)`` over an `SVector` of
coefficients, with the temperature powers and the single ``\log T`` shared across
``c_p, h, s^0``. Reading all three is **~9.5× faster** than the cached mutable
path, and reading ``c_p`` alone is under a nanosecond (**~85×**). Because the
mixture is collapsed to one equivalent polynomial, the cost is **independent of
the number of constituent species**. See [Performance](@ref performance) for the
measured table.

### 2. Zero allocations on the hot path

`FrozenGas{Float64}` is `isbits`, so it lives on the stack and inside other
immutable values (`NamedTuple`s, `GasState`s) with no heap traffic. Process
verbs return new `isbits` states; an isentropic compression that allocated **49
times** on the mutable layer is **allocation-free** here.

### 3. Differentiable — gradients for (almost) free

The core is generic over `Real`, so [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl)
propagates derivatives through every property and process. The mutable layer
**cannot** do this at all: its `::Float64` field pins reject `Dual` numbers. A
package extension supplies analytic closed forms (``dh/dT = c_p``) and
implicit-function-theorem rules for the Newton inversions, so differentiation
never re-runs a solver. See [Thermodynamic derivatives](@ref derivatives).

### 4. Thread-safe by construction

An immutable substance has no mutable shared state, so the same gas can be read
from any number of threads with no locking and no risk of one task's
``T`` assignment corrupting another's read.

### 5. Value semantics — no action at a distance

Every verb returns a new state, so a station you saved stays exactly as it was.
The bug where a downstream step advances the wrong pressure or temperature "rail"
— easy to write when state lives in a mutable parcel passed by reference —
cannot be expressed: the pressure travels *inside* the value the verb returns.

### 6. Self-describing composition

A `FrozenGas` carries its source mole-fraction vector ``X`` (with
``\sum_i X_i = 1``), so it remembers what it is made of. That is what lets
combustion products be re-mixed with bypass air, or a mixture be re-burned in an
afterburner — the composition is never thrown away after lumping.

---

## API reference

### Substances and properties

```@autodocs
Modules = [IdealGasThermo]
Pages = ["frozengas.jl", "fastfrozengas.jl"]
```

### Thermodynamic state and process verbs

```@autodocs
Modules = [IdealGasThermo]
Pages = ["gasstate.jl"]
```

### Gas dynamics

```@autodocs
Modules = [IdealGasThermo]
Pages = ["flow.jl"]
```

### Combustion, mixing and humidity

```@autodocs
Modules = [IdealGasThermo]
Pages = ["vitiator.jl", "mix.jl", "humidity.jl"]
```
