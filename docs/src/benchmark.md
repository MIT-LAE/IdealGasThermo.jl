# [Performance](@id performance)

The immutable pure core ([`FrozenGas`](@ref) and friends — see [The immutable
pure core](@ref frozengas-page)) is built for speed, zero allocations, and
differentiability.

!!! note "Read the ratios, not the nanoseconds"
    The absolute times below are medians on one reference machine (Apple Silicon,
    Julia 1.12) and **will differ on another CPU or CI runner**. What is robust
    across machines is what these numbers are here to show: the **speedup ratios**
    (set by the algorithm — no cache to maintain, one shared ``\log T``, scalar
    derivative scaling) and the **allocation counts**, which are exact and
    machine-independent (0 vs. 49 allocations is a structural fact, not a timing).
    The harness that produced them lives in `benchmark/arch_comparison/`.

Two ideas do the work:

1. A fixed-composition mixture is collapsed to **one equivalent NASA-9
   polynomial** (see [Representing mixtures with fixed composition](@ref
   gas1dthermo)), so property cost is *independent of the number of species*.
2. The substance is **immutable and `isbits`**, so there is no cache to maintain
   and no heap allocation on the hot path.

## Property reads

A temperature sweep reading the properties at each point. The legacy mutable
types (`Gas{N}`, `Gas1D`) recompute a cache on every ``T`` assignment; the
`FrozenGas` path is a pure function of ``(\text{gas}, T)``.

| workload | legacy `Gas{N}` | legacy `Gas1D` | `FrozenGas` | speedup |
|---|---:|---:|---:|---:|
| read ``c_p, h, s^0`` | 87 ns | 88 ns | **9.3 ns** | ~9.5× |
| read ``c_p`` only | 81 ns | 76 ns | **0.92 ns** | ~80× |

All three are allocation-free; the win is the removed cache bookkeeping and the
shared ``\log T``.

## Inversions and processes

Inverting ``h(\text{gas}, T)`` for ``T`` (a bounded Newton solve), and an
isentropic compression (a process verb returning a new state).

| workload | legacy | `FrozenGas` | notes |
|---|---:|---:|---|
| enthalpy inversion `T_from_h` | ~270 ns | **146 ns** | ~1.8× |
| isentropic `compress` | ~790 ns, **49 allocs** | **147 ns, 0 allocs** | ~5× and allocation-free |

The mutable `compress` allocated 49 times per call; the pure-core process verb
returns an `isbits` [`GasState`](@ref) with no heap traffic.

## Derivatives

The core is generic over `Real`, so ForwardDiff differentiates everything; the
mutable layer **cannot be differentiated at all** (its `::Float64` field pins
reject dual numbers). A package extension adds analytic rules — see
[Thermodynamic derivatives](@ref derivatives) for the math.

Recall that a forward-mode dual number with ``N`` *partials* differentiates with
respect to ``N`` inputs at once (e.g. ``N`` cycle design variables); the gradient
``\nabla f \in \RR^N`` is produced in a single evaluation. The question is
how cost scales with ``N``:

| ``dh/dT`` cost | ``N=1`` | ``N=8`` | ``N=12`` |
|---|---:|---:|---:|
| generic dual arithmetic | 11.8 ns | 26.3 ns | 31.7 ns |
| analytic rule (extension) | **11.2 ns** | **11.3 ns** | **11.5 ns** |

Generic dual arithmetic carries the length-``N`` partials tuple through every
operation, so its cost grows with ``N``. The analytic rule computes the scalar
``c_p`` once and scales the whole tuple
(``h(\text{gas}, x + \dot{x}) = h + c_p\,\dot{x}``), so it is **flat in ``N``**.

For the Newton inversions, differentiating the loop versus applying the
implicit-function-theorem rule (at ``N = 8`` partials):

| `T_from_h`, ``N=8`` | time |
|---|---:|
| dual through the Newton loop | 290 ns |
| IFT rule (extension) | **154 ns** |

The IFT rule is ~2× faster *and* exact — the solver tolerance does not enter the
derivative.

## Independent of the number of species

Because a fixed-composition mixture is represented by a single equivalent
polynomial, property cost does **not** grow with the number of constituent
species — the central trick the `FrozenGas`/`Gas1D` lumped representation shares.
The per-species `Gas{N}` path, by contrast, scales with the species count.

We sweep temperature 100 times, reading ``c_p, h, \phi`` and ``dc_p/dT`` each
time:

```julia
TT = rand(200.:600., 100)

function benchmark_Gas(TT::AbstractVector, gas::AbstractGas)
    @views for i in eachindex(TT)
        gas.T = TT[i]
        gas.cp; gas.ϕ; gas.h; gas.cp_T
    end
end
```

For a small mixture, the lumped (`Gas1D`) path is already ~3× the per-species
(`Gas`) path:

```julia-repl
julia> @benchmark benchmark_Gas($TT, $gas1D)   # lumped, 1 equivalent species
  median: 2.884 μs, 0 allocations

julia> @benchmark benchmark_Gas($TT, $gas)     # per-species mixture
  median: 7.667 μs, 0 allocations
```

Grow the mixture to 10 components and the per-species path scales with it, while
the lumped path is unchanged:

```julia-repl
julia> @benchmark benchmark_Gas($TT, $gas)     # 10 components
  median: 17.875 μs, 0 allocations

julia> @benchmark benchmark_Gas($TT, $gas1D)   # still 1 equivalent species
  median: 2.903 μs, 0 allocations
```

The `FrozenGas` pure core uses the same equivalent-polynomial representation, so
it inherits this species-count independence — on top of the ~9.5× constant
factor and zero allocations shown above.
