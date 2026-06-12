# ADR-0001: Pure immutable gas core instead of mutable cached state

Date: 2026-06-11
Status: accepted

## Context

The original design centers on mutable `Gas`/`Gas1D` structs: setting `gas.T`
triggers a hidden recomputation of cp/h/ϕ into cached fields, supported by a
`Tarray::MVector{8,Float64}` cache of temperature powers. The stated rationale
was speed (avoid recomputing T powers and `log(T)` across property reads).

Benchmarked prototypes (2026-06-11, Apple Silicon, Julia 1.11/1.12; code kept
in `claude_sandbox/proto_props/` and `claude_sandbox/proto_deriv/`) showed:

- Full mutable update + 3 property reads: **~80 ns** (views over heap `Vector`
  coefficients, `dot`, MVector write-back, always-computed `cp_T`).
- Pure immutable struct (`SVector` coefficients), `props(gas, T)` returning
  (cp, h, ϕ) with shared powers and a single `log(T)`: **~11 ns** — 5–6×
  faster. `log(T)` alone is ~4 ns (the hardware floor).
- The cache only wins when the same property is re-read ≥ ~10× per
  temperature change (cached field read ~0.7 ns).
- ForwardDiff: analytic Dual rules (`Dual(h, cp·∂T)`) are ~4× faster than
  generic Dual arithmetic at 12 partials; implicit-function-theorem rules for
  Newton inversions are ~2× faster than differentiating the loop and agree to
  machine precision.
- The `MVector{8,Float64}` cache and `::Float64` pins hard-block
  `ForwardDiff.Dual`, which downstream solvers (PowerCycles.jl) require.

## Decision

1. Introduce `FrozenGas`: immutable, `isbits`, mass-specific equivalent NASA-9
   coefficients in `SVector{9,Float64}`s, constructed once from a
   species/composite species (construction may consult the global species
   database; property calls may not).
2. All property functions are **pure functions of (gas, T)** — generic over
   `Real`, zero-allocation, thread-safe by construction. `props(gas, T)`
   returns (cp, h, ϕ) sharing the temperature powers and the single `log`.
3. Inversions (`T_of_h`, `T_isentropic`) are deterministic bounded Newton
   solves with documented tolerance, pure in `(gas, args)`.
4. Derivatives are analytic, delivered through a **ForwardDiff package
   extension** (`weakdeps`): Dual dispatches use closed forms
   (dh/dT = cp, dϕ/dT = cp/T) and implicit-function-theorem rules for
   inversions. The base package stays ForwardDiff-free; generic code is the
   always-correct fallback.
5. Mutable `Gas`/`Gas1D` remain unchanged (additive change) as the stateful
   convenience layer.

## Consequences

- Hot-path performance improves ~5–6× while gaining differentiability and
  thread safety; there is no speed/purity trade-off to manage.
- Julia compat moves to ≥ 1.10 (package extensions; also required by
  downstream).
- Future architecture reviews should not re-propose caching temperature
  powers in mutable state for performance: it was measured slower than pure
  recomputation for realistic call patterns.
- Likewise, do not re-propose hand-fusing the three Horner kernels in
  `props` into a shared power-basis kernel (measured 2026-06-11, prototype
  in `claude_sandbox/proto_fused/`): every fused/prescaled variant was
  4–16% *slower* because LLVM already CSEs the shared `1/T` divides across
  the inlined kernels and SLP-vectorizes the h/s0 chains into `<2 x double>`
  SIMD lanes, which explicit power-basis code defeats. The per-call floor is
  `log(T)` (~5 of the 9 ns), untouchable by polynomial rearrangement.
