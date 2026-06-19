# ADR-0002: FrozenGas is the architecture; Gas1D deprecated; no name recycling

Date: 2026-06-11
Status: accepted

## Context

ADR-0001 introduced the pure immutable `FrozenGas` core additively. Measured
head-to-head (benchmark/arch_comparison/): 7.8× on set-T-read-all, 4.3× on
isentropic compression (49→0 allocs), 1.7× on enthalpy inversion, and
ForwardDiff capability the mutable path structurally cannot have. The
question this ADR settles: do the mutable `Gas{N}`/`Gas1D` types survive,
and does `FrozenGas` take over the name `Gas`?

## Decision

1. **`Gas1D` is deprecated** (depwarn on construction as of now; removal in
   v2.0). It fails the deletion test against `FrozenGas`: same equivalent-
   polynomial concept, strictly worse interface (mutable cache, Float64
   pins, ~8× slower, AD-blocked). Migration: `Gas1D(sp)` + `g.T = T` +
   `g.cp` → `FrozenGas(sp)` + `cp(g, T)` (or `props`).
2. **`Gas{N}` is retained for now** — not as a thermodynamics hot path but
   as the *composition workspace* behind combustion, mixing, and humidity.
   It becomes removable only when pure FrozenGas-producing constructors
   cover those jobs (`products(sys, FAR)`, frozen-gas mixing, humid-air
   construction). Once that layer lands, `Gas{N}` is demoted to an
   unexported interactive convenience or deprecated in turn.

   > **Update (ADR-0008, 2026-06-19):** that condition is now met. The
   > pure constructors cover all three jobs (`Vitiator`/`products`,
   > `mix`, `humid_air`), and the composition-workspace role is served by
   > the immutable `FrozenGas.X` (a mole-fraction `SVector` the gas now
   > carries). `Gas{N}` therefore retires **together with `Gas1D`**; the
   > deletion itself is the ADR-0007 v2.0-final step.
3. **The name `FrozenGas` is permanent. It does not become `Gas` in v2.0.**
   Rationale: (a) recycling a name across a semantics flip
   (mutable → immutable) makes every pre-2.0 snippet, notebook, and doc
   silently wrong — a new name fails loudly instead; (b) "frozen
   composition" is the physics contract (no dissociation/equilibrium), and
   the name documents it. An equilibrium gas model is not expected soon
   (explicitly out of scope for the PowerCycles ladder); the no-recycling
   argument stands on its own.

## Consequences

- The first release PowerCycles depends on presents `FrozenGas` + pure
  property functions + inversions as *the* API; `Gas`/`Gas1D` are
  documented as legacy.
- `Base.depwarn` fires in `Gas1D` constructors; existing tests keep passing
  (depwarn is non-fatal and quiet by default). Legacy `Gas1D` tests remain
  until removal so the deprecated path stays correct while it exists.
- Future architecture reviews should not re-propose renaming `FrozenGas` to
  `Gas`, nor resurrect a mutable hot-path gas type.
