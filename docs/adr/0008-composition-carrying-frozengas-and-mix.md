# ADR-0008: A `FrozenGas` carries its composition; `mix` replaces `Mixer`/`mixed`; the species table is consolidated

Date: 2026-06-19
Status: accepted (extends ADR-0001/0002's pure core; supersedes the `Mixer`/`mixed`
pair added in the `2.0.0-beta1` work; refines the mixer-scope note left open in
CONTEXT.md)

## Context

The `2.0.0-beta1` pure core could **burn** (`Vitiator`/`products`, then named
`Combustor`) and **mix**
(`Mixer`/`mixed`) compositions, but the two could not be **composed**:
`products` returns a `FrozenGas`, and `Mixer` could not consume one — so
combustion products could not be mixed with bypass air (the internally-mixed-flow
turbofan), nor re-burned in an afterburner. The blocker was structural: a
`FrozenGas` stored only the *lumped* equivalent NASA-9 coefficients and discarded
the mole-fraction vector `X` that `products`/`mixed`/`generate_composite_species`
all compute. Two lumped gases cannot be re-merged correctly, because the entropy
of mixing `−Σ Xᵢ ln Xᵢ` folded into b₂ is nonlinear in composition and
unrecoverable once `X` is gone.

Three further facts shaped the decision:

1. **`mixed` did not do the energy balance.** The mixed temperature is a
   mass-averaged enthalpy — intrinsic to mixing — but the beta API deferred it to
   the caller (the "Mixer = composition only" note in CONTEXT.md). That was an
   incomplete API, not a deliberate boundary.
2. **The dense species "table" was triplicated.** The `9×Nspecies` matrix
   `reduce(hcat, spdict.alow)` (used as `A·X` to lump a composition) was rebuilt on
   every `generate_composite_species` call and cached as identical struct fields in
   *both* the combustion system (`Vitiator`) and `Mixer`. There was no single home.
3. **Momentum / mass flow are out of scope here** (ADR-0005): a flow station that
   carries a flow rate, velocity, or area belongs to PowerCycles, not this package.

## Decision

1. **A `FrozenGas` carries its source composition** as a new field
   `X::SVector{Nspecies,TF}` (mole fractions, spdict order, Σ = 1). The gas becomes
   *self-describing* — it remembers what it is made of — which is what makes
   re-mixing and re-burning possible. This extends, but does not violate, the
   substance-vs-state split (ADR-0001): `X` is an immutable constant of the
   composition, like the coefficients/MW/R/Hf; it carries no T or P. The struct
   stays `isbits` for `FrozenGas{Float64}` (now 328 B), and the property hot path
   (`cp`/`h`/`s0`, which never read `X`) is unchanged in speed and allocation.
   `products`/`mixed`/`humid_air` already computed `X` and threw it away; they now
   store it. For the AD path, `X` is naturally `Dual`-valued when produced by
   `products(sys, FAR::Dual)` — the eltype `TF` widens with it, exactly as the
   other fields already do.

2. **Mixing is a free function `mix`, not a precomputed system.** Because each gas
   carries `X`, no `Mixer` object is needed:
   - `mix(a::FrozenGas, b::FrozenGas, mratio) -> FrozenGas` blends the two
     compositions at mass ratio `mratio = mass_b/mass_a` and rebuilds the merged
     gas from the blended `X` (`mix(a, b, 0) == a`; `b` is the `mratio → ∞` limit).
   - `mix(a::GasState, b::GasState, mratio) -> GasState` does the **energy
     balance** as well: the mixed temperature is the mass-averaged total enthalpy
     `h = (h_a + mratio·h_b)/(1 + mratio)` inverted on the merged gas. It requires
     equal stream pressures (an isobaric mixer); a non-isobaric mix needs a
     momentum closure and throws an `ArgumentError` directing the caller to the
     flow layer.

   `Mixer` and `mixed` are **deleted** (they were new in `2.0.0-beta1` and never
   released, so no deprecation shim is owed).

3. **The energy balance is parametrized by mass *ratio*, not mass *flow*.** Mass
   ratio is dimensionless and already the mixing currency; mass-flow bookkeeping
   (computing `mratio` from `ṁ₁`, `ṁ₂`, tracking total `ṁ`) and the momentum
   closure (the mixed *pressure*) stay in PowerCycles, per ADR-0005. Energy and
   momentum are independent conservation laws: the energy balance fixes the mixed
   (stagnation) temperature regardless of momentum.

4. **The dense table is consolidated to one home and one kernel.** A single set of
   module constants in `readThermo.jl` (`SPALOW`/`SPAHIGH`/`SPMW`/`SPHF`) holds the
   static-array repacking of `spdict`, and one internal kernel `_lump_molar(X)`
   forms the equivalent molar coefficients (with the entropy-of-mixing fold). All
   producers — `generate_composite_species`, `FrozenGas(X)`, `products`, `mix` —
   route through it. `products` collapses to `FrozenGas(X)` and `Vitiator` loses
   its four cached table fields; the entropy-of-mixing fold, previously duplicated,
   now lives in exactly one place.

5. **`mix` rebuilds the merged gas from `X_mix` via the shared kernel** (rather
   than recombining the two gases' lumped coefficients). Both are mathematically
   exact (verified bit-identical to the canonical build to ~2.5e-16); rebuilding
   from `X` keeps a single construction path and is the obvious choice now that the
   const table exists.

6. **`Gas{N}` retires together with `Gas1D`** (closing ADR-0002 §2). ADR-0002 kept
   `Gas{N}` only as the "composition workspace behind combustion, mixing, and
   humidity" until pure constructors covered those jobs. They now do
   (`Vitiator`/`products`, `mix`, `humid_air`), and the composition-carrying role
   is served by the immutable `FrozenGas.X` / mole-fraction `SVector`, so the
   condition for retiring `Gas{N}` is met. The actual deletion remains the
   ADR-0007 v2.0-final step. As cleanup that de-risks that deletion,
   `fuelbreakdown` and `reaction_change_molar_fraction` — pure functions reachable
   from the pure-core `Vitiator` — were re-homed out of the legacy `combustion.jl`
   into `vitiator.jl`.

7. **The combustion precompute type is named `Vitiator`, not `Combustor`**
   (`src/vitiator.jl`, renamed from `combustor.jl`). This is the *composition*
   model of combustion — it turns a FAR into a vitiated `FrozenGas` via
   `products`. The noun `Combustor` is reserved for the **hardware component** a
   cycle deck (e.g. PowerCycles) models one abstraction level up — pressure drop,
   efficiency map, bleed, geometry. Exporting `Combustor` from this foundational
   layer would force every downstream file that also imports a cycle library to
   qualify the name; `Vitiator` (an agent noun, parallel to the deleted `Mixer`,
   and tied to the legacy `vitiated_*` vocabulary) leaves the hardware noun free.
   The constructor is split into two methods — a `species` method does the work,
   the `AbstractString` method resolves the fuel name and forwards — replacing the
   `Union{AbstractString,species}` + runtime `isa` branch.

## Consequences

- The internally-mixed-flow turbofan composes directly:
  `mix(GasState(products(comb, FAR), T_core, P), GasState(FrozenGas(DryAir),
  T_bypass, P), BPR)` returns the energy-balanced mixed stagnation state, and the
  mixed gas's `X` feeds an afterburner `Vitiator` or another `mix`. The whole
  chain is zero-allocation and ForwardDiff-through-`FAR`/`mratio`.
- `FrozenGas` is now the single composition currency; there is no separate
  composition object and no `Mixer`. Building a gas from a mole vector
  (`FrozenGas(X)`) is zero-allocation (it was not before — it routed through
  `generate_composite_species`).
- Do not re-introduce a precomputed mixing struct, and do not give `mix` a
  mass-flow or momentum argument — that boundary belongs to PowerCycles (ADR-0005).
- Do not re-duplicate the species table: `SPALOW`/`SPAHIGH`/`SPMW`/`SPHF` and
  `_lump_molar` are the one home.
