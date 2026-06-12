# IdealGasThermo.jl — Domain Context

Vocabulary used in code, tests, docs, and architecture discussions. Use these
terms exactly.

## Domain terms

- **Species** — a single chemical species with NASA-9 polynomial coefficients
  (`alow`/`ahigh` over two temperature intervals split at **Tmid** = 1000 K),
  molecular weight `MW` [g/mol], and formation enthalpy `Hf` [J/mol] at
  298.15 K. Loaded from `data/thermo.inp` (NASA ThermoBuild format).
- **Composite species** — a fixed-composition gas mixture represented as a
  *pseudo-species*: one equivalent NASA-9 coefficient set, precomputed at
  construction as the mole-fraction-weighted sum of constituent coefficients
  (entropy of mixing folded into the integration constant b₂). The package's
  core idea; see `docs/src/idealgasthermo.md`.
- **FrozenGas** — the immutable, `isbits` form of a composite species: the
  *pure property core*. Mass-specific (J/kg-based) equivalent coefficients in
  `SVector`s. All property functions of a `FrozenGas` are pure functions of
  `(gas, T)` — no state, no globals on the hot path, generic over `Real`.
- **Entropy complement (ϕ / s0)** — φ(T) = ∫cp/T dT from the standard state;
  `s(T, P) = ϕ(T) − R·ln(P/Pstd)`. The temperature-only part of entropy.
- **Enthalpy datum** — enthalpies are **formation-inclusive** (CEA-style):
  `h(gas, 298.15 K)` equals the mixture's mass-specific formation enthalpy,
  not zero. Sensible enthalpy from 298.15 K is `h(gas, T) − h(gas, 298.15)`.
- **Vitiated mixture** — combustion products of a fuel + oxidizer at a given
  **FAR** (fuel–air mass ratio), frozen composition, complete combustion.
- **Combustor** — a precomputed fuel + oxidizer combustion system: dense
  per-species `SVector`/`SMatrix` data built once from the species database.
  The pure, allocation-free replacement for the Dict-based
  `vitiated_species` path on the hot path.
- **products** — `products(sys::Combustor, FAR) -> FrozenGas`: the
  combustion-product gas at a given FAR. Pure, zero-allocation, smooth in
  FAR (ForwardDiff through FAR works; `FrozenGas{TF}` widens its eltype).
- **Inversion** — solving a property relation backwards for temperature:
  `T_of_h` (given h) and `T_isentropic` (given T1 and pressure ratio, along
  an isentrope with optional polytropic efficiency).

## Architecture terms (see docs/adr/)

- The **pure core** (`FrozenGas` + property functions + inversions) is the
  deep module everything else composes over.
- The mutable `Gas`/`Gas1D` types are the legacy *stateful convenience layer*;
  they are kept API-stable but are not the hot path. **`Gas1D` is deprecated**
  (ADR-0002): use `FrozenGas`. `Gas{N}` remains only as the composition
  workspace until pure FrozenGas-producing constructors replace it.
- `FrozenGas` keeps its name permanently — it is not renamed to `Gas` in
  v2.0 (ADR-0002: no name recycling across a semantics flip; "frozen" names
  the no-dissociation contract).
- Derivatives are **analytic first**: ForwardDiff `Dual` support is provided by
  a package extension that dispatches to closed-form derivatives
  (dh/dT = cp, dϕ/dT = cp/T) and implicit-function-theorem rules for
  inversions — never by differentiating a Newton loop.
