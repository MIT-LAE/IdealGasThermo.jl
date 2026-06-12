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
- **Mixer** — a precomputed two-stream mixing system (per-stream
  mass-fraction `SVector`s + dense species data), built once from two
  composition-bearing inputs. The pure replacement for the Dict-based
  composition step of the legacy `gas_mixing`.
- **mixed** — `mixed(sys::Mixer, mratio) -> FrozenGas`: the merged gas at
  mass ratio `mratio = mass₂/mass₁`, via the mass-fraction law of mixtures
  `Y = (Y₁ + mratio·Y₂)/(1 + mratio)` with the entropy of mixing recomputed
  for the merged composition. Pure, zero-allocation, smooth in `mratio`.
  `mixed(sys, 0)` is stream 1; stream 2 is the `mratio → ∞` limit. Does
  *not* do the energy balance (outlet temperature) of legacy `gas_mixing` —
  that is the caller's job (`temperature(gas, h = ...)` on the
  mass-weighted enthalpy).
- **humid air** — dry air (`Xair`) plus water vapor:
  `humid_air(; SH, RH, T, P) -> FrozenGas`, a constructor (not a hot path)
  taking either the specific humidity ω [kg water/kg dry air] or relative
  humidity converted via the legacy August–Roche–Magnus
  `saturation_vapor_pressure`. Same composition logic as the legacy
  `generate_humid_air` (water at `ω/ε` moles per mole dry air, renormalized).
- **temperature (the inversion verb)** — solving a property relation
  backwards for temperature; the keyword names what is known:
  `temperature(gas, h = ...)` (given enthalpy) and
  `temperature(gas, T1 = ..., PR = ...; ηp)` (along an isentrope with
  optional polytropic efficiency). One verb for every gas flavor — the
  *type* selects the algorithm and tier, never the function name.
  (Internal positional engines `T_of_h`/`T_isentropic` are unexported.)
- **FastFrozenGas{mode}** — a `FrozenGas` plus two precomputed cubic-Hermite
  *inverse* tables (h → T and s0 → T): `FastFrozenGas(gas; mode, N, Tmin,
  Tmax)`. Accelerates only the inversions; the forward functions forward to
  the wrapped gas unchanged. Modes: `:seeded` (default) uses the table as a
  Newton *seed* (exact answers, same convergence contract as `FrozenGas`;
  out-of-range targets fall back to the cold-start solve); `:fast` is pure
  table lookup (|ΔT/T| ≲ 2e-9 at N = 256; `DomainError` out of range —
  never silent extrapolation).

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
