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
  backwards for temperature: `temperature(gas, h = ...)` (given enthalpy).
  One verb for every gas flavor — the *type* selects the algorithm and
  tier, never the function name. (Internal positional engines
  `T_of_h`/`T_isentropic` are unexported.) The former isentrope form
  (`T1 = ..., PR = ...; ηp`) is removed (ADR-0004): a polytropic change of
  state is a *process*, not an inversion — use the process verbs
  `compress`/`expand`.
- **process verbs** — the three-process taxonomy on the pure core
  (ADR-0004), each pure and allocation-free, each with the direction in
  the verb, never in the number:
  - *ratio-specified*: `compress(gas, T1, PR; ηp)` and
    `expand(gas, T1, PR; ηp)` — scalar kernels `T1 -> T2`, **both with
    PR ≥ 1** (`ArgumentError` otherwise); `expand` uses the expansion ηp
    convention `s0(T2) = s0(T1) + R·ηp·ln(1/PR)`, matching the legacy
    `expand(gas, 1/PR, ηp)`. State-layer methods on `GasState` update the
    pressure rail too; `expand_to(st, P2; ηp)` is the nozzle convenience
    (target pressure instead of ratio, requires P2 ≤ st.P).
  - *work-specified*: `add_work(st, w; ηp)` / `extract_work(st, w; ηp)`,
    `w ≥ 0` [J/kg]; enthalpy ±w with the pressure on the polytrope
    `P2 = P1·exp(K/R·Δs0)`, K = ηp adding / 1/ηp extracting (the legacy
    `set_Δh!` conventions, owned by the verbs).
  - *heat at constant pressure*: `add_heat(st, q)`, signed `q` [J/kg].
- **GasState** — `GasState(gas, T, P)`: an immutable (substance, T, P)
  *value* record — ergonomics, not architecture (ADR-0004). The substance
  stays a pure set of property curves; the record only makes the caller's
  (T, P) pair travel together through a process chain so the T-rail and
  P-rail cannot diverge. `isbits` for `FrozenGas{Float64}`; never mutated —
  every process verb returns a NEW state. Read-only accessor functions
  (no getproperty magic): `cp/h/s0/gamma/R` at `st.T`, plus
  `entropy(st) = s0(T) − R·ln(P/Pstd)` and `density(st) = P/(R·T)`
  (exported full words; `s`/`rho` are unexported aliases). Stores no
  derived properties — that would be the caching ADR-0001 forbids.
- **FastFrozenGas{mode}** — a `FrozenGas` plus two precomputed cubic-Hermite
  *inverse* tables (h → T and s0 → T): `FastFrozenGas(gas; mode, N, Tmin,
  Tmax)`. Accelerates only the inversions; the forward functions forward to
  the wrapped gas unchanged. Modes: `:seeded` (default) uses the table as a
  Newton *seed* (exact answers, same convergence contract as `FrozenGas`;
  out-of-range targets fall back to the cold-start solve); `:fast` is pure
  table lookup (|ΔT/T| ≲ 2e-9 at N = 256; `DomainError` out of range —
  never silent extrapolation).

## Architecture terms (see docs/adr/)

- **Substance vs state**: a `FrozenGas` is a *set of property curves*, not a
  parcel — it holds only constants of the composition (coefficients, MW, R,
  Hf), never T or P. Temperature is an **argument, not an attribute**:
  `h(gas, T)` reads as h_gas(T). The only thermodynamic state in the system
  lives with the caller (in a cycle solver: the solver's own unknown
  vector). Corollary: caching immutable facts about the *curves* (e.g.
  `FastFrozenGas` inverse tables, `gas.R`) is fine; caching facts about
  "the current state" (the old `Tarray`/`gas.cp` pattern) is what the
  architecture forbids.
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
- **Dual-carrying gas** — a `FrozenGas{<:Dual}`, i.e. a substance whose
  *coefficients themselves* carry a tangent, as produced by
  `products(sys, FAR::Dual)` (the product composition depends on FAR, so the
  mass-scaled NASA-9 coefficients, MW, R, and Hf all carry the FAR-derivative).
  The parametric eltype `FrozenGas{TF<:Real}` is what makes this legal: a
  Dual-valued argument simply widens the gas. Forward property reads through a
  Dual-carrying gas are intrinsically cheap because every property is **linear
  in the coefficients** (the lone nonlinearity, `log T`, rides the temperature
  rail), so the tangent propagates by scale-and-add with no transcendentals.
  The inversions (`temperature`/`T_of_h`, `T_isentropic`) need the **full
  three-term IFT rule**: the constant-substance rules account only for the
  *target* moving and silently drop the *composition moves* term, which when
  the gas is Dual-typed produces a nested-Dual result instead of a number. The
  extension's substance-Dual rules dispatch on `FrozenGas{<:Dual}` and add that
  term —
  `∂T = (partials(h_spec) − partials(h(gas, T*))) / cp(gas₀, T*)` for
  `T_of_h` — while keeping the Newton loop on the value rail (strip all
  tangents, solve once in `Float64`, attach the closed-form tangent at `T*`
  via one forward evaluation). This preserves the split-rule speed: the
  substance-Dual inversion is ~36× faster than differentiating through the
  loop and zero-allocation, within ~13% of the constant-substance baseline.
