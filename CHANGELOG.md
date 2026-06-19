# Changelog

All notable changes to IdealGasThermo.jl are documented here. The format loosely
follows [Keep a Changelog](https://keepachangelog.com/); the project follows
[Semantic Versioning](https://semver.org/).

## [1.1.0] — 2026-06-19

A purely **additive** minor release over `1.0.0`: the immutable pure core becomes
the headline API, and the legacy mutable layer is **loudly deprecated but fully
working**. Nothing exported is removed or renamed — `1.1.0` is backward-compatible
with `1.0.0` and **safe to adopt and migrate against**; downstream packages can pin
this tag (or its SHA). The breaking removal of the legacy layer is deferred to a
future **`2.0.0`** (ADR-0002, ADR-0007).

### Added

The entire immutable pure core is **net-new** in this release (it did not exist in
`1.0.0`). `using IdealGasThermo` now makes the following callable without
qualification:

- **`FrozenGas`** — the immutable, `isbits` substance: mass-specific NASA-9 property
  curves that are pure functions of `(gas, T)` (no state, no globals on the hot
  path, generic over `Real`). It is **self-describing** — it carries its source
  mole-fraction vector `X`, so a gas remembers what it is made of and can be
  re-mixed or re-burned (ADR-0008).
- **Property accessors:** `cₚ` / `c_p` (specific heat; interchangeable aliases),
  `h`, `s0`, `gamma` / `γ` (ratio of specific heats), `R`, and the `props(gas, T)`
  NamedTuple. The bare name `cp` is **not** exported because it would shadow
  `Base.cp` (file copy); it stays reachable as `IdealGasThermo.cp`, and the `props`
  NamedTuple keeps the field name `cp` (a field accessor cannot collide).
- **Inversion:** `T_from_h(gas, hspec)` — the inverse of `h(gas, T)`
  (enthalpy → temperature) — and `pressure_ratio(gas, T1, T2)`.
- **`GasState`** — an immutable `(substance, T, P)` value record (ergonomics, not a
  parcel): it keeps the caller's T-rail and P-rail travelling together through a
  process chain. Accessors `cp` / `h` / `s0` / `gamma` / `R` at `st.T`, plus
  `entropy(st)` and `density(st)`.
- **Process verbs** (each pure, allocation-free, direction in the verb):
  `compress`, `expand`, `expand_to` (nozzle, target-pressure form), `add_heat`,
  `add_work`, `extract_work`. The loss model is **either** polytropic `ηp` (default)
  **or** isentropic `ηs`, never both (ADR-0004 / ADR-0005).
- **Flow / gas-dynamics verbs** (ADR-0005): `speed_of_sound`, `mach`,
  `stagnation_state`, `static_state` — built from the enthalpy / entropy curves,
  never the constant-γ relations.
- **Combustion:** `Vitiator` (the *composition* model of combustion — a precomputed
  fuel + oxidizer system) and `products` —
  `products(Vitiator(fuel, oxidizer), FAR) -> FrozenGas`. Pure, zero-allocation,
  smooth in FAR (ForwardDiff through FAR widens the gas's eltype).
- **Mixing:** `mix`, a **free function** — `mix(a, b, mratio)` for two `FrozenGas`
  (composition blend) or two `GasState` (composition blend **plus** the
  mass-averaged enthalpy energy balance). Each gas already carries its `X`, so no
  precomputed mixing object is needed (ADR-0008).
- **`FastFrozenGas`** — a `FrozenGas` plus precomputed cubic-Hermite inverse tables
  (h → T and s0 → T) that accelerate only the inversions; the forward functions are
  unchanged. Modes `:seeded` (exact; table as Newton seed) and `:fast` (pure table
  lookup) (ADR-0003).
- **`humid_air(; SH, RH, T, P) -> FrozenGas`** — dry air plus water vapor, by
  specific or relative humidity.
- **ForwardDiff support** via a package extension: analytic, closed-form derivatives
  (`dh/dT = cp`, `ds0/dT = cp/T`) and implicit-function-theorem rules for the
  inversions — including the substance-`Dual` rule for a `FrozenGas{<:Dual}` produced
  by `products(sys, FAR::Dual)` — never by differentiating a Newton loop.

### Deprecated

- **`Gas{N}`, `Gas1D`, and the functions they back** (the Dict-combustion in
  `combustion.jl`, the mutable-turbo in `turbo.jl`, `thermoProps.jl`, and the
  `Gas`-based `print_thermo_table`). They still work; constructing a `Gas` or
  `Gas1D` now emits a **loud, once-per-session** deprecation warning. Scheduled for
  removal in a future `2.0.0` (ADR-0002, ADR-0007).

  **Migration to the pure core:**

  | Legacy (deprecated)                                 | Pure core (use instead)                                              |
  | --------------------------------------------------- | -------------------------------------------------------------------- |
  | `Gas()` / `Gas(Y)` (composition + properties)       | `FrozenGas(...)` — immutable substance; `c_p(gas, T)`, `h(gas, T)`, … |
  | `Gas1D(sp)`                                          | `FrozenGas(sp)`                                                       |
  | a mutable `Gas` / `Gas1D` parcel (state in the gas)  | `GasState(gas, T, P)` — an explicit (gas, T, P) point                |
  | `set_TP!(gas, T, P)`                                | `GasState(gas, T, P)`                                                |
  | `set_h!` / `set_hP!` / `set_Δh!`                    | `compress` / `expand` / `add_heat` / `add_work` / `extract_work`     |
  | enthalpy → temperature inversion                    | `T_from_h(gas, hspec)`                                               |
  | `gas_burn` / `vitiated_species` (Dict-combustion)   | `products(Vitiator(fuel, oxidizer), FAR)`                           |
  | mutable-turbo mixing                                | `mix(a, b, mratio)` (free function; `FrozenGas` or `GasState`)       |
  | `gas_Mach!`                                         | `static_state` / `stagnation_state`                                  |

### Changed

- **`DryAir` is now built directly from the `Xair` mole-fraction table**
  (`generate_composite_species(Xidict2Array(Xair), "Dry Air")`) instead of through
  `Gas()`. This cuts the last tie between the pure core and the legacy layer, so the
  layer can be removed in the future `2.0.0` without touching any live path. The result is
  identical to the previous construction to machine precision (molecular weight
  bit-equal; `cp` / `h` / `s0` agree to ~1e-14). No user-facing change.

### Removed (internal)

- The **unexported** mass-fraction combustion helper
  `IdealGasThermo.reaction_change_fraction` was deleted, and its docs reference was
  removed in lockstep (no dangling docs). Its molar sibling,
  `reaction_change_molar_fraction`, is preserved (re-homed into the `Vitiator`
  module). Internal only — no exported API is affected.

### Notes

- **Nothing exported was removed or renamed.** Every symbol the prior release
  exported is still exported (`Gas`, `set_h!`, `set_hP!`, `set_TP!`, `set_Δh!`,
  `AbstractSpecies`, `species`, `composite_species`, `generate_composite_species`,
  `readThermo`, `species_in_spdict`, `Gas1D`, `print_thermo_table`, `X2Y`, `Y2X`,
  `DryAir`). Exports are **purely additive**.
- The deprecation warnings are plain `@warn` (always shown), bounded to one emission
  per type per session. They are **not** `Base.depwarn`, so they are unaffected by
  the `--depwarn` flag (neither hidden by its default nor escalated to errors by
  `--depwarn=error`).
- The `1.0.0` legacy test suite (the `Gas` / `Gas1D` suite, 514 tests) runs unchanged
  against this source and passes — no hard breaks, no numeric drift.

## [1.0.0]

The mutable, stateful `Gas` / `Gas1D` convenience layer for NASA-9 thermodynamics —
the baseline this `1.1.0` extends. The immutable pure core did not yet exist. This is
the **sole version registered in the Julia General registry**, so `Pkg.add("IdealGasThermo")`
installs it today. (General records `1.0.0` against a commit; the matching `v1.0.0`
git tag was never pushed — TagBot did not run — so the repository shows only the
earlier, hand-made `v0.1`, `v0.1.1`, `v0.1.2` tags, which were never registered.)

[1.1.0]: https://github.com/MIT-LAE/IdealGasThermo.jl/releases/tag/v1.1.0
[1.0.0]: https://github.com/MIT-LAE/IdealGasThermo.jl
