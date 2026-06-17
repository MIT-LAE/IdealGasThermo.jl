# Changelog

All notable changes to IdealGasThermo.jl are documented here. The format loosely
follows [Keep a Changelog](https://keepachangelog.com/); the project follows
[Semantic Versioning](https://semver.org/).

## [2.0.0-beta1] — 2026-06-17

Pre-release of `2.0.0`. The immutable pure core is now the headline API and the
legacy mutable layer is **loudly deprecated** ahead of its removal in `2.0.0` final.
Nothing is removed yet — this release is safe to adopt and migrate against.

### Added

- **The pure-core property accessors are now exported**, so `using IdealGasThermo`
  makes them callable without qualification: `cₚ` / `c_p` (specific heat), `h`,
  `s0`, `gamma` / `γ`, `R`, `T_of_h`, `T_isentropic`, `pressure_ratio`. Previously
  these were internal (`IdealGasThermo.cp(gas, T)`). The ratio of specific heats is
  exported under both `gamma` and the Unicode alias `γ`.
  - Specific heat is exported as **`cₚ` and `c_p`** (interchangeable aliases of one
    function). The bare name `cp` is *not* exported because it would shadow
    `Base.cp` (file copy); it remains reachable as `IdealGasThermo.cp`. The
    `props(gas, T)` NamedTuple still uses the field name `cp` (a field accessor is
    reached only as `.cp`, so it cannot collide with `Base.cp`).

### Deprecated

- **`Gas{N}`, `Gas1D`, and the functions they back** (the Dict-combustion in
  `combustion.jl`, the mutable-turbo in `turbo.jl`, `thermoProps.jl`, and the
  `Gas`-based `print_thermo_table`). Constructing a `Gas` or `Gas1D` now emits a
  loud, once-per-session deprecation warning. These are scheduled for deletion in
  `2.0.0` final (ADR-0002, ADR-0007).

  **Migration to the pure core:**

  | Legacy (deprecated)                              | Pure core (use instead)                                             |
  | ------------------------------------------------ | ------------------------------------------------------------------- |
  | `Gas()` / `Gas(Y)` (composition + properties)    | `FrozenGas(...)` — immutable substance; `cp(gas, T)`, `h(gas, T)`, … |
  | `Gas1D(sp)`                                       | `FrozenGas(sp)`                                                     |
  | a mutable `Gas`/`Gas1D` parcel (state in the gas) | `GasState(gas, T, P)` — an explicit (gas, T, P) point               |
  | `set_TP!(gas, T, P)`                              | `GasState(gas, T, P)`                                              |
  | `set_h!` / `set_hP!` / `set_Δh!`                  | `compress` / `expand` / `add_heat` / `add_work` / `extract_work`    |
  | `gas_burn` / `vitiated_species` (Dict-combustion) | `products(Combustor(fuel, oxidizer), FAR)`                          |
  | mutable-turbo mixing                              | `mixed(Mixer(...), ...)`                                            |
  | `gas_Mach!`                                       | `static_state` / `stagnation_state`                                 |

### Changed

- **`DryAir` is now built directly from the `Xair` mole-fraction table**
  (`generate_composite_species(Xidict2Array(Xair), "Dry Air")`) instead of through
  `Gas()`. This cuts the last tie between the pure core and the legacy layer; the
  result is identical to the previous construction to machine precision (molecular
  weight is bit-equal; `cp`/`h`/`s0` agree to ~3e-16). No user-facing change.

### Notes

- The deprecation warnings are plain `@warn` (always shown), bounded to one emission
  per type per session. They are **not** `Base.depwarn`, so they are unaffected by the
  `--depwarn` flag (neither hidden by its default nor escalated to errors by
  `--depwarn=error`).

## [1.0.0]

First stable release: immutable pure core (`FrozenGas`, `FastFrozenGas`, `GasState`,
`Combustor`/`products`, `Mixer`/`mixed`, `humid_air`, gas-dynamics flow verbs) shipped
alongside the (now deprecated) mutable `Gas`/`Gas1D` layer.

[2.0.0-beta1]: https://github.com/MIT-LAE/IdealGasThermo.jl/releases/tag/v2.0.0-beta1
[1.0.0]: https://github.com/MIT-LAE/IdealGasThermo.jl
