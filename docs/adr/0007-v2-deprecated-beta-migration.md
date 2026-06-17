# ADR-0007: Retire the mutable legacy layer via a loudly-deprecated `2.0.0-beta` series, then delete at `2.0.0` final

Date: 2026-06-17
Status: accepted

## Context

ADR-0001/0002 established the immutable pure core (`FrozenGas`, `GasState`,
`Combustor`/`products`, `Mixer`/`mixed`, the flow verbs) as *the* architecture and
declared the mutable legacy layer — `Gas{N}`, `Gas1D`, the Dict-combustion in
`combustion.jl`, the mutable-turbo in `turbo.jl`, `thermoProps.jl`, and the
`Gas`-based `print_thermo_table` — deprecated, to be removed at v2.0. ADR-0006 then
re-grounded the test suite so the pure core is validated by property tests + CEA
anchors, **not** by agreement with the legacy types. With that done, the only thing
the live (non-legacy) code still needed from the legacy layer was a single tie:
`DryAir` was built through `Gas()`.

Two facts shape *how* we remove the layer:

1. **A live downstream consumer (PowerCycles.jl) must migrate against a real,
   Pkg-installable release.** A hard delete in one breaking commit gives consumers
   nothing to migrate *against* — they would jump from "legacy present" to "legacy
   gone" with no overlap window. The package is not in the General registry, so a
   consumer installs it by `url = ` + `rev = ` (a tag or branch).
2. **Removing exported API (`Gas`, `Gas1D`) is what *defines* the major version
   bump.** Under semver, the whole `1.x → 2.0` transition is one breaking step.

## Decision

Retire the legacy layer in **two phases**, separated by a migration window:

1. **Beta phase (this branch → `2.0.0-betaN` tags).** The legacy layer stays present
   but is made **loudly deprecated**, and the pure core is cut free of it:
   - **Untie `DryAir` from `Gas()`** (`IdealGasThermo.jl`): build it directly from the
     `Xair` mole-fraction table via `generate_composite_species(Xidict2Array(Xair),
     "Dry Air")`. This is identical to the old `let g = Gas(); g.X = Xair; … end`
     form to machine precision (MW bit-equal; cp/h/s0 agree to ~3e-16). After this,
     **no non-legacy code path can reach `Gas`/`Gas1D`** — the deletion test passes
     and the eventual deletion is purely subtractive.
   - **Loud, once-per-use deprecation warnings** (`src/deprecation.jl`,
     `_legacy_warn`). The warning is a plain `@warn` — *always shown*, unlike
     `Base.depwarn`, which Julia hides unless `--depwarn=yes` — bounded to one
     emission per entry point per session (`maxlog = 1`, keyed by `_id`). The warns
     live **only on the `Gas` and `Gas1D` constructors**: every legacy public
     function (`set_TP!`, `set_h!`, `set_hP!`, `set_Δh!`, `print_thermo_table`,
     `gas_burn`, `vitiated_species`, the mutable-turbo verbs) takes a `Gas`/`Gas1D`,
     so constructing one is the single choke point that signals a caller is on the
     legacy path. Putting a warn on `set_TP!`/`set_hP!` would cascade through internal
     plumbing (they are called all over `combustion.jl`/`turbo.jl`/`setproperty!`).
     Those downstream symbols instead carry docstring deprecation notes.
   - **Version → `2.0.0-beta1`** (a semver pre-release of the same 2.0.0 that will
     delete the layer).
2. **Deletion phase (later branch → `2.0.0` final).** Delete `Gas.jl`, `Gas1D.jl`, the
   legacy halves of `combustion.jl`/`turbo.jl`/`io.jl`, `thermoProps.jl`, and the
   legacy-only smoke test files (`unit_test_turbo`, `unit_test_vitiated`,
   `unit_test_composite`, the legacy parts of `unit_test_mixthermo`). Replace the
   remaining `vitiated_species(...)` *fixtures* in `unit_test_mixing`/
   `unit_test_properties` with `products(Combustor(...), FAR)` (ADR-0006 §75). This
   touches no pure-core oracle. Because the betas are pre-releases of `2.0.0`, `Gas`
   only ever exists in pre-release versions; the first *stable* 2.x release has it gone.
   **Before** deleting the legacy test files, two coverage items must be re-homed into
   pure-core test files, because they cover functions reachable from the pure core
   (verified by call-graph trace): (a) the coefficient-level correctness check of
   `generate_composite_species` (composite vs the fitted "Air" species: MW / Hf /
   `alow` / `ahigh`), currently only in `test/unit_test_composite.jl` —
   `generate_composite_species` is reachable from `DryAir`, `FrozenGas(X, name)`, and
   `humid_air`; (b) the component-value pin of `reaction_change_molar_fraction`,
   currently only in `test/unit_test_vitiated.jl` — it is called directly by the pure
   `Combustor` constructor (`src/combustor.jl`). The other two legacy combustion
   functions (`stoich_molar_FOR`, `vitiated_mixture`) are legacy-only and their tests
   die correctly with the layer.

PowerCycles.jl pins to a published `2.0.0-betaN` tag on
`github.com/MIT-LAE/IdealGasThermo.jl`
(`Pkg.add(url = "https://github.com/MIT-LAE/IdealGasThermo.jl", rev = "v2.0.0-beta1")`;
`using IdealGasThermo`), migrates off `Gas`/`Gas1D` while the deprecation warns guide
it, then moves to `2.0.0` final.

## Consequences

- The beta is **not itself a breaking removal** — it adds warnings and a version
  pre-release tag; no exported symbol disappears. The breaking event is the `2.0.0`
  final deletion, and the single `1.x → 2.0` major bump covers it (semver-clean).
- Because the loud warns are `@warn`, **not** `Base.depwarn`, they are shown
  regardless of the `--depwarn` flag and are *not* escalated to errors by
  `--depwarn=error`. The one constraint on the test suite: do **not** wrap legacy
  construction in `@test_nowarn`/`@test_logs`-without-`:warn`; the legacy smoke tests
  emit one `Gas` and one `Gas1D` warning (bounded by `maxlog = 1`) and stay green.
- Do not add runtime deprecation warnings to the downstream legacy functions; the
  constructor choke point already covers every legacy path, and per-verb warns would
  cascade through internal plumbing.
- This does not re-litigate ADR-0002 (the layer *is* going away); it records the
  *mechanism and sequencing*, and the reason the deletion was split from the
  deprecation: a real consumer needs an installable migration window.
