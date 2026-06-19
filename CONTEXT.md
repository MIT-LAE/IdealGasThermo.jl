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
  `(gas, T)` — no state, no globals on the hot path, generic over `Real`. It is
  **self-describing**: it also carries its source mole-fraction vector
  `X::SVector{Nspecies}` (spdict order, Σ = 1), so a gas remembers what it is
  made of and can be re-mixed ([`mix`](@ref)) or re-burned ([`Vitiator`](@ref))
  — two lumped gases alone cannot, since the entropy of mixing `−Σ Xᵢ ln Xᵢ` is
  unrecoverable from the lumped coefficients (ADR-0008). `X` is an immutable
  constant of the composition, not state — it carries no T or P. All
  composition→`FrozenGas` construction routes through one zero-allocation kernel
  (`FrozenGas(X)`) over the shared species table `SPALOW`/`SPAHIGH`/… and the
  `_lump_molar` lumping function.
- **Entropy complement (ϕ / s0)** — φ(T) = ∫cp/T dT from the standard state;
  `s(T, P) = ϕ(T) − R·ln(P/Pstd)`. The temperature-only part of entropy.
- **Enthalpy datum** — enthalpies are **formation-inclusive** (CEA-style):
  `h(gas, 298.15 K)` equals the mixture's mass-specific formation enthalpy,
  not zero. Sensible enthalpy from 298.15 K is `h(gas, T) − h(gas, 298.15)`.
- **Vitiated mixture** — combustion products of a fuel + oxidizer at a given
  **FAR** (fuel–air mass ratio), frozen composition, complete combustion.
- **Vitiator** — a precomputed fuel + oxidizer combustion system, built once
  from a fuel and an oxidizer; stores **only the reaction stoichiometry**
  (`Xin`, `ΔX`, `massratio`, `ηburn`) — the NASA-9 lumping uses the shared
  module-const basis. The pure, allocation-free replacement for the Dict-based
  `vitiated_species` path on the hot path. Deliberately **not** named
  `Combustor`: that noun is reserved for the hardware component (pressure drop,
  efficiency, geometry) one abstraction level up in a cycle deck (ADR-0008).
  Construction is two methods: a `species` method does the work; the
  `AbstractString` method resolves the fuel name and forwards.
- **products** — `products(sys::Vitiator, FAR) -> FrozenGas`: the
  combustion-product gas at a given FAR. Pure, zero-allocation, smooth in
  FAR (ForwardDiff through FAR works; `FrozenGas{TF}` widens its eltype).
- **mix** — the merge of two gases, a free function (no precomputed system —
  each `FrozenGas` already carries its composition `X`); ADR-0008.
  - `mix(a::FrozenGas, b::FrozenGas, mratio) -> FrozenGas`: the merged
    *composition* at mass ratio `mratio = mass_b/mass_a`. The merged mole
    fractions are `X = (n_a·a.X + n_b·b.X)/(n_a + n_b)` with molar amounts
    `n_a = 1/a.MW`, `n_b = mratio/b.MW` (equivalent to the mass-fraction law of
    mixtures), and the merged gas is rebuilt from `X` with the entropy of mixing
    of the *merged* composition. Pure, zero-allocation, smooth in `mratio`.
    `mix(a, b, 0)` is `a`; `b` is the `mratio → ∞` limit.
  - `mix(a::GasState, b::GasState, mratio) -> GasState`: the above **plus the
    energy balance** — the mixed temperature is the mass-averaged total enthalpy
    `h = (h_a + mratio·h_b)/(1 + mratio)` inverted on the merged gas (the energy
    balance is *intrinsic* to mixing, not the caller's job). Requires equal
    stream pressures (an isobaric mixer); a non-isobaric mix needs a **momentum
    closure**, which — like mass-flow bookkeeping — lives in the flow layer
    (PowerCycles), not here (ADR-0005). Returns the merged stagnation state.
  - Replaces the `Mixer`/`mixed` pair (added then deleted during the `1.1.0`
    work — never released).
- **humid air** — dry air (`Xair`) plus water vapor:
  `humid_air(; SH, RH, T, P) -> FrozenGas`, a constructor (not a hot path)
  taking either the specific humidity ω [kg water/kg dry air] or relative
  humidity converted via the legacy August–Roche–Magnus
  `saturation_vapor_pressure`. Same composition logic as the legacy
  `generate_humid_air` (water at `ω/ε` moles per mole dry air, renormalized).
- **T_from_h (the inversion verb)** — solving a property relation backwards
  for temperature: `T_from_h(gas, hspec)` (the inverse of `h(gas, T)`). One
  verb for every gas flavor — the *type* selects the algorithm and tier
  (`FrozenGas` Newton, `FastFrozenGas{:seeded}`/`{:fast}`), never the function
  name. Named in the direction of the computation (`T_from_h`); an analogous
  `T_from_s0` would invert entropy. Replaces the former `temperature(gas; h)`
  keyword facade (removed, ADR-0004 update): keyword args don't participate in
  Julia dispatch, so a clean inversion is an explicitly-named function, not a
  `nothing`-checking facade. The isentrope form was never an inversion at all —
  a polytropic change of state is a *process*, expressed by `compress`/`expand`.
  The internal isentropic/polytropic temperature engine is **`_T_polytropic`**
  (unexported; renamed from the misleading `T_isentropic`, which claimed
  isentropic even with `ηp ≠ 1`).
- **process verbs** — the three-process taxonomy on the pure core
  (ADR-0004), each pure and allocation-free, each with the direction in
  the verb, never in the number:
  - *ratio-specified*: `compress(gas, T1, PR; ηp)` and
    `expand(gas, T1, PR; ηp)` — scalar kernels `T1 -> T2`, **both with
    PR ≥ 1** (`ArgumentError` otherwise); `expand` uses the expansion ηp
    convention `s0(T2) = s0(T1) + R·ηp·ln(1/PR)`, matching the legacy
    `expand(gas, 1/PR, ηp)`. State-layer methods on `GasState` update the
    pressure rail too; `expand_to(st, P2; ηp)` is the nozzle convenience
    (target pressure instead of ratio, requires P2 ≤ st.P). The loss model
    is **either** polytropic `ηp` (default) **or** isentropic `ηs`, never
    both (`ArgumentError`; ADR-0005): `ηs` degrades the ideal enthalpy
    change to the *same* PR (compressor `Δh = Δh_ideal/ηs`, turbine
    `Δh = ηs·Δh_ideal`), so the outlet pressure is identical to the `ηp`
    case and only the temperature differs.
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
- **speed of sound** — `speed_of_sound(gas, T) = √(γ·R·T)` [m/s]: a *pure
  property* of `(gas, T)` (alongside `cp`/`gamma`), needing no pressure;
  forwards through `FastFrozenGas` and has a `GasState` accessor
  `speed_of_sound(st)`. `mach(gas, T, V) = V/a` (and `mach(st, V)`) is the
  Mach number of a flow of speed `V`.
- **stagnation_state / static_state** — the gas-dynamics pair on `GasState`
  (ADR-0005), built from the enthalpy/entropy curves, **not** the constant-γ
  relations `1 + ½(γ−1)M²` (which ADR-0004 rejected the analog of).
  `stagnation_state(st, M)` brings the static flow (speed `V = M·a`)
  **isentropically** to rest: total enthalpy `h_t = h(st) + ½V²` (Tt by the
  h-inversion), entropy preserved (`Pt = P·exp((s0(Tt) − s0(T))/R)`).
  `static_state(st, M)` is the inverse — the static state at Mach `M` whose
  stagnation state is `st`, by a bounded Newton solve of
  `h(T) + ½(M·a(T))² = h(st)`. They are exact inverses to the inversion
  tolerance and reproduce the legacy `gas_Mach!(gas, 0, M, 1)`. A stagnation
  state is a **loss-free reference** — it carries no efficiency; a lossy ram
  is composed separately (recovery factor or `expand`). Named `…_state` (not
  bare `static`/`stagnation`) to avoid shadowing common identifiers and to say
  they return a `GasState`. Supersede the orphaned, never-`include`d
  `FlowStations.jl`.

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
- The mutable `Gas`/`Gas1D` types are the legacy *stateful convenience layer*.
  As of `1.1.0` the **whole layer is loudly deprecated** and on a
  scheduled-deletion path (ADR-0002, ADR-0007): `Gas{N}`, `Gas1D`, and the
  Dict-combustion / mutable-turbo functions they back are removed in a future
  `2.0.0`, after a `1.x` migration window. The pure core no longer touches
  them — `DryAir` is now built directly from `Xair` (`generate_composite_species(
  Xidict2Array(Xair), …)`), not through `Gas()` — so the deletion is purely
  subtractive. **Migration:** `Gas{N}` → `FrozenGas` (composition + properties);
  a mutable `Gas`/`Gas1D` parcel → a `GasState` (gas, T, P) point with the process
  verbs; `set_TP!`/`set_h!`/`set_hP!`/`set_Δh!` → `GasState` constructors +
  `compress`/`expand`/`add_heat`/`add_work`. Deprecation is signalled by a loud,
  once-per-session `@warn` on the `Gas`/`Gas1D` **constructors** only (the choke
  point every legacy path goes through; `src/deprecation.jl`, ADR-0007).
- `FrozenGas` keeps its name permanently — it is not renamed to `Gas` in
  v2.0 (ADR-0002: no name recycling across a semantics flip; "frozen" names
  the no-dissociation contract).
- **Exported property API.** The pure-core accessors are exported for
  unqualified use by consumers (e.g. PowerCycles): `cₚ`/`c_p`, `h`, `s0`,
  `gamma`/`γ`, `R`, `T_from_h`, `pressure_ratio` (plus `props`,
  `entropy`, `density`, `speed_of_sound`). The isentropic/polytropic engine
  `_T_polytropic` is **un**exported (the public process API is `compress`/
  `expand`), and the old `temperature` keyword facade is gone. Specific heat is
  exported as **`cₚ` and `c_p`** — interchangeable aliases (`const cₚ = cp`,
  `const c_p = cp`) of the internal function `cp`, which is *not* exported
  because the bare name collides with `Base.cp` (file copy). `cp` stays usable
  as `IdealGasThermo.cp`, and the `props` NamedTuple keeps the field name `cp`
  (field accessors are reached only as `.cp`, never as a bare identifier, so no
  collision). Internal code and tests continue to use `cp`.
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
  The inversions (`T_from_h`, `_T_polytropic`) need the **full
  three-term IFT rule**: the constant-substance rules account only for the
  *target* moving and silently drop the *composition moves* term, which when
  the gas is Dual-typed produces a nested-Dual result instead of a number. The
  extension's substance-Dual rules dispatch on `FrozenGas{<:Dual}` and add that
  term —
  `∂T = (partials(h_spec) − partials(h(gas, T*))) / cp(gas₀, T*)` for
  `T_from_h` — while keeping the Newton loop on the value rail (strip all
  tangents, solve once in `Float64`, attach the closed-form tangent at `T*`
  via one forward evaluation). This preserves the split-rule speed: the
  substance-Dual inversion is ~36× faster than differentiating through the
  loop and zero-allocation, within ~13% of the constant-substance baseline.

## Testing vocabulary (see docs/adr/0006)

- **property (metamorphic) test** — a check that holds over the whole space of
  gases the package can build, independent of any reference implementation:
  `∫cp` vs `h`/`s0` by quadrature, second-law entropy generation, process
  round-trips, `s(T,kP) = s(T,P) − R·ln k`, mixing symmetry, combustion mass
  conservation, AD-vs-FD. The permanent oracle (`test/unit_test_properties.jl`).
- **reference anchor** — an absolute-value check against an *external
  authority*, not an in-package twin. The CEA anchors
  (`test/unit_test_cea_reference.jl`, data in `test/CEA_output.txt`) evaluate
  the *same* NASA-9 coefficients the package loads, so `cp`/`s0`/`h` are
  compared in CEA's printed molar units to its printed precision
  (`atol = 1e-3` = 1 unit in the last place: 5e-4 CEA rounding + the small
  `Runiv = 8.3145` vs CEA `8.31451` difference, heaviest on `s0`). Anchors the
  one thing metamorphic tests cannot — a systematic datum/units/scale error.
- **migration test** — a temporary "new core agrees with the legacy mutable
  layer" check (`≈ Gas1D`, `≈ vitiated_species`, `≈ set_Δh!`). Scaffolding
  only: it is a same-NASA-9-kernel twin (circular) and dies with the mutable
  layer at v2.0. The pure core is **not** validated this way (ADR-0006);
  `Gas`/`Gas1D` keep a minimal smoke test in the legacy-only files until v2.0.
