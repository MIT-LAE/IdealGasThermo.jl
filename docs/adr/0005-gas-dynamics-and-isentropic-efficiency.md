# ADR-0005: Gas dynamics on the pure core (speed of sound, Mach, stagnation/static) and the isentropic-efficiency form of the ratio verbs

Date: 2026-06-15
Status: accepted (extends ADR-0004's process layer; supersedes the orphaned
`FlowStations.jl` and the legacy `gas_Mach!`)

## Context

Two capabilities of the legacy stateful layer had no equivalent on the pure
core, blocking feature parity for compressible-flow cycle code:

- **Mach-number state-setting.** The legacy `gas_Mach!(gas, M0, M, ηp)`
  (`src/turbo.jl`) advanced a mutable gas across a change in Mach number by an
  energy balance; the stagnation↔static helpers `isenTR`/`isenPR`/`get_static`
  lived in `src/FlowStations.jl` — a file that **was never `include`d in the
  module** (orphaned dead code) and which used the **constant-γ** relations
  `1 + ½(γ−1)M²` and `(…)^(γ/(γ−1))`. ADR-0004 already rejected the analog of
  those shortcuts (the `temp_ratio(gas, PR)` proposal) because for a
  thermally-perfect gas the ratio drifts ~10% with inlet temperature at fixed
  ratio. A faithful port therefore cannot copy `isenTR`/`isenPR`.

- **Isentropic efficiency.** Every legacy verb (`PressureRatio`,
  `gas_Mach!`, `set_Δh!`) and every ADR-0004 process verb is **polytropic-`ηp`
  only**. Turbomachinery is just as often specified by *isentropic* (adiabatic)
  efficiency `ηs`. This is net-new, not a regression.

## Decision

1. **Speed of sound is a property of `(gas, T)`**, not of a state:
   `speed_of_sound(gas, T) = √(γ·R·T)` lives at the `FrozenGas` level
   (alongside `gamma`), forwards through `FastFrozenGas`, and has a
   one-argument `GasState` accessor `speed_of_sound(st)`. Composition and
   temperature are all it needs — no pressure, consistent with "temperature is
   an argument, not an attribute" (ADR-0001). `mach(gas, T, V) = V/a` and
   `mach(st, V)` likewise.

2. **Stagnation is an isentropic reference state, computed from the enthalpy
   and entropy curves — never constant-γ.** `stagnation_state(st, M)` brings
   the static flow (speed `V = M·a`) to rest at constant total enthalpy
   `h_t = h(st) + ½V²` (`T_t` by the enthalpy inversion) and constant entropy
   (`P_t = P·exp((s0(T_t) − s0(T))/R)`). `static_state(st, M)` is the inverse:
   the static state at Mach `M` whose isentropic stagnation state is `st`,
   solving `h(T) + ½(M·a(T))² = h(st)` by a bounded Newton iteration (the
   residual is strictly increasing in `T`; the constant-γ value seeds it).
   `static_state(stagnation_state(st, M), M) == st` to the inversion
   tolerance, verified to M = 3. `static_state(st, 0)` and
   `stagnation_state(st, 0)` are the identity. The pair reproduces the legacy
   `gas_Mach!(gas, 0, M, 1)` to ~1e-8.

3. **Stagnation carries no efficiency.** A stagnation state is a *reference*
   (isentropic by definition), unlike the legacy `gas_Mach!`'s `ηp` argument.
   A lossy ram/diffuser is the caller's job — compose with a pressure-recovery
   factor or an `expand`. Keeping the reference loss-free is what makes
   `stagnation_state`/`static_state` exact inverses.

4. **Isentropic efficiency is a keyword on the ratio verbs, mutually
   exclusive with `ηp`.** `compress`/`expand`/`expand_to` accept *either*
   `ηp` (polytropic, the ADR-0004 default — entropy distributed along the
   path) *or* `ηs` (isentropic — the ideal enthalpy change to the **same
   pressure ratio**, degraded at the outlet): both efficiencies given is an
   `ArgumentError`. Compressor: `h(T2) = h(T1) + (h(T2s) − h(T1))/ηs`;
   turbine: `h(T2) = h(T1) − ηs·(h(T1) − h(T2s))`, with `T2s` the loss-free
   outlet. Both conventions land on the **same outlet pressure** (`P·PR` /
   `P/PR`) — only the temperature differs — so the state-layer pressure rail
   is unchanged. `ηs = 1` reproduces the isentrope. The conventions are
   reconcilable: measuring the isentropic efficiency of a polytropic outlet
   and feeding it back as `ηs` reproduces that outlet exactly (tested both
   directions).

5. **ForwardDiff.** `speed_of_sound`, `mach`, `stagnation_state`, and the
   `ηs` verbs are built on the property curves and the IFT-ruled engines
   (`T_of_h`, `T_isentropic`), so `Dual`s propagate analytically and
   allocation-free — no new extension methods. `static_state`'s Newton runs
   on the value rail; differentiating it is supported but is not a closed-form
   IFT rule (no AD-through-Mach contract is claimed, mirroring the "don't
   differentiate the loop" stance — a substance-`static_state` rule can be
   added later if a caller needs it).

6. **Names: `_state`-suffixed full words.** `speed_of_sound`, `mach`,
   `stagnation_state`, `static_state` are exported. The pair is named
   `…_state` rather than the bare `stagnation`/`static` for two reasons: a
   bare `static` is a common identifier that would shadow (or be shadowed by)
   caller bindings, and the suffix says what the functions *return* — a
   `GasState`, echoing the `GasState`/`entropy`/`density` state vocabulary.
   `speed_of_sound`/`mach` stay bare (they return scalars, and neither
   collides).

## Consequences

- Compressible-flow cycle code expresses a station as
  `stagnation_state`/`static_state`/`mach` over `GasState`, with the same
  exactness and zero-allocation guarantees as the rest of the pure core. No
  constant-γ approximation enters anywhere.
- `src/FlowStations.jl` (orphaned, constant-γ) is **deleted** (2026-06-16):
  its `isenTR`/`isenPR` were the rejected constant-γ approximation and its
  `get_static` is subsumed by `static_state` — nothing was lost. The one
  concept it held that the pure core does *not* have — a **flow station that
  carries geometry and a flow rate** (cross-sectional area, velocity, mass
  flux ρ·A·V) — is intentionally out of scope: a `GasState` is `(gas, T, P)`
  only. That flow-station record belongs upstream in **PowerCycles.jl**; do
  not re-add an area/velocity-carrying state to this package. The legacy
  mutable `gas_Mach!` remains only as legacy until the v2.0 removal of the
  mutable layer.
- The ratio verbs now answer both efficiency questions a turbomachinery
  designer asks, without a second verb name — the convention lives in the
  keyword, the direction still lives in the verb (ADR-0004).
- Do not re-introduce constant-γ stagnation ratios on the pure core, and do
  not give `stagnation_state`/`static_state` an efficiency argument: losses
  compose, the reference does not.
