# ADR-0004: Process verbs (compress/expand/add_heat/add_work/extract_work) and the GasState record

Date: 2026-06-12
Status: accepted (supersedes the isentrope-as-inversion form an earlier draft
hid inside a `temperature(gas; …)` verb)

> **Inversion naming (settled 2026-06-19).** An early draft routed the
> enthalpy → temperature inversion through a `temperature(gas; h = …)` keyword
> facade. That facade is **dropped**: Julia does not dispatch on keyword arguments,
> so a facade meant to grow across more inversions could only become a
> `nothing`-checking branch. The inversion is the explicitly-named, exported
> **`T_from_h(gas, hspec)`** (the inverse of `h(gas, T)`; an analogous `T_from_s0`
> would invert entropy). The internal isentropic/polytropic engine is
> **`_T_polytropic`** (unexported) — named for what it computes, since with `ηp ≠ 1`
> the outlet is not isentropic. The public process API is `compress`/`expand` (this
> ADR); §2 and §5 below are stated in these final names.

## Context

Cycle code written against the pure core had two recurring frictions,
studied and measured in `claude_sandbox/proto_ergo/` (single-spool turbojet
design point, three styles, identical numerics to 1.6e-14):

- **Process hiding inside an inversion verb.** An early isentrope form
  `temperature(gas, T1 = ..., PR = ...; ηp)` answered "what T comes out of
  this polytropic process" — but that verb was meant to *invert a property
  relation*. The `ηp` keyword was the tell: an inversion has no
  efficiency. A polytropic change of state is a **process**, and a process
  has a direction; hiding the direction in whether `PR` is above or below 1
  produced the legacy `compress`/`expand` pair on the mutable layer and
  sign-convention bugs around it.
- **Two parallel rails of scalars.** Functional cycle code carries a T-rail
  and a P-rail per station, advanced in lockstep by hand. The prototype's
  deliberately seeded one-character bug — nozzle expansion ratio computed
  from `P4` instead of `P5` — ran silently and shifted specific thrust by
  **16.5%**. The state-record style cannot express that bug: the pressure
  travels inside the value the verb returns.

Prototype measurements (`proto_ergo/study_out.log`, `bench_out.log`):

- GasState chain vs bare-scalar functional chain: **1608 vs 1588 ns**
  median (≈1.3% overhead), both 0 allocations; legacy mutable Gas1D chain:
  ~30 µs, 957 allocations.
- Line count for the same cycle: **11 lines (GasState) vs 15 (functional)
  vs 32 (legacy)**; collecting all stations at the end is free in the
  record style (the stations *are* the bindings).
- A `temp_ratio(gas, PR)` convenience **cannot exist** for a thermally
  perfect gas: at PR = 12 the temperature ratio varies from 2.024 to 1.834
  (**10.4%**) over T1 ∈ [250, 900] K — TR is a function of (T1, PR), not of
  PR alone, so the verb must take the inlet state.

## Decision

1. **Three-process taxonomy**, exported as verbs on the pure core
   (`FrozenGas` and `FastFrozenGas`, scalar kernels `(gas, T1, ...) -> T2`):
   - *Ratio-specified* (isentrope family): `compress(gas, T1, PR; ηp)` and
     `expand(gas, T1, PR; ηp)`. **Both take PR ≥ 1** (`ArgumentError`
     otherwise): the direction lives in the verb, never in the number.
     `expand` matches the legacy expansion convention exactly:
     `s0(T2) = s0(T1) + R·ηp·ln(1/PR)` (= legacy `expand(gas, 1/PR, ηp)`).
   - *Work-specified*: `add_work(st, w; ηp)` / `extract_work(st, w; ηp)`
     with `w ≥ 0`; enthalpy change ±w (T from the h-inversion), pressure
     along the polytrope `P2 = P1·exp(K/R·Δs0)` with K = ηp (adding) or
     1/ηp (extracting) — the two caller-side conventions of the legacy
     `set_Δh!(gas, ±Δh, ηp-or-1/ηp)`, now owned by the verbs.
   - *Heat at constant pressure*: `add_heat(st, q)`, signed `q`.
2. **Inversion is separate from process.** The enthalpy → temperature
   inversion is its own named function, `T_from_h(gas, hspec)`; the
   isentrope/polytrope is a *process*, expressed by `compress`/`expand`, never
   by an inversion verb carrying a `PR`/`ηp`. Process ≠ inversion.
3. **`GasState{G,F}`: an immutable (gas, T, P) value record** — ergonomics,
   not architecture. The substance stays a pure set of property curves
   (ADR-0001/0002); the record only lets the caller's (T, P) pair travel
   together through a process chain. It is `isbits` for
   `G = FrozenGas{Float64}`, never mutated — every verb returns a NEW
   state. State-layer verbs: `compress(st, PR; ηp)`, `expand(st, PR; ηp)`,
   `expand_to(st, P2; ηp)` (nozzle convenience, P2 ≤ st.P), `add_heat`,
   `add_work`, `extract_work`. Accessor *functions* (no getproperty
   magic): `cp/h/s0/gamma/R` at `st.T`, plus the two quantities only a
   (T, P) pair enables: `entropy(st) = s0 − R·ln(P/Pstd)` and
   `density(st) = P/(R·T)`. Full words are exported; `s`/`rho` remain
   unexported aliases (too collision-prone to export).
4. **Name coexistence**: the legacy unexported
   `compress/expand(gas::AbstractGas, PR, ηp)` methods on the mutable layer
   keep working under the now-exported names — dispatch on the first
   argument separates them; no behavior change for legacy callers.
5. ForwardDiff: the state verbs carry `Dual`s in T/P/PR/w/q through the
   parametric `F` and the existing extension IFT rules on
   `T_from_h`/`_T_polytropic`; **no new extension methods were needed**.

## Consequences

- Cycle code reads as the process chain it is
  (`compress → add_heat → extract_work → expand_to`), each station an
  immutable value; the wrong-P-rail class of bug is structurally gone for
  ~1.3% overhead and zero allocations.
- Both-ratios-≥-1 is a hard convention: translating an old `PR < 1` isentrope
  call means switching verb and inverting the ratio (a `PR = 0.25` expansion
  becomes `expand(g, T1, 4.0)`). At ηp = 1 the
  numbers are identical; with ηp ≠ 1 the expansion verb applies the
  expansion convention (R·ηp·ln(1/PR)), which is what the legacy `expand`
  always did.
- Do not re-propose `temp_ratio(gas, PR)` or constant-γ shortcuts on the
  pure core: measured 10.4% TR variation over inlet temperature at fixed
  PR (see Context).
- `GasState` deliberately stores no derived properties — caching "the
  current state's cp" is the pattern ADR-0001 forbids; accessors recompute
  from the curves (~ns scale).
