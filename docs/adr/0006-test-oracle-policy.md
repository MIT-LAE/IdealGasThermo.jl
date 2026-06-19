# ADR-0006: The pure core is validated by property tests + external anchors, not by agreement with the legacy mutable types

Date: 2026-06-17
Status: accepted

## Context

While the pure core (`FrozenGas`, `GasState`, `Vitiator`/`products`,
`mix`, `humid_air`, `static_state`/`stagnation_state`) was being
built, its tests were written as **agreement with the legacy mutable layer**:
`FrozenGas.cp ≈ Gas1D.cp`, `entropy(st) ≈ Gas1D.s`, `compress(st) ≈
PressureRatio(Gas1D)`, `products ≈ vitiated_species`, and so on. That was the
right scaffolding *during migration* — it proved the new core reproduced the
established behaviour — but it has two problems as a permanent test surface:

1. **It is circular.** `Gas1D`/`Gas{N}` and `FrozenGas` evaluate the *same*
   NASA-9 polynomials from the *same* `data/thermo.inp` coefficients through the
   *same* entropy-of-mixing fold. "FrozenGas agrees with Gas1D" is two arrangements
   of one implementation agreeing with itself — it catches a transcription slip,
   not a physics error, and it cannot survive the deletion of `Gas1D`/`Gas{N}`
   (candidate 3 / v2.0, ADR-0002).
2. **It hid quantity-over-quality.** An adversarial audit (2026-06-17) found the
   suite's green count overstated its distinct-claim count ~6× — a CEA loop
   sampled ~67 temperatures where ~7 carry the same power, `GasState` accessor
   `===` forwards and "wrapper equals the body it calls" assertions, allocation
   checks padding a *correctness* count, and duplicated agreement testsets.
   Worse, it surfaced two **blind spots** that no amount of agreement-testing
   closes: a global enthalpy-datum shift, and a global entropy-of-mixing sign
   flip (invisible to every pure-species check, where `X·ln X = 0`).

## Decision

1. **The permanent oracle is reference-implementation-independent.**
   `test/unit_test_properties.jl` holds the load-bearing checks: metamorphic
   (`∫cp` vs `h`/`s0` by quadrature), conservation (combustion mass balance,
   heat-release ≈ `f·LHV`), invariants (second-law entropy generation, process
   round-trips, `s(T,kP) = s(T,P) − R·ln k`, mixing symmetry), and AD-vs-FD.
   Each is a distinct physical claim, not a coverage statistic.
2. **Absolute values are anchored to an external authority, not to a twin.**
   `test/unit_test_cea_reference.jl` checks `cp`/`s0`/`h` of CO2/N2/O2 and the
   dry-air pseudo-species against **CEA** (`test/CEA_output.txt`), which
   evaluates the *same* NASA-9 coefficients (verified identical to
   `data/thermo.inp`). The comparison is in CEA's printed molar units to its
   printed precision: the error model is **absolute** (CEA prints 3 decimals),
   so `atol = 1 unit in the last place (1e-3)` — the 5e-4 half-ULP is CEA's
   rounding, the small remainder is the `Runiv = 8.3145` vs CEA's `8.31451`
   difference (heaviest on `s0`). Curated temperatures straddle the 1000 K
   coefficient seam; dense sampling of a degree-≤4 polynomial is padding.
3. **The pure core is NOT validated by agreement with `Gas`/`Gas1D`.** Those
   assertions are deleted from the pure-core test files
   (`unit_test_frozengas`, `unit_test_gasstate`, `unit_test_products`,
   `unit_test_mixing`, `unit_test_humidair`, `unit_test_flow`) — dropped when a
   property/CEA check already covers the claim, or converted to a self-contained
   contract (e.g. `entropy` at `Pstd` equals `s0`; `add_work` energy balance and
   polytropic-pressure convention). Wrapper-equals-body and accessor-`===`-forward
   tautologies are removed; allocation checks are kept but understood as
   *performance*, not correctness, contracts.
4. **Blind spots get a dedicated test, not more agreement tests.** The
   entropy-of-mixing sign is pinned in `unit_test_properties.jl` against a
   hand-written `Σ Xᵢ·s0ᵢ − R·Σ Xᵢ·ln Xᵢ` (not via `generate_composite_species`).
   The enthalpy datum is anchored by the CEA datum check plus the independent
   `release ≈ f·LHV` cross-oracle.
5. **Legacy types keep a minimal smoke test until v2.0** (ADR-0002): `Gas`/`Gas1D`
   remain exercised in the legacy-only files (`unit_test_turbo`,
   `unit_test_vitiated`, `unit_test_composite`, `unit_test_mixthermo`) so the
   deprecated path stays correct while it ships. These do **not** validate the
   pure core; they test the legacy types against pinned/legacy values and are
   deleted with the mutable layer at v2.0.

## Consequences

- The pure core's correctness no longer depends on the types it replaces; the
  v2.0 deletion of `Gas`/`Gas1D`/legacy combustion/turbo removes the legacy
  smoke files without touching the pure-core oracle.
- Do not re-introduce "≈ Gas1D"/"≈ vitiated_species"/"≈ set_Δh!" agreement
  assertions in the pure-core files. If a new pure feature needs an absolute
  check, anchor it to CEA (same-polynomial, tight) or to a documented external
  reference — not to another in-package implementation.
- Tolerances state their basis: CEA's printed precision (`atol`), the NASA-9
  fit residual, or machine epsilon — never a loose round number chosen to make
  a test pass. The TASOPT pins in `unit_test_turbo` are a *loose golden-master
  of `Gas{N}`* (legacy), not an authority, and are not re-pointed onto the pure
  core.
- Green count is not a goal. A new test must catch a bug a named source change
  would introduce; if an existing test already catches it, do not add another.
