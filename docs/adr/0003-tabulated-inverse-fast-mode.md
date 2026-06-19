# ADR-0003: `FastFrozenGas` — tabulated inverse "fast mode", two tiers, mandatory bounds guard

Date: 2026-06-12
Status: accepted (naming settled 2026-06-12; refreshed 2026-06-19 to the final
`FastFrozenGas` / `T_from_h` names — see the note at the end of Context)

## Context

The pure Newton inversions (`T_from_h`, and the internal isentropic/polytropic
engine `_T_polytropic`) cost ~120–150 ns from a cold start (~5 iterations at ~25 ns
each: one h + one cp evaluation per step). For a fixed composition, the inverse
functions T(h) and T(s0) can be precomputed. Prototyped and measured 2026-06-11/12
(dry air, T ∈ [200, 2400] K, Apple Silicon; code in `claude_sandbox/proto_invfit/`):

- **Cubic Hermite tables** (uniform grid in h or s0, node slopes from the
  exact derivatives 1/cp and T/cp) dominate alternatives: N=256 gives
  max |ΔT/T| = 5.8e-10 (h) / 1.4e-9 (s0) at ~3 ns/call (50×), 4 KB, ~37 µs
  build. Linear tables need 8192 nodes for what Hermite does at 128 and
  have a discontinuous derivative (poison for outer Newton solvers);
  Chebyshev (deg 16) is compact but slower and stuck at ~5e-6.
- **Table-seeded Newton** (Hermite seed + polish with the exact equation)
  keeps the package's full inversion contract: 1 polish step → 3.3e-10 at
  ~13 ns (11×); convergence-checked polish → machine precision (~1e-15) at
  ~30 ns (5×).
- Interpolation error saturates at ~2e-10 regardless of N (float rounding
  floor), so a pure table can honestly advertise ≤1e-9 but never the
  package's 1e-12 Newton tolerance.
- **Consistency**: only the inverses are approximated — the forward
  functions are untouched, so both tables sample the same NASA-9
  polynomials and there is no second gas model. Measured cross-consistency
  of the h-table and s0-table inverses at the same physical state: 1.7e-9
  (≈ 2 µK at 1000 K); isentropic round trip closes to 5.7e-9.
- **Trap (measured)**: silent extrapolation outside the table range
  (T1 = 1500 K, PR = 40 → T2 ≈ 4300 K beyond Tmax = 2400) degrades the
  round trip to 1.9e-3 with no warning.

> **Naming history (folded in 2026-06-19).** This ADR originally named the wrapper
> `TabulatedGas`/`tabulate` and exposed a separate `*_interp` function pair and a
> `temperature(gas; …)` inversion facade. None of those names survived: the wrapper
> is `FastFrozenGas` ("tabulated" described the implementation, not the contract, and
> the default tier is exact), the interp tier is selected by the `:fast` *mode*
> rather than by a `*_interp` name, and the inversion the wrapper accelerates is the
> exported `T_from_h` (plus the internal `_T_polytropic`), the `temperature` facade
> having been dropped entirely (ADR-0004). The Decision below is stated in the final
> names.

## Decision

1. **The tier is selected by the type, not by the function name**, so an accelerated
   gas drops into existing call sites unchanged. `FastFrozenGas(gas; mode = :seeded |
   :fast, N = 256, Tmin = 200, Tmax = 2400)` wraps a `FrozenGas` plus the Hermite
   inverse tables (h → T and s0 → T). The forward property functions delegate to the
   wrapped gas unchanged; the inversions `T_from_h` and `_T_polytropic` dispatch on
   the `FastFrozenGas{mode}` type parameter.
2. **`:seeded` tier — exact (default).** The table provides only the Newton *seed*;
   the solve then polishes with the same convergence criterion as `FrozenGas` — the
   same documented contract, ~5–11× faster. ForwardDiff support via the same IFT
   extension rules on `T_from_h`/`_T_polytropic`.
3. **`:fast` tier — opt-in and loud.** Pure Hermite lookup (no polish), documented
   ≤1e-9 (N = 256), throwing `DomainError` outside the table range. This *is* the
   interpolation tier, selected by the `:fast` mode — there is no separately named
   `*_interp` function.
4. **Bounds guard is mandatory.** The `:seeded` tier falls back to the plain
   cold-start Newton solve when the target lies outside the table range; no code path
   may silently extrapolate the tables (the `:fast` tier throws rather than
   extrapolate).

## Consequences

- Hot loops doing repeated inversions on a fixed composition get 5–50×
  depending on tier, with the `:seeded` tier changing no accuracy contract.
- Round-trip identities through the `:fast` tier hold to ~1e-9, not 1e-12;
  code asserting machine-precision round trips must use the `:seeded` tier.
- Do not re-propose linear interpolation tables or global Chebyshev fits
  for these inverses (measured inferior, see Context).
