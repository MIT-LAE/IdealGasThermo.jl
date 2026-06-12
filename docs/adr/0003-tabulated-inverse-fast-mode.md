# ADR-0003: Tabulated inverse "fast mode" — two tiers, mandatory bounds guard

Date: 2026-06-12
Status: accepted

## Context

The pure Newton inversions (`T_of_h`, `T_isentropic`) cost ~120–150 ns from a
cold start (~5 iterations at ~25 ns each: one h + one cp evaluation per
step). For a fixed composition, the inverse functions T(h) and T(s0) can be
precomputed. Prototyped and measured 2026-06-11/12 (dry air, T ∈ [200,
2400] K, Apple Silicon; code in `claude_sandbox/proto_invfit/`):

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

## Decision

1. `TabulatedGas` wraps a `FrozenGas` plus Hermite inverse tables
   (`tabulate(gas; N=256, Tmin=200, Tmax=2400)`). Forward property
   functions delegate to the wrapped gas unchanged.
2. **Default tier is exact**: `T_of_h`/`T_isentropic` on a `TabulatedGas`
   use the table only as the Newton seed, then polish with the same
   convergence criterion as `FrozenGas` — same documented contract,
   ~5–11× faster. ForwardDiff support via the same IFT extension rules.
3. **Interp tier is opt-in and loud**: separately named `*_interp`
   functions, documented ≤1e-9 (N=256), throwing `DomainError` outside the
   table range.
4. **Bounds guard is mandatory**: the default tier falls back to the plain
   cold-start Newton solve when the target lies outside the table range —
   no code path may silently extrapolate the tables.

## Consequences

- Hot loops doing repeated inversions on a fixed composition get 5–50×
  depending on tier, with the default tier changing no accuracy contract.
- Round-trip identities through the interp tier hold to ~1e-9, not 1e-12;
  code asserting machine-precision round trips must use the default tier.
- Do not re-propose linear interpolation tables or global Chebyshev fits
  for these inverses (measured inferior, see Context).
