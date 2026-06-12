# Provenance: `data/thermo.inp`

## Source

NASA Glenn thermodynamic database in **NASA-9 polynomial format**, obtained
via the NASA ThermoBuild web interface
(<https://cearun.grc.nasa.gov/ThermoBuild/>). The file header carries the
database release date `9/09/04` (September 9, 2004), i.e. the NASA Glenn
coefficients described in McBride, Zehe & Gordon, *NASA Glenn Coefficients
for Calculating Thermodynamic Properties of Individual Species*,
NASA/TP-2002-211556 (2002), as distributed with CEA.

Retrieval history (this repository):

- 2022-03-31 — initial species set added (commit `061535b`).
- 2023-10-01 — additional species appended from the same ThermoBuild
  database (commit `0869a44`).

## Format and contents

- Two temperature intervals per species, split at **Tmid = 1000 K**
  (enforced by `readThermo`; it errors on any other Tmid). Typical ranges
  200–1000 K and 1000–6000 K.
- Nine coefficients per interval (`a1…a7` plus integration constants
  `b1`, `b2`) for the dimensionless forms cp/R, H/(RT), S/R.
- Per the format spec, the species header line carries the molecular
  weight [g/mol] and the heat of formation at 298.15 K [J/mol].
  Format reference: <https://shepherd.caltech.edu/EDL/PublicResources/sdt/formats/nasa.html>

## Units and conversions applied by this package

- `readThermo` stores coefficients and header values as given in the file
  (no rescaling at read time): MW [g/mol], Hf [J/mol].
- Molar → mass-specific conversion happens at `FrozenGas`/composite-species
  construction: properties are scaled by `Runiv / MW` to J/kg-based units.
- `Runiv = 8.3145 J/K/mol` (`src/constants.jl`), the value also used by the
  NASA-9 fits' era; not the full-precision CODATA 2018 value
  (8.31446261815324).
- Enthalpy datum is **formation-inclusive** (CEA-style): `h(gas, 298.15 K)`
  equals the mixture's mass-specific formation enthalpy, not zero.

## Standard-state note

The NASA-9 entropy integration constants define s° at a standard pressure
of **1 bar** (100 000 Pa). This package evaluates pressure-dependent
entropy as `s(T, P) = s0(T) − R·ln(P/Pstd)` with `Pstd = 101 325 Pa`
(1 atm, `src/constants.jl`). The two references differ by
`R·ln(101325/100000) ≈ 0.0132·R` per unit mass; this is inherited from the
package's original convention and is documented here rather than silently
changed.
