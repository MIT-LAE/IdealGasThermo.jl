# External absolute anchor: per-species cp, h, s0 of CO2 / N2 / O2 against CEA.
#
# test/CEA_output.txt is CEA's own evaluation of the SAME NASA-9 coefficients
# (from https://cearun.grc.nasa.gov/ThermoBuild/) that
# this package loads: the fitted coefficients CEA prints in that file
# are identical to data/thermo.inp (e.g. CO2 alow = 49436.5054, -626.411601,
# 5.30172524, ...). So this is a TIGHT check that FrozenGas evaluates those
# polynomials correctly — the mass-scaling (×1000/MW), the formation-inclusive
# enthalpy datum, and the entropy integration constant — against an independent,
# authoritative evaluator.
#
# CEA columns (OPTIONS: joules), molar basis:
#   T[K]  Cp[J/mol/K]  H-H298[kJ/mol]  S[J/mol/K]  -(G-H298)/T  H[kJ/mol]  ...
# Mass-basis conversions (MW in g/mol):
#   cp = Cp·1000/MW [J/kg/K];  s0 = S·1000/MW [J/kg/K];
#   sensible h = (H−H298)·1e6/MW [J/kg]  (the formation datum is checked separately)

@testset "CEA per-species reference (CO2/N2/O2, same NASA-9 coeffs)" begin

    # parse CEA_output.txt -> species => Vector of (T, Cp, Hsens, S), finite rows only
    function parse_cea(path)
        data = Dict{String,Vector{NTuple{4,Float64}}}()
        cur = ""
        for line in eachline(path)
            m = match(r"CALCULATED FROM COEFFICIENTS FOR\s+(\S+)", line)
            if m !== nothing
                cur = String(m.captures[1])
                data[cur] = NTuple{4,Float64}[]
                continue
            end
            isempty(cur) && continue
            toks = split(strip(line))
            length(toks) >= 4 || continue
            T, Cp, Hs, S = tryparse.(Float64, toks[1:4])
            any(isnothing, (T, Cp, Hs, S)) && continue   # header/blank/INFINITE rows
            push!(data[cur], (T, Cp, Hs, S))
        end
        data
    end

    cea = parse_cea(joinpath(@__DIR__, "CEA_output.txt"))
    @test Set(keys(cea)) == Set(["CO2", "N2", "O2"])

    # Compare in CEA's own printed (molar) units to its printed precision. CEA
    # prints 3 decimals (Cp, S in J/mol/K; H in kJ/mol), so the error model is
    # ABSOLUTE: atol = one unit in the last printed place (1e-3). The half-ULP
    # (5e-4) is CEA's rounding; the small remainder is the universal-gas-constant
    # difference (Runiv = 8.3145 here vs CEA's 8.31451), heaviest on s0. Measured
    # max abs dev: cp 5.6e-4, s0 8.7e-4, sensible-h 6.1e-4 — all inside 1e-3.
    # (rtol would hide this structure; atol names it.)
    #
    # Curated temperatures, NOT the whole table: NASA-9 cp is a degree-≤4
    # polynomial per interval, so dense sampling adds count, not bug-catching
    # power (test audit, 2026-06-17). These bracket the 1000 K coefficient seam
    # (950 uses alow, 1000/1050 use ahigh) and hit both interval interiors and
    # the endpoints of the tested range — the points where a coefficient or
    # interval-selection bug actually shows.
    targets = (250.0, 500.0, 950.0, 1000.0, 1050.0, 1500.0, 3000.0)
    for spp in ("CO2", "N2", "O2")
        sp = species_in_spdict(spp)
        gas = FrozenGas(sp)
        h298 = IdealGasThermo.h(gas, 298.15)
        rows = Dict(T => (Cp, Hs, S) for (T, Cp, Hs, S) in cea[spp])

        # formation-inclusive datum: h(298.15) reproduces the stored formation
        # enthalpy, which equals CEA's H(298.15) column (ΔHf): −393.510 kJ/mol
        # for CO2, 0 for the elements N2/O2. The CO2 residual (4.7e-4 kJ/mol) is
        # the NASA-9 header-vs-coefficient consistency, not a CEA difference.
        @test h298 * sp.MW / 1e6 ≈ sp.Hf / 1000 atol = 1e-3   # kJ/mol

        for T in targets
            haskey(rows, T) ||
                error("CEA table for $spp is missing the curated T = $T K")
            Cp, Hs, S = rows[T]
            @test IdealGasThermo.cp(gas, T) * sp.MW / 1000 ≈ Cp atol = 1e-3         # J/mol/K
            @test IdealGasThermo.s0(gas, T) * sp.MW / 1000 ≈ S atol = 1e-3          # J/mol/K
            @test (IdealGasThermo.h(gas, T) - h298) * sp.MW / 1e6 ≈ Hs atol = 1e-3  # kJ/mol, sensible
        end
    end
end
