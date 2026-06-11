# Head-to-head: legacy mutable architecture (Gas{N}, Gas1D — as on public
# main) vs the new pure FrozenGas core, across realistic call patterns.
# Run:  julia --project=claude_sandbox/proto_deriv claude_sandbox/walkthrough/bench_compare.jl
# Emits claude_sandbox/walkthrough/benchdata.json

using IdealGasThermo, BenchmarkTools, ForwardDiff, Printf
const IGT = IdealGasThermo

# ---- fixtures ---------------------------------------------------------------
using Random
Random.seed!(42)
const TT = collect(range(250.0, 2200.0, length = 1000))[randperm(1000)] # unpredictable branch
const NCALL = length(TT)

frozen_air() = FrozenGas(DryAir)

results = Vector{Tuple{String,Float64,Float64,Float64}}() # key, min ns, median ns, allocs

function record!(key, trial; n = 1)
    t_min = minimum(trial).time / n
    t_med = median(trial).time / n
    al = minimum(trial).allocs / n
    push!(results, (key, t_min, t_med, al))
    @printf "%-46s %10.2f ns  %10.2f ns  %8.3f allocs\n" key t_min t_med al
end

println("=== scenario benchmarks (per-call, n = $NCALL sweep) ===")

# ---- 1. set T, read cp + h + ϕ (the cycle-deck inner pattern) ---------------
function sweep_all_gasN(g, TT)
    s = 0.0
    for T in TT
        g.T = T
        s += g.cp + g.h + g.ϕ
    end
    s
end
function sweep_all_frozen(air, TT)
    s = 0.0
    for T in TT
        p = props(air, T)
        s += p.cp + p.h + p.s0
    end
    s
end
record!("T-sweep, read cp+h+phi | legacy Gas{N}",
    @benchmark(sweep_all_gasN(g, $TT), setup = (g = Gas()), evals = 1); n = NCALL)
record!("T-sweep, read cp+h+phi | legacy Gas1D",
    @benchmark(sweep_all_gasN(g, $TT), setup = (g = Gas1D()), evals = 1); n = NCALL)
record!("T-sweep, read cp+h+phi | FrozenGas props",
    @benchmark(sweep_all_frozen(air, $TT), setup = (air = frozen_air()), evals = 1); n = NCALL)

# ---- 2. set T, read cp only -------------------------------------------------
function sweep_cp_gas(g, TT)
    s = 0.0
    for T in TT
        g.T = T
        s += g.cp
    end
    s
end
sweep_cp_frozen(air, TT) = sum(T -> IGT.cp(air, T), TT)
record!("T-sweep, read cp only  | legacy Gas{N}",
    @benchmark(sweep_cp_gas(g, $TT), setup = (g = Gas()), evals = 1); n = NCALL)
record!("T-sweep, read cp only  | legacy Gas1D",
    @benchmark(sweep_cp_gas(g, $TT), setup = (g = Gas1D()), evals = 1); n = NCALL)
record!("T-sweep, read cp only  | FrozenGas cp",
    @benchmark(sweep_cp_frozen(air, $TT), setup = (air = frozen_air()), evals = 1); n = NCALL)

# ---- 3. enthalpy inversion (burner / work-balance pattern) ------------------
# targets chosen so every solve is in-range; legacy Gas{N} works in J/mol,
# Gas1D and FrozenGas in J/kg
const T_targets = collect(range(300.0, 2100.0, length = 200))[randperm(200)]
function sweep_seth!(g, htargets)
    s = 0.0
    for ht in htargets
        set_h!(g, ht)
        s += g.T
    end
    s
end
sweep_Tofh(air, htargets) = sum(ht -> T_of_h(air, ht), htargets)

let gN = Gas(), g1 = Gas1D(), air = frozen_air()
    h_molar = map(T -> (gN.T = T; gN.h), T_targets)
    h_mass1 = map(T -> (g1.T = T; g1.h), T_targets)
    h_mass2 = map(T -> IGT.h(air, T), T_targets)
    record!("enthalpy inversion     | legacy Gas{N} set_h!",
        @benchmark(sweep_seth!(g, $h_molar), setup = (g = Gas()), evals = 1); n = 200)
    record!("enthalpy inversion     | legacy Gas1D set_h!",
        @benchmark(sweep_seth!(g, $h_mass1), setup = (g = Gas1D()), evals = 1); n = 200)
    record!("enthalpy inversion     | FrozenGas T_of_h",
        @benchmark(sweep_Tofh(air, $h_mass2), setup = (air = frozen_air()), evals = 1); n = 200)
end

# ---- 4. isentropic compression (PR = 12 from 288.15 K) ----------------------
function one_compress!(g, PR)
    set_TP!(g, 288.15, 101325.0)
    IGT.compress(g, PR)
    g.T
end
record!("isentropic compression | legacy Gas{N} compress",
    @benchmark(one_compress!(g, 12.0), setup = (g = Gas()), evals = 1))
record!("isentropic compression | legacy Gas1D compress",
    @benchmark(one_compress!(g, 12.0), setup = (g = Gas1D()), evals = 1))
record!("isentropic compression | FrozenGas T_isentropic",
    @benchmark(T_isentropic(air, 288.15, 12.0), setup = (air = frozen_air())))

# ---- 5. derivatives: dh/dT with N partials ----------------------------------
# legacy: structurally impossible (::Float64 pins + MVector{8,Float64} cache)
legacy_ad_supported = try
    g = Gas1D()
    ForwardDiff.derivative(t -> (g.T = t; g.h), 600.0)
    true
catch
    false
end
println("legacy path ForwardDiff support: $legacy_ad_supported")

# generic fallback (bypasses the extension rules by calling kernels directly)
h_generic(air, T) = IGT.Runiv * IGT.poly_h_R(IGT.coeffs(air, T), T, log(T))

mkdual(T, ::Val{N}) where {N} = ForwardDiff.Dual(T, ntuple(_ -> 1.0, Val(N)))
sweep_dual(f, air, TT, v) = sum(T -> ForwardDiff.value(f(air, mkdual(T, v))), TT)

for N in (1, 8, 12)
    record!("dh/dT, $(lpad(N,2)) partials     | FrozenGas generic Dual",
        @benchmark(sweep_dual(h_generic, air, $TT, $(Val(N))),
            setup = (air = frozen_air()), evals = 1); n = NCALL)
    record!("dh/dT, $(lpad(N,2)) partials     | FrozenGas analytic rule",
        @benchmark(sweep_dual(IGT.h, air, $TT, $(Val(N))),
            setup = (air = frozen_air()), evals = 1); n = NCALL)
end

# ---- 6. inversion under Dual (solver pattern, 8 partials) -------------------
function Tofh_generic(air, hd; Tguess = 500.0)
    T = one(hd) / oneunit(hd) * Tguess
    for _ = 1:30
        dT = (hd - h_generic(air, T)) / (IGT.Runiv * IGT.poly_cp_R(IGT.coeffs(air, T), T))
        T += dT
        abs(dT) <= 1e-12 * abs(T) && return T
    end
    error("no convergence")
end
let air = frozen_air()
    hvals = map(T -> IGT.h(air, T), T_targets)
    record!("T_of_h, 8 partials     | Dual through Newton loop",
        @benchmark(sum(hh -> ForwardDiff.value(Tofh_generic(air, mkdual(hh, Val(8)))), $hvals),
            setup = (air = frozen_air()), evals = 1); n = 200)
    record!("T_of_h, 8 partials     | IFT rule (extension)",
        @benchmark(sum(hh -> ForwardDiff.value(T_of_h(air, mkdual(hh, Val(8)))), $hvals),
            setup = (air = frozen_air()), evals = 1); n = 200)
end

# ---- emit JSON ---------------------------------------------------------------
out = joinpath(@__DIR__, "benchdata.json")
open(out, "w") do io
    println(io, "{")
    println(io, "  \"julia\": \"$(VERSION)\",")
    println(io, "  \"legacy_forwarddiff_supported\": $legacy_ad_supported,")
    println(io, "  \"scenarios\": [")
    for (i, (k, tmin, tmed, al)) in enumerate(results)
        comma = i == length(results) ? "" : ","
        @printf io "    {\"key\": \"%s\", \"min_ns\": %.3f, \"median_ns\": %.3f, \"allocs\": %.4f}%s\n" k tmin tmed al comma
    end
    println(io, "  ]")
    println(io, "}")
end
println("\nwrote $out")
