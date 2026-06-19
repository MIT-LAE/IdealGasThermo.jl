# PkgBenchmark regression suite for the pure, immutable core (v2 architecture):
# FrozenGas / FastFrozenGas / GasState property reads, inversions, process and
# flow verbs, combustion, mixing, humidity, and the analytic-derivative
# (ForwardDiff extension) paths.
#
# The legacy mutable layer (Gas{N}/Gas1D and the Dict-combustion) is deprecated
# and removed at v2.0 (ADR-0002/0007), so it is intentionally NOT tracked here;
# the historical old-vs-new comparison lives in benchmark/arch_comparison/.
#
# Run:  julia --project=benchmark -e 'using PkgBenchmark; benchmarkpkg("IdealGasThermo")'
# All fixtures are built once and `$`-interpolated into each @benchmarkable so the
# measured work is the verb itself, not global lookup (BenchmarkTools #269).

using BenchmarkTools
using IdealGasThermo
using ForwardDiff
using Random

const IGT = IdealGasThermo
const SUITE = BenchmarkGroup()

# ---- fixtures ---------------------------------------------------------------
Random.seed!(42)
# unpredictable branch order that straddles the 1000 K NASA-9 coefficient seam
const TT = collect(range(250.0, 2200.0, length = 256))[randperm(256)]
# in-range enthalpy targets for the inversion sweeps
const HT = collect(IGT.h(FrozenGas(DryAir), T) for T in range(300.0, 2100.0, length = 200))

const AIR       = FrozenGas(DryAir)
const FG_SEEDED = FastFrozenGas(AIR)                 # exact, table-seeded Newton
const FG_FAST   = FastFrozenGas(AIR, mode = :fast)   # pure Hermite lookup
const ST        = GasState(AIR, 288.15, 101325.0)
const COMB      = Vitiator("CH4", DryAir)
# `mix` blends two FrozenGas directly (each now carries its composition `gas.X`,
# so no precomputed mixing system is needed). Here: dry air + a CO2 stream
# (EGR-like).
const CO2_G     = FrozenGas(species_in_spdict("CO2"))

# ---- 1. property reads (the cycle-deck inner loop) --------------------------
sweep_props(gas, TT) = (s = 0.0; for T in TT; p = props(gas, T); s += p.cp + p.h + p.s0; end; s)
sweep_cp(gas, TT)    = sum(T -> IGT.cp(gas, T), TT)
SUITE["properties"]["props cp+h+s0 (T-sweep)"] = @benchmarkable sweep_props($AIR, $TT)
SUITE["properties"]["cp (T-sweep)"]            = @benchmarkable sweep_cp($AIR, $TT)
SUITE["properties"]["gamma"]                   = @benchmarkable IGT.gamma($AIR, 1200.0)

# ---- 2. inversions: enthalpy and isentrope ----------------------------------
sweep_Tofh(gas, HT) = sum(hh -> T_from_h(gas, hh), HT)
SUITE["inversions"]["T_from_h | FrozenGas (Newton)"]          = @benchmarkable sweep_Tofh($AIR, $HT)
SUITE["inversions"]["T_from_h | FastFrozenGas :seeded"]       = @benchmarkable sweep_Tofh($FG_SEEDED, $HT)
SUITE["inversions"]["T_from_h | FastFrozenGas :fast"]         = @benchmarkable sweep_Tofh($FG_FAST, $HT)
SUITE["inversions"]["T_polytropic | FrozenGas"]             = @benchmarkable IGT._T_polytropic($AIR, 288.15, 12.0)
SUITE["inversions"]["T_polytropic | FastFrozenGas :seeded"] = @benchmarkable IGT._T_polytropic($FG_SEEDED, 288.15, 12.0)

# ---- 3. process verbs (GasState) --------------------------------------------
SUITE["verbs"]["compress (ηp)"]   = @benchmarkable compress($ST, 12.0; ηp = 0.9)
SUITE["verbs"]["compress (ηs)"]   = @benchmarkable compress($ST, 12.0; ηs = 0.85)
SUITE["verbs"]["expand (ηp)"]     = @benchmarkable expand($ST, 4.0; ηp = 0.92)
SUITE["verbs"]["add_heat"]        = @benchmarkable add_heat($ST, 3.0e5)
SUITE["verbs"]["add_work (ηp)"]   = @benchmarkable add_work($ST, 2.0e5; ηp = 0.9)

# ---- 4. gas dynamics (flow verbs) -------------------------------------------
SUITE["flow"]["speed_of_sound"]   = @benchmarkable speed_of_sound($AIR, 288.15)
SUITE["flow"]["stagnation_state"] = @benchmarkable stagnation_state($ST, 0.8)
SUITE["flow"]["static_state"]     = @benchmarkable static_state($ST, 0.8)

# ---- 5. combustion, mixing, humidity ----------------------------------------
SUITE["combustion"]["products (CH4/air, FAR=0.03)"] = @benchmarkable products($COMB, 0.03)
SUITE["mixing"]["mix (mratio=0.25)"]                = @benchmarkable mix($AIR, $CO2_G, 0.25)
SUITE["humidity"]["humid_air (RH=0.6)"]             = @benchmarkable humid_air(RH = 0.6, T = 288.15, P = 101325.0)

# ---- 6. automatic differentiation (analytic IFT extension) ------------------
# Forward properties are linear in the coefficients (cheap, flat in #partials);
# inversions use the implicit-function-theorem rules, never loop differentiation.
mkdual(x, ::Val{N}) where {N} = ForwardDiff.Dual(x, ntuple(_ -> 1.0, Val(N)))
sweep_dh(gas, TT, v) = sum(T -> ForwardDiff.value(IGT.h(gas, mkdual(T, v))), TT)
SUITE["autodiff"]["dh/dT, 1 partial (T-sweep)"]   = @benchmarkable sweep_dh($AIR, $TT, $(Val(1)))
SUITE["autodiff"]["dh/dT, 8 partials (T-sweep)"]  = @benchmarkable sweep_dh($AIR, $TT, $(Val(8)))
SUITE["autodiff"]["d(T_from_h), 8 partials (IFT)"]  =
    @benchmarkable sum(hh -> ForwardDiff.value(T_from_h($AIR, mkdual(hh, Val(8)))), $HT)
# Dual-carrying substance: derivative of an outlet property w.r.t. FAR, where the
# product composition (hence the gas's own coefficients) carry the FAR-tangent.
SUITE["autodiff"]["d(products)/dFAR"] =
    @benchmarkable ForwardDiff.derivative(far -> IGT.h(products($COMB, far), 1600.0), 0.03)
