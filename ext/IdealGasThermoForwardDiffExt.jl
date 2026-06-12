"""
Analytic-derivative fast paths for `FrozenGas` under ForwardDiff.

The base package is generic over `Real`, so ForwardDiff works without this
extension — these rules replace Dual arithmetic through the polynomials with
closed forms (dh/dT = cp, ds0/dT = cp/T, dcp/dT analytic) and replace
differentiation of the Newton inversions with implicit-function-theorem
rules. Benchmarked ~4x faster at 12 partials for properties and ~2x for
inversions (see docs/adr/0001-pure-immutable-gas-core.md).
"""
module IdealGasThermoForwardDiffExt

using IdealGasThermo
using IdealGasThermo: FrozenGas, TabulatedGas, coeffs, Runiv, cp, h, s0
import IdealGasThermo: cp, h, s0, gamma, props, T_of_h, T_isentropic
using ForwardDiff: Dual, value, partials

# dcp/dT [J/kg/K²], closed form from the NASA-9 cp polynomial
@inline function cp_dT(gas::FrozenGas, T)
    a = coeffs(gas, T)
    Runiv *
    ((-2 * a[1] / T - a[2]) / T / T + a[4] + T * (2 * a[5] + T * (3 * a[6] + T * (4 * a[7]))))
end

cp(gas::FrozenGas, d::Dual{Tag}) where {Tag} =
    Dual{Tag}(cp(gas, value(d)), cp_dT(gas, value(d)) * partials(d))

h(gas::FrozenGas, d::Dual{Tag}) where {Tag} =
    Dual{Tag}(h(gas, value(d)), cp(gas, value(d)) * partials(d))

s0(gas::FrozenGas, d::Dual{Tag}) where {Tag} =
    Dual{Tag}(s0(gas, value(d)), cp(gas, value(d)) / value(d) * partials(d))

function props(gas::FrozenGas, d::Dual{Tag}) where {Tag}
    T = value(d)
    p = props(gas, T)
    ∂T = partials(d)
    (
        cp = Dual{Tag}(p.cp, cp_dT(gas, T) * ∂T),
        h = Dual{Tag}(p.h, p.cp * ∂T),
        s0 = Dual{Tag}(p.s0, p.cp / T * ∂T),
    )
end

# Inversions by the implicit function theorem: solve on the primal in plain
# Float64, attach exact partials at the solution — never differentiate the
# Newton loop.

# h(T) = hspec  ⟹  dT = dh / cp(T)
function T_of_h(gas::FrozenGas, hd::Dual{Tag}; Tguess = 500.0) where {Tag}
    T = T_of_h(gas, value(hd); Tguess = Tguess)
    Dual{Tag}(T, partials(hd) / cp(gas, T))
end

# s0(T2) = s0(T1) + R·ln(PR)/ηp  ⟹
#   cp(T2)/T2 · dT2 = cp(T1)/T1 · dT1 + R/(ηp·PR) · dPR
function T_isentropic(gas::FrozenGas, T1::Dual{Tag}, PR::Real; ηp = 1.0) where {Tag}
    T2 = T_isentropic(gas, value(T1), PR; ηp = ηp)
    ∂ = cp(gas, value(T1)) / value(T1) * partials(T1)
    Dual{Tag}(T2, T2 / cp(gas, T2) * ∂)
end

function T_isentropic(gas::FrozenGas, T1::Real, PR::Dual{Tag}; ηp = 1.0) where {Tag}
    T2 = T_isentropic(gas, T1, value(PR); ηp = ηp)
    ∂ = gas.R / (ηp * value(PR)) * partials(PR)
    Dual{Tag}(T2, T2 / cp(gas, T2) * ∂)
end

function T_isentropic(gas::FrozenGas, T1::Dual{Tag}, PR::Dual{Tag}; ηp = 1.0) where {Tag}
    T2 = T_isentropic(gas, value(T1), value(PR); ηp = ηp)
    ∂ = cp(gas, value(T1)) / value(T1) * partials(T1) +
        gas.R / (ηp * value(PR)) * partials(PR)
    Dual{Tag}(T2, T2 / cp(gas, T2) * ∂)
end

# TabulatedGas: identical implicit-function-theorem rules — the primal solve
# goes through the table-seeded (or out-of-range fallback) Newton path on
# plain Float64; the partials are the same closed forms evaluated with the
# wrapped gas's exact polynomials.

function T_of_h(tg::TabulatedGas, hd::Dual{Tag}) where {Tag}
    T = T_of_h(tg, value(hd))
    Dual{Tag}(T, partials(hd) / cp(tg.gas, T))
end

function T_isentropic(tg::TabulatedGas, T1::Dual{Tag}, PR::Real; ηp = 1.0) where {Tag}
    gas = tg.gas
    T2 = T_isentropic(tg, value(T1), PR; ηp = ηp)
    ∂ = cp(gas, value(T1)) / value(T1) * partials(T1)
    Dual{Tag}(T2, T2 / cp(gas, T2) * ∂)
end

function T_isentropic(tg::TabulatedGas, T1::Real, PR::Dual{Tag}; ηp = 1.0) where {Tag}
    gas = tg.gas
    T2 = T_isentropic(tg, T1, value(PR); ηp = ηp)
    ∂ = gas.R / (ηp * value(PR)) * partials(PR)
    Dual{Tag}(T2, T2 / cp(gas, T2) * ∂)
end

function T_isentropic(tg::TabulatedGas, T1::Dual{Tag}, PR::Dual{Tag}; ηp = 1.0) where {Tag}
    gas = tg.gas
    T2 = T_isentropic(tg, value(T1), value(PR); ηp = ηp)
    ∂ = cp(gas, value(T1)) / value(T1) * partials(T1) +
        gas.R / (ηp * value(PR)) * partials(PR)
    Dual{Tag}(T2, T2 / cp(gas, T2) * ∂)
end

end # module
