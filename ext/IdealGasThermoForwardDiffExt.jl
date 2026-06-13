"""
Analytic-derivative fast paths for `FrozenGas` under ForwardDiff.

The base package is generic over `Real`, so ForwardDiff works without this
extension — these rules replace Dual arithmetic through the polynomials with
closed forms (dh/dT = cp, ds0/dT = cp/T, dcp/dT analytic) and replace
differentiation of the Newton inversions with implicit-function-theorem
rules. Benchmarked ~4x faster at 12 partials for properties and ~2x for
inversions (see docs/adr/0001-pure-immutable-gas-core.md).

Inversions also cover a *Dual-carrying gas* (`FrozenGas{<:Dual}`, produced by
`products(sys, FAR::Dual)`): the full IFT rule adds the "composition moves"
term, while keeping the Newton loop on the value rail (the composition
tangent is one forward evaluation, since properties are linear in the
coefficients). See CONTEXT.md ("Dual-carrying gas").
"""
module IdealGasThermoForwardDiffExt

using IdealGasThermo
using IdealGasThermo: FrozenGas, FastFrozenGas, coeffs, Runiv, cp, h, s0
import IdealGasThermo: cp, h, s0, gamma, props, T_of_h, T_isentropic
using ForwardDiff: Dual, value, partials
using StaticArrays: SVector

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

# Dual-carrying substance (FrozenGas{<:Dual}) — the gas's own coefficients
# carry the tangent, as produced by `products(sys, FAR::Dual)` (composition
# depends on FAR). The rules above dispatch on the TARGET being a Dual and
# assume a constant substance; these dispatch on the GAS being Dual-typed and
# add the "composition moves" IFT term the constant-substance rules drop.
# Same split-rule speed: the Newton loop runs entirely on the value rail, and
# the composition tangent is one forward property evaluation (properties are
# linear in the coefficients) — never differentiation of the loop. See ADR-0004
# and CONTEXT.md ("Dual-carrying gas").

# value-strip a FrozenGas{<:Dual} to its plain Float64 value-part gas
@inline function _value_gas(gas::FrozenGas)
    FrozenGas(
        SVector{9,Float64}(value.(gas.alow)),
        SVector{9,Float64}(value.(gas.ahigh)),
        value(gas.MW), value(gas.R), value(gas.Hf),
    )
end

# Tangent extraction that also handles a plain-Float argument as identically
# zero. `false` is the additive identity: `false - p::Partials == -p`, with no
# allocation and no need to probe the partials length.
@inline _tangent(x::Dual) = partials(x)
@inline _tangent(::Real) = false
@inline _minus(a, b) = a - b
@inline _minus(::Bool, b) = -b   # zero (plain-Float) tangent minus a Partials
@inline _gas_tag(::FrozenGas{<:Dual{Tag}}, args...) where {Tag} = Tag

# T_of_h: F(T, p) = h(gas(p), T) - hspec(p) = 0  ⟹
#   ∂T = ( partials(hspec) - partials(h(gas, T*)) ) / cp(gas₀, T*)
function _T_of_h_dualgas(gas::FrozenGas, hspec; Tguess)
    Tag = _gas_tag(gas, hspec)
    gas0 = _value_gas(gas)
    Tstar = T_of_h(gas0, value(hspec); Tguess = Tguess) # expensive solve, pure Float64
    ∂comp = _tangent(h(gas, Tstar))                     # one Dual eval at float T*
    ∂targ = _tangent(hspec)
    Dual{Tag}(Tstar, _minus(∂targ, ∂comp) / cp(gas0, Tstar))
end
# Dual gas + plain-Float target
T_of_h(gas::FrozenGas{<:Dual{Tag}}, hspec::Real; Tguess = 500.0) where {Tag} =
    _T_of_h_dualgas(gas, hspec; Tguess = Tguess)
# Dual gas + Dual target (same tag) — more specific in arg 1 than the
# constant-substance rule above, so this resolves the would-be ambiguity.
T_of_h(gas::FrozenGas{<:Dual{Tag}}, hspec::Dual{Tag}; Tguess = 500.0) where {Tag} =
    _T_of_h_dualgas(gas, hspec; Tguess = Tguess)

# T_isentropic: G(T2, p) = s0(gas, T2) - [s0(gas, T1) + gas.R·ln(PR)/ηp] = 0  ⟹
#   ∂T2 = T2/cp(gas₀, T2) · ( partials(target) - partials(s0(gas, T2)) )
function _T_isentropic_dualgas(gas::FrozenGas, T1, PR; ηp)
    Tag = _gas_tag(gas, T1, PR)
    gas0 = _value_gas(gas)
    T2 = T_isentropic(gas0, value(T1), value(PR); ηp = value(ηp)) # pure Float64 solve
    target = s0(gas, T1) + gas.R * log(PR) / ηp                   # Dual eval, no Newton
    ∂ = _minus(_tangent(target), _tangent(s0(gas, T2)))
    Dual{Tag}(T2, T2 / cp(gas0, T2) * ∂)
end
T_isentropic(gas::FrozenGas{<:Dual{Tag}}, T1::Real, PR::Real; ηp = 1.0) where {Tag} =
    _T_isentropic_dualgas(gas, T1, PR; ηp = ηp)
T_isentropic(gas::FrozenGas{<:Dual{Tag}}, T1::Dual{Tag}, PR::Real; ηp = 1.0) where {Tag} =
    _T_isentropic_dualgas(gas, T1, PR; ηp = ηp)
T_isentropic(gas::FrozenGas{<:Dual{Tag}}, T1::Real, PR::Dual{Tag}; ηp = 1.0) where {Tag} =
    _T_isentropic_dualgas(gas, T1, PR; ηp = ηp)
T_isentropic(gas::FrozenGas{<:Dual{Tag}}, T1::Dual{Tag}, PR::Dual{Tag}; ηp = 1.0) where {Tag} =
    _T_isentropic_dualgas(gas, T1, PR; ηp = ηp)

# FastFrozenGas: identical implicit-function-theorem rules — the primal
# solve goes through the mode-appropriate tier on plain Float64 (:seeded =
# table-seeded Newton with out-of-range fallback; :fast = pure Hermite
# lookup, which approximates the same inverse to ~1e-9 so the exact IFT
# partials remain the principled tangent); the partials are the same closed
# forms evaluated with the wrapped gas's exact polynomials.

function T_of_h(fg::FastFrozenGas, hd::Dual{Tag}) where {Tag}
    T = T_of_h(fg, value(hd))
    Dual{Tag}(T, partials(hd) / cp(fg.gas, T))
end

function T_isentropic(fg::FastFrozenGas, T1::Dual{Tag}, PR::Real; ηp = 1.0) where {Tag}
    gas = fg.gas
    T2 = T_isentropic(fg, value(T1), PR; ηp = ηp)
    ∂ = cp(gas, value(T1)) / value(T1) * partials(T1)
    Dual{Tag}(T2, T2 / cp(gas, T2) * ∂)
end

function T_isentropic(fg::FastFrozenGas, T1::Real, PR::Dual{Tag}; ηp = 1.0) where {Tag}
    gas = fg.gas
    T2 = T_isentropic(fg, T1, value(PR); ηp = ηp)
    ∂ = gas.R / (ηp * value(PR)) * partials(PR)
    Dual{Tag}(T2, T2 / cp(gas, T2) * ∂)
end

function T_isentropic(fg::FastFrozenGas, T1::Dual{Tag}, PR::Dual{Tag}; ηp = 1.0) where {Tag}
    gas = fg.gas
    T2 = T_isentropic(fg, value(T1), value(PR); ηp = ηp)
    ∂ = cp(gas, value(T1)) / value(T1) * partials(T1) +
        gas.R / (ηp * value(PR)) * partials(PR)
    Dual{Tag}(T2, T2 / cp(gas, T2) * ∂)
end

end # module
