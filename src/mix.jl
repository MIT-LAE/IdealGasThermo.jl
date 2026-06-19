# Energy-conserving, composition-carrying mixing of two gases. Replaces the
# Mixer/mixed pair: because a FrozenGas now carries its composition (`gas.X`),
# two gases mix directly — no precomputed mixing system is needed. See ADR-0008.

"""
    mix(a::FrozenGas, b::FrozenGas, mratio) -> FrozenGas

Composition of the gas formed by blending `mratio = mass_b / mass_a` parts of
`b` into one part of `a`. Pure, zero-allocation, smooth in `mratio` (and in the
inputs' own composition, for the Dual-carrying case). `mix(a, b, 0)` is `a`; `b`
is the `mratio → ∞` limit.

The merged mole fractions are `X = (n_a·a.X + n_b·b.X)/(n_a + n_b)` with molar
amounts `n_a = 1/a.MW`, `n_b = mratio/b.MW` (the mass ratio expressed per mole),
and the merged [`FrozenGas`](@ref) is rebuilt from `X` — identical to building
any gas from that composition, entropy of mixing included.
"""
function mix(a::FrozenGas, b::FrozenGas, mratio)
    na = 1 / a.MW
    nb = mratio / b.MW
    X = (na * a.X + nb * b.X) / (na + nb)
    FrozenGas(X)
end

"""
    mix(a::GasState, b::GasState, mratio) -> GasState

Adiabatic mix of two streams at mass ratio `mratio = ṁ_b/ṁ_a`: the merged
composition (as above) plus the energy balance — the mixed temperature is the
mass-averaged total enthalpy `h = (h_a + mratio·h_b)/(1 + mratio)` inverted on
the merged gas. The streams must be at the same pressure (an isobaric mixer);
a non-isobaric mix needs a momentum closure, which belongs to the flow layer
(ADR-0005), not here. Returns the merged stagnation state.
"""
function mix(a::GasState, b::GasState, mratio)
    a.P ≈ b.P || throw(
        ArgumentError(
            "mix requires equal stream pressures (got $(a.P) and $(b.P)); a " *
            "non-isobaric mix needs a momentum closure, which lives in the flow " *
            "layer, not IdealGasThermo.",
        ),
    )
    gm = mix(a.gas, b.gas, mratio)
    hmix = (h(a.gas, a.T) + mratio * h(b.gas, b.T)) / (1 + mratio)
    GasState(gm, T_of_h(gm, hmix), a.P)
end
