# [Humid air and saturation](@id tutorial-humid-air)

!!! note "Scope: frozen composition, no phase change"
    `IdealGasThermo.jl` models **frozen-composition ideal gases** — there is no
    condensation, latent heat, or vapour–liquid equilibrium. "Humid air" here is
    dry air plus a *fixed* amount of water **vapour**; the saturation-pressure
    model only tells you *when* that vapour would begin to condense (the dew
    point), not what happens afterwards.

## Moist air as a `FrozenGas`

[`humid_air`](@ref) builds dry air plus water vapour, specified by **either** the
specific humidity ``\omega`` (`SH`) **or** the relative humidity (`RH`):

```@example humid
using IdealGasThermo

wet = humid_air(RH = 0.60, T = 300.0, P = 101325.0)   # 60% RH at 27 °C
dry = FrozenGas(DryAir)

(R_wet = R(wet), R_dry = R(dry))
```

Water (``M \approx 18`` g/mol) is lighter than dry air (``\approx 29`` g/mol), so
adding vapour *raises* the specific gas constant ``R = \Ru/M``. The
result is an ordinary `FrozenGas`, so all the usual accessors apply:

```@example humid
(cp = c_p(wet, 300.0), h = h(wet, 300.0))
```

## Saturation and the dew point

The water saturation pressure ``P_\text{sat}(T)`` sets the maximum vapour the air
can hold; relative humidity is ``\mathrm{RH} = p_v / P_\text{sat}(T)``:

```@example humid
T = 300.0
Psat = IdealGasThermo.saturation_vapor_pressure(T)   # Pa
```

The specific humidity at saturation (``\mathrm{RH} = 1``) at this temperature and
pressure is the most vapour the parcel could carry before condensing:

```@example humid
ω_sat = IdealGasThermo.specific_humidity(1.0, T, 101325.0)
```

Cooling the 60%-RH parcel raises its relative humidity (``P_\text{sat}`` falls
while ``p_v`` is fixed) until it reaches ``\mathrm{RH} = 1`` — the **dew point**.
Below it, the frozen-composition model no longer applies, because water would
leave the gas phase.
