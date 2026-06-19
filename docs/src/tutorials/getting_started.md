# [Getting started: a gas and a process chain](@id tutorial-getting-started)

This tutorial creates an immutable [`FrozenGas`](@ref), reads its thermodynamic
properties, and pushes it through a small turbomachinery process chain. **Every
code block below is executed when the documentation is built**, so the printed
outputs are real.

A `FrozenGas` is a *substance*: a fixed-composition set of property curves
``c_p(T)``, ``h(T)``, ``s^0(T)``. Temperature is an **argument**, instead of a field
stored on the gas — see [Substance vs. state](@ref frozengas-substance-state). This lets 
us keep most of the structures immutable that has significant performance benefits. 

## Create a gas

```@example gs
using IdealGasThermo

air = FrozenGas(DryAir)
```

`DryAir` is a built-in composite of ``\mathrm{N_2}/\mathrm{O_2}/\mathrm{Ar}/\mathrm{CO_2}``.
A gas can equally be built from a single database species or from a
mole-fraction vector:

```@example gs
n2 = FrozenGas(species_in_spdict("N2"))
nothing # hide
```

## Read properties

The property accessors are pure functions of ``(\mathrm{gas}, T)``. At
``T = 600`` K:

```@example gs
T = 600.0   # K
(cp = c_p(air, T), h = h(air, T), s0 = s0(air, T), γ = γ(air, T))
```

`props` returns ``c_p``, ``h`` and ``s^0`` together (sharing the temperature
powers and the single ``\log T``), which is what you want in a hot loop:

```@example gs
props(air, T)
```

The specific gas constant ``R``, ratio of specific heats ``\gamma`` and speed of
sound ``a = \sqrt{\gamma R T}`` are also one call away:

```@example gs
(R = R(air), γ300 = γ(air, 300.0), a300 = speed_of_sound(air, 300.0))
```

## Invert: temperature from enthalpy

`T_from_h` inverts ``h(\mathrm{gas}, T)`` — given an enthalpy it returns the
temperature, by a bounded Newton solve:

```@example gs
h_target = h(air, 900.0)
T_from_h(air, h_target)   # ≈ 900.0
```

## Do something: a process chain

A [`GasState`](@ref) pairs a substance with a ``(T, P)`` point. Process *verbs*
map a state to a **new** state (the substance is never mutated), so a cycle reads
as the chain it is. Here is a single-spool core: compress, burn (as a
constant-pressure heat addition), then extract the compressor work back out
through the turbine.

```@example gs
inlet = GasState(air, 288.15, 101325.0)          # 15 °C, 1 atm

comp  = compress(inlet, 12.0; ηp = 0.90)         # 12:1, polytropic η = 0.90
burn  = add_heat(comp, 8.0e5)                     # add 800 kJ/kg at constant P
turb  = extract_work(burn, h(comp.gas, comp.T) - h(inlet.gas, inlet.T))

(T_comp = comp.T, T_burn = burn.T, T_turb = turb.T, P_turb = turb.P)
```

Two derived quantities need the full ``(T, P)`` pair — entropy and density:

```@example gs
(s_inlet = entropy(inlet), ρ_inlet = density(inlet),
 s_burn  = entropy(burn),  ρ_burn  = density(burn))
```

Each station above is an immutable value; collecting them is just keeping the
bindings. Because the pressure travels *inside* the state, the "advanced the
wrong pressure rail" class of bug cannot be written.

## Where to go next

- [Combustion and mixing](@ref tutorial-combustion) — `Vitiator`/`products` and `mix`.
- [Humid air and saturation](@ref tutorial-humid-air).
- [Thermodynamic derivatives](@ref derivatives) — differentiate any of the above
  with ForwardDiff, analytically.
- [The immutable pure core](@ref frozengas-page) — why immutability buys all of this.
