# [Combustion and mixing](@id tutorial-combustion)

An executable walkthrough of ideal combustion ([`Vitiator`](@ref) /
[`products`](@ref)) and stream mixing ([`mix`](@ref)). As before, every block is
run when the docs are built.

## Vitiated products

A `Vitiator` is the *composition* model of combustion: build it **once** for a
given fuel and oxidizer, then call `products(vit, FAR)` for any fuel–air ratio
``\mathrm{FAR}`` (by mass). It returns an immutable `FrozenGas` — the burned-gas
substance — with no allocation. (We choose `Vitiator` because something like `Combustor` is 
very likely to conflict with the use cases of this package)

```@example comb
using IdealGasThermo

vit    = Vitiator("CH4", DryAir)     # methane burned in dry air, built once
burned = products(vit, 0.03)         # FAR = 0.03 (lean)

(R = R(burned), cp_1500 = c_p(burned, 1500.0))
```

The burned gas is a `FrozenGas` like any other, so every property and process
verb applies to it. For example, its enthalpy rise from the unburned reference:

```@example comb
h(burned, 1500.0) - h(burned, 298.15)
```

Because `products` is a smooth, allocation-free function of ``\mathrm{FAR}``, you
can differentiate through it — see [Thermodynamic derivatives](@ref derivatives).

## Mixing two streams

`mix(a, b, mratio)` blends two gases at mass ratio
``\mathrm{mratio} = \dot m_b / \dot m_a``. On bare `FrozenGas`es it blends the
**composition** (each gas remembers its mole fractions):

```@example comb
core   = products(vit, 0.03)         # hot combustion products
bypass = FrozenGas(DryAir)           # cold bypass air

blend = mix(core, bypass, 5.0)       # 5 parts bypass per part core
c_p(blend, 1000.0)
```

On [`GasState`](@ref)s, `mix` additionally closes the **energy balance** — the
mixed temperature is the mass-averaged total enthalpy

```math
h_\text{mix} = \frac{h_a + \mathrm{mratio}\cdot h_b}{1 + \mathrm{mratio}},
```

inverted on the merged gas. It requires equal stream pressures (an isobaric
mixer):

```@example comb
P    = 1.0e5
hot  = GasState(core,   1500.0, P)
cold = GasState(bypass,  800.0, P)

mixed = mix(hot, cold, 5.0)          # bypass ratio 5:1
(T_mixed = mixed.T, P_mixed = mixed.P)
```

The mixed temperature lands between the two streams, weighted toward the bypass 
flow. The merged gas keeps its blended composition, so it can feed
an afterburner `Vitiator` or a further `mix` — the whole chain stays immutable
and allocation-free.
