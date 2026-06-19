# Representing ideal gases

## Pure (single-component) gases
`IdealGasThermo.jl` exports the `species` type which stores relevant thermodynamic
information about a single element/compound. See [`readThermo`](@ref).

```@docs
species
```

## Composite species

```@docs
composite_species
```
!!! warning "Deprecated — use `FrozenGas`"
    The mutable `Gas` / `Gas1D` types and their in-place `set_*!` mutators
    documented below are **deprecated** as of `v1.1.0` and will be removed in a
    future `v2.0.0`. New code should use the immutable [`FrozenGas`](@ref) pure
    core — see [The immutable pure core](@ref frozengas-page). (The `species` and
    `composite_species` types above are *not* deprecated.)

## Gas mixtures
`IdealGasThermo.jl` exports the `Gas` type which stores relevant thermodynamic
information about the gas mixture.

```@docs
Gas
Gas()
Base.setproperty!(gas::Gas, sym::Symbol, val::Float64)
```

## Single component gases

`Gas1D` type objects are a subtype of `AbstractGas` which allows us to use most of the functions that work with [`Gas`](@ref). `Gas1D` types additionally store a representation of the composite species ([`composite_species`](@ref)). See [here](@ref gas1dthermo) for the theory of representing fixed composition multi-component mixtures as single component mixtures.

```@docs
Gas1D
Gas1D()
Gas1D(sp::composite_species)
```
## Setting the thermodynamic state of the gas

The following functions let you set the thermodynamic state of the gas. 
These functions change the state of the gas *in place* i.e., the gas object
is modified and no new copy is created.

```@docs
set_h!
set_hP!
set_TP!
set_Δh!
```