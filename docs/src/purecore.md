# Pure-core API (v2)

The immutable, allocation-free pure core. A substance is a set of property
curves ([`FrozenGas`](@ref), or [`FastFrozenGas`](@ref) with tabulated
inverses); a thermodynamic point pairs a substance with `T` and `P`
([`GasState`](@ref)); process and flow *verbs* map states to states. Combustion,
mixing and humidity all produce immutable gases. See the architecture decision
records under `docs/adr/` for the rationale.

!!! note "Specific heat is exported as `cₚ` / `c_p`"
    The specific-heat accessor is exported under the interchangeable aliases
    `cₚ` and `c_p` (the bare name `cp` would shadow `Base.cp`); the ratio of
    specific heats is exported as `gamma` and the Unicode alias `γ`.

## Substances and properties

```@autodocs
Modules = [IdealGasThermo]
Pages = ["frozengas.jl", "fastfrozengas.jl"]
```

## Thermodynamic state and process verbs

```@autodocs
Modules = [IdealGasThermo]
Pages = ["gasstate.jl"]
```

## Gas dynamics

```@autodocs
Modules = [IdealGasThermo]
Pages = ["flow.jl"]
```

## Combustion, mixing and humidity

```@autodocs
Modules = [IdealGasThermo]
Pages = ["combustor.jl", "mixer.jl", "humidity.jl"]
```
