# IdealGasThermo.jl

[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://mit-lae.github.io/IdealGasThermo.jl/dev/) [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://mit-lae.github.io/IdealGasThermo.jl/stable/) ![Codecov](https://img.shields.io/codecov/c/github/MIT-LAE/IdealGasThermo.jl) [![CI](https://github.com/MIT-LAE/IdealGasThermo.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/MIT-LAE/IdealGasThermo.jl/actions/workflows/CI.yml)


A simple package to efficiently calculate ideal gas properties based on NASA polynomials. This assumes that the specifc heat
of the gas/mixture is only a function of its temperature, i.e., $c_p(T)$ , $h(T)$, and $s(T,p)$ (note the entropy is a function of both pressure and temperature).

This package is primarily intended for use in gas turbine cycle calculations with a focus on computational efficiency (and gradient calculations for downstream calculations). As such, there are a number of functions here that are relevant for thermodynamic cycle analysis (e.g., compressing a gas at a given polytropic efficiency, setting a flow Mach number, mixing, and ideal combustion etc). 

This package also introduces a method for fast calculations of thermodynamic properties of mixtures with fixed composition by calculating an equivalent polynomial representation of the mixture. See docs for theory and performance details. 

## Installation

The package is not registered in the Julia General registry, so install it directly from the repository.

```julia
using Pkg
Pkg.add(url = "https://github.com/MIT-LAE/IdealGasThermo.jl")
```

`using IdealGasThermo` to load it (the module name is `IdealGasThermo`, independent of the repository name).

### v2.0 pre-release (beta)

The `2.0.0-beta` series promotes the immutable pure core (`FrozenGas`, `GasState`, and the process/flow verbs) to the headline API and **loudly deprecates** the legacy mutable layer (`Gas`, `Gas1D`, and the functions they back), which is removed in `2.0.0` final. To pin a downstream project to a specific beta during migration, install by tag:

```julia
using Pkg
Pkg.add(url = "https://github.com/MIT-LAE/IdealGasThermo.jl", rev = "v2.0.0-beta1")
```

See [`CHANGELOG.md`](CHANGELOG.md) for the migration mapping and [`docs/adr/0007-v2-deprecated-beta-migration.md`](docs/adr/0007-v2-deprecated-beta-migration.md) for the rationale.

 
