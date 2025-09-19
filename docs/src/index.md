# IdealGasThermo.jl Documentation

`IdealGasThermo.jl` is a tool for fast thermodynamic properties calculations for ideal gases. This assumes that the specifc heat
of the gas/mixture is only a function of its temperature, i.e., $c_p(T)$ , 
$h(T)$, and $s(T,p)$ (note the entropy is a function of both pressure and temperature).

One of the important features of `IdealGasThermo.jl` is the ability to represent
a mixture of gases as a single [`composite_species`](@ref vitiated). This package uses parametric typing to make it compatible with forward-mode automatic differentiation.

## Getting started

There are several workflows that are possible to use `IdealGasThermo.jl`. We outline here the most common few.

### Simple install

The easiest way to run `IdealGasThermo.jl` would be to add the package using the julia package manager using the github repository.

You can do this by starting a Julia session and then activating the package manager by typing `]` and then entering:
```julia-repl
pkg> add "https://github.com/MIT-LAE/IdealGasThermo.jl.git"
```

You can then import `IdealGasThermo` as you would with any Julia package:
```julia-repl
julia> using IdealGasThermo
```
### Local development

If you are going to develop the source code of `IdealGasThermo.jl` you might benefit from a local clone of the git repository which
can then fit into a workflow using [`Revise.jl`](https://timholy.github.io/Revise.jl/stable/) for example.

Step 1: Clone the git repo locally
```bash
git clone https://github.com/MIT-LAE/IdealGasThermo.jl.git
```

Step 2: `cd` to the folder where IdealGasThermo is cloned

Step 3: Use `Pkg` to install/ develop the package

```julia
pkg> dev .
```

You should now be able to import IdealGasThermo from within any Julia script in your base environment.

!!! tip "Tip"
    If you are using `Revise.jl` be sure to first import `Revise` before importing `IdealGasThermo`

    ```julia
    using Revise
    using IdealGasThermo
    ```
