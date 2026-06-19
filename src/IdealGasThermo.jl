# """
# Thermally-perfect gas thermodynamics based on NASA polynomials
# """
module IdealGasThermo

const __Gasroot__ = dirname(@__DIR__)
const default_thermo_path = joinpath(__Gasroot__, "data/thermo.inp")
using LinearAlgebra
using StaticArrays
using Printf

export Gas, set_h!, set_hP!, set_TP!, set_Δh!

include("constants.jl")
include("deprecation.jl")
include("species.jl")
export AbstractSpecies, species, composite_species, generate_composite_species

include("readThermo.jl")
export readThermo, species_in_spdict
include("Gas.jl")
include("Gas1D.jl")
export Gas1D
include("combustion.jl")
include("turbo.jl")
include("io.jl")
export print_thermo_table
include("utils.jl")
export X2Y, Y2X
include("thermoProps.jl")
include("frozengas.jl")
export FrozenGas, props
# Pure-core property accessors. The standalone `cp` function collides with
# `Base.cp` (file copy), so it is deliberately NOT exported; `cₚ` and `c_p` are
# exported aliases of the same function (same methods, including the
# ForwardDiff-extension Dual methods). `cp` stays reachable as `IdealGasThermo.cp`,
# and the `props(gas, T)` NamedTuple field keeps the name `cp` — a field accessor
# is reached only as `.cp`, never as a bare identifier, so it cannot collide with
# `Base.cp`. (ADR-0007; see CHANGELOG migration notes.)
const cₚ = cp
const c_p = cp
# `γ` is an exported Unicode alias of `gamma` (ratio of specific heats); unlike
# `cp`, the ASCII `gamma` does not collide with Base, so both names are exported.
const γ = gamma
export cₚ, c_p, h, s0, gamma, γ, R, T_from_h, pressure_ratio
include("fastfrozengas.jl")
export FastFrozenGas
include("gasstate.jl")
export GasState, entropy, density
export compress, expand, expand_to, add_heat, add_work, extract_work
include("flow.jl")
export speed_of_sound, mach, stagnation_state, static_state
include("vitiator.jl")
export Vitiator, products
include("mix.jl")
export mix
include("atmosphere.jl")

# Dry-air pseudo-species, built directly from the `Xair` mole-fraction table.
# Deliberately does NOT go through `Gas()` (deprecated, ADR-0002/0007): the pure
# core must stand entirely free of the mutable legacy layer so that layer can be
# deleted in v2.0.0 without touching any live code path. `Xidict2Array` places the
# mole fractions directly (no lossy X→Y→X round-trip), giving a composite identical
# to the old `let g=Gas(); g.X=Xair; … end` form to machine precision (MW bit-equal).
const DryAir = generate_composite_species(Xidict2Array(Xair), "Dry Air")
export DryAir
include("humidity.jl")
export humid_air

end
