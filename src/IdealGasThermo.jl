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
export FrozenGas, props, temperature
include("fastfrozengas.jl")
export FastFrozenGas
include("gasstate.jl")
export GasState, entropy, density
export compress, expand, expand_to, add_heat, add_work, extract_work
include("combustor.jl")
export Combustor, products
include("mixer.jl")
export Mixer, mixed
include("atmosphere.jl")

const DryAir = let g = Gas()
    g.X = Xair
    generate_composite_species(g.X, "Dry Air")
end
export DryAir
include("humidity.jl")
export humid_air

end
