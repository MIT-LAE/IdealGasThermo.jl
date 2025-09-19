# """
# Thermally-perfect gas thermodynamics based on NASA polynomials
# """
module IdealGasThermo

const __Gasroot__ = dirname(@__DIR__)
const default_thermo_path = joinpath(__Gasroot__, "data/thermo.inp")
using LinearAlgebra
using StaticArrays
using Printf
using ForwardDiff

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
include("atmosphere.jl")

gas = Gas()
gas.X = Xair
const DryAir = generate_composite_species(gas.X, "Dry Air")
export DryAir
include("humidity.jl")
end
