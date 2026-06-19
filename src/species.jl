
abstract type AbstractSpecies end

"""
    species <: AbstractSpecies
species is a structure that holds the NASA 9 polynomial coefficients `alow` and `ahigh` 
for the two temprature regions separated by `Tmid` 
(here we only work with temperature less than 6000 K so typically only 2 T intervals required)
the molecular weight `MW` and the heat of formation `Hf` (J/mol) for a given chemical species (at 298.15 K).

See [here](https://shepherd.caltech.edu/EDL/PublicResources/sdt/formats/nasa.html) for typical data format
"""
struct species <: AbstractSpecies
    name::String
    Tmid::Float64
    alow::Array{Float64,1}
    ahigh::Array{Float64,1}
    MW::Float64
    Hf::Float64
    formula::AbstractString
end

"""
    composite_species <: AbstractSpecies

Represents a gas mixture of multiple components as a 
psuedo-species by calculating an equivalent set of polynomials 
defining ``c_p``, ``h``, and ``s``.  

See [here](@ref gas1dthermo) for a more detailed explanation.
"""
struct composite_species <: AbstractSpecies
    name::String
    Tmid::Float64
    alow::Array{Float64,1}
    ahigh::Array{Float64,1}
    MW::Float64
    Hf::Float64
    composition::Dict{String,Float64}
end
"""
    generate_composite_species(Xi::AbstractVector, name::AbstractString="composite species")

Generates a composite psuedo-species to represent a gas mixture given the
mole fraction `Xi` of its constitutents.
"""
function generate_composite_species(
    Xi::AbstractVector,
    name::AbstractString = "composite species",
)
    if !(sum(Xi) ≈ 1.0)
        error("Gas mixture composition is not well defined. Sum of Xi = $(sum(Xi)) != 1.0")
    end
    if any(Xi .< 0.0)
        error("Composition has negative values.")
    end
    # Equivalent molar coefficients (entropy of mixing folded into b₂) via the
    # shared lumping kernel — same `A * X` over the species basis used by the
    # pure-core constructors. Tmid is always 1000 K (checked at read time).
    alow, ahigh, MW, Hf = _lump_molar(Xi)
    d = Dict{String,Float64}()
    for i in eachindex(Xi)
        if Xi[i] != 0
            push!(d, spdict.name[i] => Xi[i])
        end
    end
    return composite_species(name, 1000.0, Vector(alow), Vector(ahigh), MW, Hf, d)
end  # function generate_composite_species
