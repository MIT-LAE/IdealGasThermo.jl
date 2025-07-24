
abstract type AbstractSpecies end

"""
    species <: AbstractSpecies
species is a structure that holds the NASA 9 polynomial coefficients `alow` and `ahigh` 
for the two temprature regions separated by `Tmid` 
(here we only work with temperature less than 6000 K so typically only 2 T intervals required)
the molecular weight `MW` and the heat of formation `Hf` (J/mol) for a given chemical species (at 298.15 K).

See [here](https://shepherd.caltech.edu/EDL/PublicResources/sdt/formats/nasa.html) for typical data format
"""
struct species{R<:Real} <: AbstractSpecies
    name::String
    Tmid::R
    alow::Array{R,1}
    ahigh::Array{R,1}
    MW::R
    Hf::R
    formula::AbstractString
end

"""
    composite_species <: AbstractSpecies

Represents a gas mixture of multiple components as a 
psuedo-species by calculating an equivalent set of polynomials 
defining ``c_p``, ``h``, and ``s``.  

See [here](@ref gas1dthermo) for a more detailed explanation.
"""
struct composite_species{R<:Real} <: AbstractSpecies
    name::String
    Tmid::R
    alow::Array{R,1}
    ahigh::Array{R,1}
    MW::R
    Hf::R
    composition::Dict{String,R}
end
"""
    generate_composite_species(Xi::AbstractVector, name::AbstractString="composite species")

Generates a composite psuedo-species to represent a gas mixture given the
mole fraction `Xi` of its constitutents.
"""
function generate_composite_species(
    Xi::AbstractVector{R},
    name::AbstractString = "composite species",
) where R<:Real
    ALOW = reduce(hcat, spdict.alow)
    AHIGH = reduce(hcat, spdict.ahigh)
    alow = ALOW * Xi
    ahigh = AHIGH * Xi
    MW = dot(spdict.MW, Xi)
    Hf = dot(spdict.Hf, Xi)
    # Need to account for the entropy of mixing:
    Δs_mix = 0.0
    for i in eachindex(Xi)
        if Xi[i] != 0.0
            Δs_mix = Δs_mix + Xi[i] * log(Xi[i])
        end
    end
    # This is independent of temperature represented by the constant of integration:
    alow[end] = alow[end] - Δs_mix
    ahigh[end] = ahigh[end] - Δs_mix
    Tmid = 1000.0 #Assumed to always be the mid (this is checked for all input thermo when read)
    # comp
    d = Dict{String, R}()
    if !(sum(Xi) ≈ 1.0)
        error("Gas mixture composition is not well defined. Sum of Xi = $(sum(Xi)) != 1.0")
    end
    if any(Xi .< 0.0)
        error("Composition has negative values.")
    end
    for i in eachindex(Xi)
        if (Xi[i] != 0)
            push!(d, spdict.name[i] => Xi[i])
        end
    end
    return composite_species(name, Tmid, alow, ahigh, MW, Hf, d)

end  # function generate_composite_species
