
"""
    fuelbreakdown(fuel::String)

Returns the number of C, H, O, and N atoms that the fuel is composed of.

# Examples
```julia-repl

julia> IdealGasThermo.fuelbreakdown("CH4")' #transpose is simply to save space in the docs
1×4 adjoint(::Vector{Float64}) with eltype Float64:
 1.0  4.0  0.0  0.0

julia> IdealGasThermo.fuelbreakdown("C12H23.5")'
1×4 adjoint(::Vector{Float64}) with eltype Float64:
 12.0  23.5  0.0  0.0

julia> IdealGasThermo.fuelbreakdown("CH3COOH")'
1×4 adjoint(::Vector{Float64}) with eltype Float64:
 2.0  4.0  2.0  0.0

julia> IdealGasThermo.fuelbreakdown("CH3CH2OH")'
1×4 adjoint(::Vector{Float64}) with eltype Float64:
 2.0  6.0  1.0  0.0

```
"""
function fuelbreakdown(fuel::String)
    C, H, O, N = 0.0, 0.0, 0.0, 0.0
    if !isempty(findall(r"[^cChHoOnN.^[0-9]", fuel))
        try
            fuel = species_in_spdict(fuel).formula
        catch e
            if isa(e, ArgumentError)
                error("""The input fuel string $fuel is not found in 
                the thermo database and contains
                elements other than C,H,O, and N.\n""")
            end
        end
    end
    chunks = [fuel[idx] for idx in findall(r"[a-zA-Z][a-z]?\d*\.?\d*", fuel)]
    for chunk in chunks
        element, number = match(r"([a-zA-Z][a-z]?)(\d*\.?\d*)", chunk).captures
        element = uppercase(element)
        if isempty(number)
            number = 1
        else
            number = parse(Float64, number)
        end
        if element == "C"
            C = C + number
        elseif element == "H"
            H = H + number
        elseif element == "O"
            O = O + number
        elseif element == "N"
            N = N + number
        else
            error("Fuel can only contain C, H, O or N atoms!")
        end
    end
    return ([C, H, O, N])

end

fuelbreakdown(fuel::species) = fuelbreakdown(fuel.formula)

"""
    reaction_change_fraction(fuel::String)

Returns the mass fraction change due to complete combustion

Assume fuel of type  CᵢHⱼOₖNₗ , then
```
  CᵢHⱼOₖNₗ + n(O2) * O2 ---> n(CO2)*CO2 + n(H2O)*H2O + n(N2)*N2
⟹CᵢHⱼOₖNₗ              ---> n(CO2)*CO2 + n(H2O)*H2O + n(N2)*N2 - n(O2)*O2 
```
# Examples
```julia-repl
julia> reaction_change_fraction("CH4")
Dict{String, Float64} with 4 entries:
  "O2"  => -3.98926
  "H2O" => 2.24595
  "CO2" => 2.74331
  "N2"  => 0.0
```
"""
function reaction_change_fraction(fuel::String)
    CHON = fuelbreakdown(fuel) # Returns number of C, H, O, and N atoms in fuel
    element_mass_array = [12.01070, 1.00794, 15.99940, 14.00670]

    Mfuel = dot(CHON, element_mass_array)

    M_CO2 = 44.00950
    M_N2 = 28.01340
    M_H2O = 18.01528
    M_O2 = 31.99880

    # Calculate the number of moles of products + O₂
    nCO2 = CHON[1] * 1.0
    nN2 = CHON[4] / 2
    nH2O = CHON[2] / 2
    nO2 = -(CHON[1] + CHON[2] / 4 - CHON[3] / 2) # Oxygen is used up/ lost

    ΔY = [M_CO2, M_N2, M_H2O, M_O2] .* [nCO2, nN2, nH2O, nO2] ./ Mfuel

    # return as dict to make it easier to set Gas mass fractions
    names = ["CO2", "N2", "H2O", "O2"]
    Ydict = Dict(zip(names, ΔY))

    return Ydict

end

"""
    reaction_change_molar_fraction(fuel::AbstractString)

Returns the mole fraction change due to complete combustion of one mole of
the specified fuel

Assume fuel of type  CᵢHⱼOₖNₗ , then
```
    CᵢHⱼOₖNₗ + n(O2) * O2 ---> n(CO2)*CO2 + n(H2O)*H2O + n(N2)*N2
   ⟹CᵢHⱼOₖNₗ              ---> n(CO2)*CO2 + n(H2O)*H2O + n(N2)*N2 - n(O2)*O2 
```

# Examples
```julia-repl
julia> IdealGasThermo.reaction_change_molar_fraction("CH4")
4-element Vector{Float64}:
  1.0
  0.0
  2.0
 -2.0
```
"""
function reaction_change_molar_fraction(fuel::AbstractString)
    CHON = fuelbreakdown(fuel) # Returns number of C, H, O, and N atoms in fuel
    # Calculate the number of moles of products + O₂
    nCO2 = CHON[1] * 1.0
    nN2 = CHON[4] / 2
    nH2O = CHON[2] / 2
    nO2 = -(CHON[1] + CHON[2] / 4 - CHON[3] / 2) # Oxygen is used up/ lost

    X = [nCO2, nN2, nH2O, nO2]
    return X
end  # function reaction_change_molar_fraction

"""
    stoich_molar_fuel_oxy_ratio(fuel::AbstractString)

Calculates the molar fuel-**oxygen** ratio for stoichiometric combustion.

# Examples
```julia-repl
julia> using IdealGasThermo

julia> IdealGasThermo.stoich_molar_fuel_oxy_ratio("CH4")
0.5

julia> IdealGasThermo.stoich_molar_fuel_oxy_ratio("C12H23")
0.056338028169014086
```
"""
function stoich_molar_fuel_oxy_ratio(fuel::AbstractString)
    X = reaction_change_molar_fraction(fuel)
    molFuelOxyRatio = abs(1 / X[end])
    return molFuelOxyRatio
end  # function stoich_molar_FAR

"""
    stoich_molar_FOR(fuel::AbstractSpecies, oxidizer::AbstractSpecies)

Calculates the **molar** fuel-oxidizer ratio for stoichiometeric combustion for 
and arbitrary fuel and oxidizer.

# Examples
```julia-repl
julia> CH4 = species_in_spdict("CH4");

julia> IdealGasThermo.stoich_molar_FOR(CH4)
0.104738
```
"""
function stoich_molar_FOR(fuel::AbstractSpecies, oxidizer::AbstractSpecies = DryAir)
    molFuelOxyRatio = stoich_molar_fuel_oxy_ratio(fuel.name)
    if oxidizer isa species
        Xin = Dict(oxidizer.name => 1.0)
        if oxidizer.name == "Air"
            Xin = Xair
        end
    else
        Xin = oxidizer.composition
    end

    return molFuelOxyRatio * Xin["O2"]

end  # function stoich_molar_FOR 

"""
    stoich_FOR(fuel::AbstractSpecies, oxidizer::AbstractSpecies=DryAir)

Calculates the **mass** fuel-oxidizer ratio for stoichiometeric combustion for 
and arbitrary fuel and oxidizer.

# Examples
```julia-repl
julia> IdealGasThermo.stoich_FOR(CH4)
0.05800961333050494
```
"""
function stoich_FOR(fuel::AbstractSpecies, oxidizer::AbstractSpecies = DryAir)
    molFOR = stoich_molar_FOR(fuel, oxidizer)
    return molFOR * fuel.MW / oxidizer.MW
end  # function stoich_FOR

stoich_FOR(fuel::AbstractString, oxi::AbstractString) =
    stoich_FOR(species_in_spdict(fuel), species_in_spdict(oxi))

"""
    vitiated_mixture(fuel, oxidizer, 
    FAR, ηburn=1.0)

Calculates the composition of a burnt gas mixture. Defaults to stoichiometric 
conditions if FAR is not specified. `vitiated_mixture` returns the number of 
moles of each species present in the burnt gas mixture after combustion at the 
specified FAR. *Note* the sum of result in general will **not** sum to 1. 

# Examples
```julia-repl
julia> CH4 = species_in_spdict("CH4");

julia> Air = IdealGasThermo.DryAir;

julia> IdealGasThermo.vitiated_mixture(CH4, Air, 0.04)
Dict{Any, Any} with 6 entries:
  "O2"  => 0.0650337
  "CH4" => 0.0
  "Ar"  => 0.009365
  "H2O" => 0.144442
  "CO2" => 0.0725401
  "N2"  => 0.78084

julia> IdealGasThermo.vitiated_mixture(CH4, Air)
Dict{Any, Any} with 6 entries:
  "O2"  => 0.0
  "CH4" => 0.0
  "Ar"  => 0.009365
  "H2O" => 0.209476
  "CO2" => 0.105057
  "N2"  => 0.78084
```
See [here](@ref vitiated) for some explanation of the background.
"""
function vitiated_mixture(
    fuel::AbstractSpecies,
    oxidizer::AbstractSpecies,
    FAR::R1,
    ηburn::R2 = 1.0,
) where {R1 <: Real, R2 <: Real}

    if oxidizer isa species
        Xin = Dict(oxidizer.name => 1.0)
        if oxidizer.name == "Air"
            Xin = Xair
        end
    else
        Xin = oxidizer.composition
    end

    molFAR = FAR * oxidizer.MW / fuel.MW

    nCO2, nN2, nH2O, nO2 = ηburn .* molFAR .* reaction_change_molar_fraction(fuel.name)

    nFuel = molFAR * (1.0 - ηburn)

    names = [fuel.name, "CO2", "H2O", "N2", "O2"]
    ΔX = [nFuel, nCO2, nH2O, nN2, nO2]

    Xdict = Dict(zip(names, ΔX))

    Xvitiated::Dict{String, R1} = mergewith(+, Xin, Xdict)

    for key in keys(Xvitiated) #Normalize such that sum(Xi) = 1
        Xvitiated[key] = Xvitiated[key] / (1 + ηburn * molFAR)
    end

    return Xvitiated
end  # function vitiated_gas

vitiated_mixture(fuel, oxidizer) =
    vitiated_mixture(fuel, oxidizer, stoich_FOR(fuel, oxidizer))

"""
    vitiated_mixture(fuel, oxidizer, 
    FAR, ηburn=1.0)

Convenience function that finds fuel and oxidizer from thermo database
"""
function vitiated_mixture(
    fuel::AbstractString,
    oxidizer::AbstractString,
    FAR::R1,
    ηburn::R2 = 1.0,
) where {R1 <: Real, R2 <: Real}

    fuel = species_in_spdict(fuel)
    oxidizer = species_in_spdict(oxidizer)

    return vitiated_mixture(fuel, oxidizer, FAR, ηburn)

end  # function vitiated_mixture

"""
    LHV(fuel::AbstractSpecies)

Returns the lower heating value (LHV) of the given fuel at standard conditions.

## Examples
```julia-repl
julia> ethanol = species_in_spdict("C2H5OH")
Species "C2H5OH"
MW = 46.06844 g/mol


julia> IdealGasThermo.LHV(ethanol)/1e6 # convert to MJ/kg
27.731522925456126

julia> CH4 = species_in_spdict("CH4");

julia> Jet = species_in_spdict("Jet-A(g)");

julia> fuels = [CH4, Jet, ethanol];

julia> IdealGasThermo.LHV.(fuels)./1e6 # Can brodcast across fuels
3-element Vector{Float64}:
 50.027364880448516
 43.35316346765443
 27.731522925456126
```

"""
function LHV(fuel::AbstractSpecies)
    nCO2, _, nH2O, _ = reaction_change_molar_fraction(fuel.name)

    HfCO2 = -393510.0
    HfH2O = -241826.0

    ΔHrxn = (nCO2 * HfCO2 + nH2O * HfH2O) - fuel.Hf

    return -ΔHrxn / fuel.MW * 1000.0 #J/kg

end  # function LHV

LHV(fuel::AbstractString) = LHV(species_in_spdict(fuel))

"""
    AFT(fuel::AbstractSpecies, oxidizer::AbstractSpecies=DryAir)

Calculates the adiabatic flame temperature at standard conditions.
"""
function AFT(fuel::AbstractSpecies, oxidizer::AbstractSpecies = DryAir)
    molFAR = stoich_molar_FOR(fuel, oxidizer)
    gas = Gas()
    dict::Dict{String,Float64} = Dict()

    try
        dict = copy(oxidizer.composition)
        push!(dict, fuel.name => molFAR)
    catch
        dict = Dict(fuel.name => molFAR, oxidizer.name => 1.0)
    end

    gas.X = dict
    gas.T = 298.15
    h0 = gas.h
    Hfr = gas.Hf
    # Hr298 = gas.h*gas.MW/1000.0

    gas.X = vitiated_mixture(fuel, oxidizer)
    # Hp298 = gas.h*gas.MW/1000.0
    Hfp = gas.Hf

    #Energy Balance: Hrf + H(298.15) - Hr298.15 = Hpf + H(T) - Hp298.15
    #                Hrf - Hpf = ΔH(T)
    hf = Hfr - Hfp + h0
    # set_h!(gas, Hpsens)
    set_h!(gas, hf)

    return gas.T
end  # function AFT
"""
    vitiated_species(fuel, oxidizer, 
    FAR, ηburn=1.0, name="vitiated species")

Returns a [`composite_species`](@ref) that represents the burnt gas mixture
at the specified FAR. If no FAR is provided stoichiometeric conditions are
assumed. 

# Examples
```julia-repl
julia> IdealGasThermo.vitiated_species(CH4, Air, 0.05, name = "CH4-Air-0.05")
Composite Species: "CH4-Air-0.05"
MW = 27.89510190126262 g/mol
with composition:
 Species        Xᵢ
      O2   0.02653
      Ar   0.00859
     H2O   0.16560
     CO2   0.08309
      N2   0.71619
------------------
       Σ   1.00000
```
"""
function vitiated_species(
    fuel::AbstractSpecies,
    oxidizer::AbstractSpecies,
    FAR::R1;
    ηburn::R2 = 1.0,
    name::AbstractString = "vitiated species",
) where {R1 <: Real, R2 <: Real}

    Xdict = vitiated_mixture(fuel, oxidizer, FAR, ηburn)

    X = zeros(R1, Nspecies)
    Xidict2Array!(Xdict, X)

    return generate_composite_species(X, name)
end  # function vitiated_species

function vitiated_species(
    fuel::AbstractString,
    oxidizer::AbstractString,
    FAR::R1;
    ηburn::R2 = 1.0,
    name::AbstractString = "vitiated species",
) where {R1 <: Real, R2 <: Real}

    fuel = species_in_spdict(fuel)
    oxidizer = species_in_spdict(oxidizer)
    return vitiated_species(fuel, oxidizer, FAR; ηburn = ηburn, name = name)

end

vitiated_species(f::AbstractSpecies, o::AbstractSpecies) =
    vitiated_species(f, o, stoich_FOR(f, o))

vitiated_species(f::AbstractString, o::AbstractString) =
    vitiated_species(f, o, stoich_FOR(f, o))

"""
    fixed_fuel_vitiated_species(fuel, oxidizer, ηburn=1.0)

Returns a function `burntgas(FAR::Real)` that is specific to the fuel and oxidizer
combination provided. This gives a highly performant function that can 
simply be called at any given FAR for that specific fuel+oxidizer combo.

This is ~4x faster than doing 
```burntgas(FAR) = vitiated_species("CH4", "Air", FAR)```

## Examples
```julia-repl
julia> burntgas = IdealGasThermo.fixed_fuel_vitiated_species(CH4, Air)
(::IdealGasThermo.var"#burntgas#52"{species, composite_species, Vector{Float64}, Vector{Float64}, Float64}) (generic function with 1 method)

julia> burntgas(0.05)
Composite Species: "burntgas(CH4 + Dry Air; 0.05)"
MW = 27.89510190126262 g/mol
with composition:
 Species        Xᵢ
      O2   0.02653
      Ar   0.00859
     H2O   0.16560
     CO2   0.08309
      N2   0.71619
------------------
       Σ   1.00000

julia> gas1 = Gas1D(burntgas(0.02)) #Returns a Gas1D intialized with the burnt gas properties
Gas1D(burntgas(CH4 + Dry Air; 0.02); MW = 28.51473501878705 g/mol)
at T = 298.15 K; P = 101.325 kPa
```
"""
function fixed_fuel_vitiated_species(fuel, oxidizer, ηburn::R = 1.0) where R <: Real

    nCO2, nN2, nH2O, nO2 = ηburn .* reaction_change_molar_fraction(fuel.name)
    nFuel = (1.0 - ηburn)
    massratio = oxidizer.MW / fuel.MW
    if oxidizer isa species
        Xin::Dict{String,Real} = Dict(oxidizer.name => 1.0)
        if oxidizer.name == "Air"
            Xin = Xair
        end
    else
        Xin = oxidizer.composition
    end
    names::Vector{String} = [fuel.name, "CO2", "H2O", "N2", "O2"]
    ΔX::Vector{R} = [nFuel, nCO2, nH2O, nN2, nO2]
    Xdict::Dict{String,Real} = Dict(zip(names, ΔX))

    ΔX_array = zeros(R, Nspecies)
    Xin_array = zeros(R, Nspecies)
    allnames = view(spdict.name, :)

    for (key, value) in Xdict
        index::Int64 = findfirst(x -> x == key, allnames)
        ΔX_array[index] = value
    end

    for (key, value) in Xin
        index::Int64 = findfirst(x -> x == key, allnames)
        Xin_array[index] = value
    end

    """
       burntgas(FAR::Float64)
    Inner function that will be returned
    """
    function burntgas(FAR::Real)
        molFAR = FAR .* massratio
        X = ΔX_array .* molFAR
        @. X = (X + Xin_array) / (1 + molFAR)
        generate_composite_species(X, "burntgas($(fuel.name) + $(oxidizer.name); $FAR)")

    end  # function burntgas

    return burntgas
end  # function fixed_fuel_vitiated_species

"""
    fuel_combustion(gas_ox, fuel, Tf, FAR, 
    ηburn = 1.0, hvap = 0.0)

This function returns an `AbstractGas` with the combustion products at the combustor exit temperature
for a given fuel type, oxidizer gas at a given enthalpy, and FAR. It includes the combustion 
efficiency and the fuel enthalpy of vaporization as optional inputs. 
"""
function fuel_combustion(
    gas_ox::AbstractGas,
    fuel::String,
    Tf::R1,
    FAR::R2,
    ηburn::R3 = 1.0,
    hvap::R4 = 0.0
) where  {R1 <: Real, R2 <: Real, R3 <: Real, R4 <: Real}
    #Create variables corresponding to the oxidizer and fuel species and mixtures
    fuel_sps = species_in_spdict(fuel)

    #Find oxidizer composition
    if gas_ox isa Gas1D
        gas_sps = gas_ox.comp_sp
    else
        if "Air" in keys(gas_ox.Xdict)
            Xin = Xair
        else
            Xin = gas_ox.Xdict
        end
        gas_sps = generate_composite_species(IdealGasThermo.Xidict2Array(Xin))
    end

    #Find the vectors with the fuel mole and mass fractions
    Xfuel = Xidict2Array(Dict([(fuel, 1.0)])) #Mole fraction
    Yfuel = X2Y(Xfuel) #Mass fraction
    gas_fuel = Gas(Tf, gas_ox.P) #Create a fuel gas to calculate enthalpy
    gas_fuel.Y = Yfuel

    #Find product composition and mole and mass fractions
    Xprod_dict = vitiated_mixture(fuel_sps, gas_sps, FAR, ηburn)
    Xprod = Xidict2Array(Xprod_dict)
    Yprod = X2Y(Xprod)

    #Initialize output 
    gas_prod = Gas(gas_ox.T, gas_ox.P) #Initialize with oxidizer properties
    gas_prod.Y = Yprod #Set correct composition

    # Combustion enthalpy calculations
    #enthalpy before reaction = enthalpy after reaction
    h0 = (gas_ox.h + FAR * (gas_fuel.h - abs(hvap))) / (1 + FAR)

    #Set the products at the correct enthalpy
    set_h!(gas_prod, h0)
    return gas_prod
end

"""
    gas_burn(gas_ox, fuel, Tf, Tburn, 
    ηburn = 1.0, hvap = 0.0)

This function returns a tuple containing the FAR and an `AbstractGas` with the combustion 
products for a given fuel type, oxidizer gas at a given enthalpy, and desired combustor exit temperature. 
It includes the combustion efficiency and the fuel enthalpy of vaporization as optional inputs. 
"""
function gas_burn(
    gas_ox::AbstractGas,
    fuel::String,
    Tf::R1,
    Tburn::R2,
    ηburn::R3 = 1.0,
    hvap::R4 = 0.0
) where {R1 <: Real, R2 <: Real, R3 <: Real, R4 <: Real}

    #Create variables corresponding to the oxidizer and fuel species and mixtures
    fuel_sps = species_in_spdict(fuel)

    #Extract composite species with oxidizer gas composition
    if gas_ox isa Gas1D
        gas_sps = gas_ox.comp_sp
    else
        if "Air" in keys(gas_ox.Xdict)
            Xin = Xair
        else
            Xin = gas_ox.Xdict
        end
        gas_sps = generate_composite_species(IdealGasThermo.Xidict2Array(Xin))
    end

    #Find the vectors with the fuel mole and mass fractions
    Xfuel = Xidict2Array(Dict([(fuel, 1.0)])) #Mole fraction
    Yfuel = X2Y(Xfuel) #Mass fraction
    gas_fuel = Gas(Tf, gas_ox.P)
    gas_fuel.Y = Yfuel #Create a fuel gas to calculate enthalpy

    #Store enthalpies of oxidizer and fuel at original temperatures
    ho = gas_ox.h
    hf = gas_fuel.h

    #Find change in gas composition for a FAR of 1
    nCO2, nN2, nH2O, nO2 = reaction_change_molar_fraction(fuel_sps.name)

    names = ["CO2", "H2O", "N2", "O2"]
    ΔX = [nCO2, nH2O, nN2, nO2]
    Xdict = Dict(zip(names, ΔX))

    Xc = Xidict2Array(Xdict)
    Yc = X2Y(Xc) #mass fraction change in combustion for FAR = 1

    gas_c = Gas(Tburn, gas_ox.P) #Create a "virtual" gas with the changes in combustion, for enthalpy
    #calculations
    gas_c.Y = Yc

    hc = gas_c.h #Enthalpy change for FAR = 1

    gas_ox_burnt = deepcopy(gas_ox) #make a copy to avoid modifying the input
    set_TP!(gas_ox_burnt, Tburn, gas_ox.P)
    ha = gas_ox_burnt.h #Enthalpy of original oxidizer at final temperature

    set_TP!(gas_fuel, Tburn, gas_ox.P)
    hff = gas_fuel.h #Enthalpy of fuel at final temperature

    #Find FAR corresponding to Tburn
    FAR = (ha - ho) / (hf - ηburn * hc - (1 - ηburn) * hff - abs(hvap)) #solve for FAR 

    #Find product composition and mole and mass fractions
    Xprod_dict = vitiated_mixture(fuel_sps, gas_sps, FAR, ηburn)
    Xprod = Xidict2Array(Xprod_dict)
    Yprod = X2Y(Xprod)

    #Initialize output 
    gas_prod = Gas(Tburn, gas_ox.P)

    #Set the correct composition
    gas_prod.Y = Yprod

    return FAR, gas_prod
end
