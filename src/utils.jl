
"""
    Tarray(T)

Function to create the required temperature array
```math
[T^-2, T^-1, 1.0, T, T^2, T^3, T^4, \\log(T)]
```

"""
function Tarray(T)
    return [T^-2, T^-1, 1.0, T, T^2, T^3, T^4, log(T)]
end


"""
    Tarray!(T, TT)

In place Tarray update that returns
[T^-2, T^-1, 1.0, T, T^2, T^3, T^4, log(T)]
"""
function Tarray!(T, TT)
    TT[1] = T^-2    #T^-2
    TT[2] = TT[1] * T #T^-1
    TT[3] = 1.0     #T^0
    TT[4] = T       #T^1
    TT[5] = T * T     #T^2
    TT[6] = T * TT[5] #T^3
    TT[7] = T * TT[6] #T^4
    TT[8] = log(float(T))
    return TT
end

"""
    thermo_table(gas::Gas; 
    Tstart::Float64=Tstd, Tend::Float64=2000.0, Tinterval::Float64=100.0)

Quickly generate a table of cp, h, and s for a gas
"""
function thermo_table(
    gas::Gas;
    Tstart::Float64 = Tstd,
    Tend::Float64 = 2000.0,
    Tinterval::Float64 = 100.0,
)

    Trange = range(Tstart, Tend, step = Tinterval)
    thermo_table(gas, Trange)
end

"""
    thermo_table(gas::Gas, Trange::AbstractVector)

Method to generate cp, h, and s for a given `Trange`.
"""
function thermo_table(gas::Gas, Trange::AbstractVector)
    cp_array = zero(Trange)
    h_array = zero(Trange)
    𝜙_array = zero(Trange)
    s_array = zero(Trange)
    for (i, T) in enumerate(Trange)
        gas.T = T
        cp_array[i] = gas.cp
        h_array[i] = gas.h
        𝜙_array[i] = gas.ϕ
        s_array[i] = gas.s
    end
    return Trange, cp_array, h_array, 𝜙_array, s_array
end

"""
    Y2X(Y::AbstractVector)

Convert from mass fraction Yi to mole fractions Xi
"""
function Y2X(Y::AbstractVector)
    MW = spdict.MW
    num = Y ./ MW
    den = dot(Y, 1 ./ MW)
    return num ./ den
end  # function Y2X

"""
    X2Y(X::AbstractVector)
Convert from mole fraction Xi to mass fractions Yi
"""
function X2Y(X::AbstractVector)
    MW = spdict.MW
    num = X .* MW
    den = dot(X, MW)
    return num ./ den
end  # function X2Y

"""
    Xidict2Array!(Xdict::Dict{AbstractString, AbstractFloat}, X::AbstractVector{AbstractFloat})
Convert a mole fraction dictonary into the given array X with the right order of
compounds.
"""
function Xidict2Array!(Xdict::Dict{String,Float64}, X::AbstractVector)
    names = spdict.name
    for (key, value) in Xdict
        index = findfirst(x -> x == key, names)
        X[index] = value
    end
    X .= X ./ sum(X)

    return X
end  # function Xidict2Array

"""
    Xidict2Array(Xdict::Dict{String, Float64})
Converts the dict into a new array with mole fractions in the right order
"""
function Xidict2Array(Xdict::Dict{String,Float64})
    names = spdict.name
    X = zeros(Float64, Nspecies)
    for (key, value) in Xdict
        index = findfirst(x -> x == key, names)
        X[index] = value
    end
    # X .= X./sum(X)

    return X

end  # function Xidict2ArrayXdict

# Mole fractions (database order, normalized) and molecular weight of a
# composition-bearing input. Accepted forms: a `composite_species`, a database
# `species` or its name (the species "Air" maps to the dry-air composition
# `Xair` but keeps the database MW, exactly mirroring
# `vitiated_mixture`/`fixed_fuel_vitiated_species`), a mole-fraction
# `Dict{String,Float64}`, or a mole-fraction vector ordered as `spdict`.
# Shared by `Combustor` (oxidizer resolution) and `Mixer` (stream resolution).
function _X_MW(input::composite_species)
    X = Xidict2Array(input.composition)
    return X ./ sum(X), input.MW
end
function _X_MW(input::species)
    Xdict = input.name == "Air" ? Xair : Dict(input.name => 1.0)
    X = Xidict2Array(Xdict)
    return X ./ sum(X), input.MW
end
_X_MW(input::AbstractString) = _X_MW(species_in_spdict(input))
function _X_MW(Xdict::Dict)
    X = Xidict2Array(convert(Dict{String,Float64}, Xdict))
    X = X ./ sum(X)
    return X, dot(spdict.MW, X)
end
function _X_MW(X::AbstractVector{<:Real})
    length(X) == Nspecies ||
        error("Mole-fraction vector must have length $Nspecies (spdict order)")
    Xn = Vector{Float64}(X) ./ sum(X)
    return Xn, dot(spdict.MW, Xn)
end
