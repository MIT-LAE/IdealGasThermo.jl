abstract type AbstractGas end
"""
    Gas{N}

A type that represents an ideal gas that is calorically perfect 
i.e. ``c_p(T)``, ``h(T)``, ``\\phi(T)`` and ``s(T,P)``.
"""
mutable struct Gas{N} <: AbstractGas
    P::Float64 # [Pa]
    T::Float64 # [K]
    Tarray::MVector{8,Float64} # Temperature array to make calcs allocation free

    cp::Float64 #[J/mol/K]
    cp_T::Float64 # derivative dcp/dT
    h::Float64  #[J/mol]
    ϕ::Float64  #[J/mol/K] Entropy complement fn ϕ(T) = ∫ cp(T)/T dT from Tref to T

    Y::MVector{N,Float64} # Mass fraction of species
    MW::Float64 # Molecular weight [g/mol]
end

"""
    Gas(Y)

Constructs `Gas` with given composition `Y`

"""
function Gas(Y::AbstractVector)
    gas = Gas()
    gas.Y = convert(Vector{Float64}, Y)
    set_TP!(gas, Tstd, Pstd) #setting temperature and pressure to recalculate thermodynamic properties
    return gas
end

"""
    Gas()

Constructor that returns a `Gas` type representing 
Air at standard conditions

See also [`Gas`](@ref).

# Examples
```julia-repl
julia> Gas()
Ideal Gas at
  T =  298.150 K
  P =  101.325 kPa
 cp =   29.102 J/K/mol
  h =   -0.126 kJ/mol
  s =    0.199 kJ/K/mol

with composition:
-----------------------------
 Species        Yᵢ  MW[g/mol]
-----------------------------
     Air     1.000     28.965
-----------------------------
     ΣYᵢ     1.000     28.965
```
"""
function Gas()
    i = findfirst(x -> x == "Air", spdict.name)
    Air = spdict[i]
    Y = zeros(Nspecies)
    Y[i] = 1.0
    Gas{Nspecies}(
        Pstd,
        Tstd,
        Tarray(Tstd),
        Cp(Tstd, Air),
        (Cp(Tstd + 1.0, Air) - Cp(Tstd - 1.0, Air)) / 2.0, #finite diff dCp/dT
        h(Tstd, Air),
        s(Tstd, Pstd, Air),
        Y,
        Air.MW,
    )

end

# Overload Base.getproperty for convenience
function Base.getproperty(gas::Gas, sym::Symbol)
    if sym === :h_T # dh/dT
        return getfield(gas, :cp)
    elseif sym === :ϕ_T # dϕ/dT
        return getfield(gas, :cp) / getfield(gas, :T)
    elseif sym === :s_T # ∂s/∂T = dϕ/dT
        return getfield(gas, :cp) / getfield(gas, :T)
    elseif sym === :hs
        return [getfield(gas, :h), getfield(gas, :s)]
    elseif sym === :TP
        return [getfield(gas, :T), getfield(gas, :P)]
    elseif sym === :s
        Xi = view(getproperty(gas, :X), :)
        Δs_mix = 0.0
        Rgas = Runiv / getfield(gas, :MW) * 1000.0
        for i in eachindex(Xi)
            if Xi[i] != 0.0
                Δs_mix = Δs_mix + Xi[i] * log(Xi[i])
            end
        end
        return getfield(gas, :ϕ) - Rgas * (log(getfield(gas, :P) / Pstd) + Δs_mix)
    elseif sym === :Hf # formation enthalpy [J/mol]
        Xi = view(getproperty(gas, :X), :)
        Hf = view(spdict.Hf, :)
        H = 0.0
        for i in eachindex(Xi)
            H += Xi[i] * Hf[i]
        end
        # H += getproperty(gas, :h) * getproperty(gas, :MW)/1000.0
        return H
    elseif sym === :R #specific gas constant
        return Runiv / getproperty(gas, :MW) * 1000.0
    elseif sym === :γ
        R = getproperty(gas, :R)
        cp = getproperty(gas, :cp)
        return cp / (cp - R)
    elseif sym === :gamma
        return getproperty(gas, :γ)
    elseif sym === :ρ
        R = getproperty(gas, :R)
        T = getfield(gas, :T)
        P = getfield(gas, :P)
        return P / (R * T)
    elseif sym === :rho
        return getproperty(gas, :ρ)
    elseif sym === :ν
        return 1 / getproperty(gas, :ρ)
    elseif sym === :nu
        return getproperty(gas, :ν)
    elseif sym === :X # Get mole fractions
        Y = getfield(gas, :Y)
        MW = spdict.MW
        num = Y ./ MW
        den = dot(Y, 1 ./ MW)
        return num ./ den
    elseif sym === :Xdict
        X = getproperty(gas, :X)
        index = X .!= 0.0
        names = view(spdict.name, :)
        return Dict(zip(names[index], X[index]))
    elseif sym === :Ydict
        Y = getproperty(gas, :Y)
        names = view(spdict.name, :)
        return Dict(zip(names, Y))
    else
        return getfield(gas, sym)
    end
end


"""
    Base.setproperty!(gas::Gas, s::Symbol, val)

"""
function Base.setproperty!(gas::Gas, sym::Symbol, val::Float64)
    ## Setting Temperature
    if sym === :T
        setfield!(gas, :T, val) # first set T
        setfield!(gas, :Tarray, Tarray!(val, getfield(gas, :Tarray))) # update Tarray
        TT = view(getfield(gas, :Tarray), :) # Just convenience
        # Next set the cp, h and s of the gas
        ## Get the right coefficients 
        ## (assumes Tmid is always 1000.0. Check performed in readThermo.jl.):
        if val < 1000.0 #Value given is the desired temperature
            A = view(spdict.alow, :)
        else
            A = view(spdict.ahigh, :)
        end
        ## Initialize temporary vars
        cptemp = 0.0
        htemp = 0.0
        ϕtemp = 0.0
        cp_Ttemp = 0.0

        P = getfield(gas, :P)
        Y = getfield(gas, :Y)
        MW = view(spdict.MW, :) # g/mol
        # Go through every species where mass fraction is not zero
        @inbounds for (Yᵢ, a, m) in zip(Y, A, MW)
            if Yᵢ != 0.0
                cptemp = cptemp + Yᵢ * Cp(TT, a) / m
                htemp = htemp + Yᵢ * h(TT, a) / m
                ϕtemp = ϕtemp + Yᵢ * 𝜙(TT, a) / m
                cp_Ttemp = cp_Ttemp + Yᵢ * dCpdT(TT, a) / m
            end
        end

        setfield!(gas, :cp, cptemp * 1000.0)
        setfield!(gas, :h, htemp * 1000.0)
        setfield!(gas, :ϕ, ϕtemp * 1000.0)
        setfield!(gas, :cp_T, cp_Ttemp * 1000.0)

        ## Setting Pressure
    elseif sym === :P
        setfield!(gas, :P, val)
        TT = view(getfield(gas, :Tarray), :) # Just convenience
        # Next set s of the gas
        ## Get the right coefficients 
        ## (assumes Tmid is always 1000.0. Check performed in readThermo.jl.):
        if TT[4] < 1000.0 #TT[4] is T
            A = view(spdict.alow, :)
        else
            A = view(spdict.ahigh, :)
        end
        ## Initialize temporary vars
        ϕtemp = 0.0

        P = val
        Y = view(getfield(gas, :Y), :)
        MW = view(spdict.MW, :) # g/mol
        # Go through every species where mass fraction is not zero
        @inbounds for (Yᵢ, a, m) in zip(Y, A, MW)
            if Yᵢ != 0.0
                ϕtemp = ϕtemp + Yᵢ * 𝜙(TT, a) / m
            end
        end

        setfield!(gas, :ϕ, ϕtemp * 1000.0)

    elseif sym === :h
        set_h!(gas, val)
    elseif sym === :TP
        set_TP!(gas, val[1], val[2])
    else
        error("You tried setting gas.", sym, " to a ", typeof(val))
    end
    # Note: intentionally not including other variables to prevent 
    # users from trying to directly set h, s, cp, MW etc.
    nothing
end

function Base.setproperty!(gas::Gas, sym::Symbol, val::AbstractVector{Float64})
    if sym === :Y # directly set mass fractions Y
        n = length(getfield(gas, :Y))
        setfield!(gas, :Y, MVector{n}(val))
        setfield!(gas, :MW, MW(gas)) # Update the MW of the gas mixture
        setfield!(gas, :T, gas.T) # Force update of h

    elseif sym === :X # directly set mole fractions Y
        n = length(getfield(gas, :Y))
        Y = X2Y(val)
        setfield!(gas, :Y, MVector{n}(Y))
        setfield!(gas, :MW, MW(gas)) # Update the MW of the gas mixture
        setfield!(gas, :T, gas.T) # Force update of h
    else
        error(
            "Only mass factions Y can be set with an array.",
            "\n       You tried to set gas.",
            sym,
            " to a ",
            typeof(val),
        )
    end
    nothing
end

function Base.setproperty!(gas::Gas, sym::Symbol, val::AbstractDict{String,Float64})
    if sym === :Y # directly set mass fractions Y
        # If dict provided set each species in the right order
        names = spdict.name
        Y = zeros(MVector{length(names),Float64})
        for (key, value) in val
            index = findfirst(x -> x == key, names)
            Y[index] = value
        end
        setfield!(gas, :Y, Y)
        setfield!(gas, :MW, MW(gas)) # Update the MW of the gas mixture
    elseif sym === :X
        names = spdict.name
        X = zeros(MVector{length(names),Float64})
        Y = zeros(MVector{length(names),Float64})
        S = 0.0
        for (key, value) in val
            index = findfirst(x -> x == key, names)
            S += value
            X[index] = value
        end
        X = X ./ S
        Y .= X2Y(X)
        setfield!(gas, :Y, Y)
        setfield!(gas, :MW, MW(gas))
    else
        error(
            "Only mass factions Y can be set with a dict.",
            "\n       You tried to set gas.",
            sym,
            " to a ",
            typeof(val),
        )
    end
    nothing
end

"""
    MW(g::Gas)

Calculates mean molecular weight of the gas
"""
@views function MW(g::Gas)
    MW = 1 / dot(g.Y, 1 ./ spdict.MW)
    return MW
end


"""
    set_h!(gas::AbstractGas, hspec::Float64)

Calculates gas temperature for a specified enthalpy via a non-linear 
Newton-Raphson method.

# Examples
```julia-repl
julia> gas = Gas();
julia> set_h!(gas, 0.0)
Ideal Gas at
  T =  302.463 K
  P =  101.325 kPa
 cp =   29.108 J/K/mol
  h =    0.000 kJ/mol
  s =    0.199 kJ/K/mol

with composition:
-----------------------------
 Species        Yᵢ  MW[g/mol]
-----------------------------
     Air     1.000     28.965
-----------------------------
     ΣYᵢ     1.000     28.965
```
"""
function set_h!(gas::AbstractGas, hspec::Float64)
    T = gas.T
    dT = T

    itermax = 20
    for i = 1:itermax # abs(dT) > ϵ
        res = gas.h - hspec # Residual
        res_t = gas.cp  # ∂R/∂T = ∂h/∂T = cp
        dT = -res / res_t # Newton step

        if abs(dT) ≤ ϵ
            break
        end
        if i > itermax/2
            dT = dT * i/itermax
        end
        T = T + dT
        gas.T = T
    end

    if abs(dT) > ϵ
        error(
            "Error: `set_h!` did not converge:\ngas=",
            print(gas),
            "\n\nabs(dT) = ",
            abs(dT),
            " > ϵ (",
            ϵ,
            ")",
        )
    end

    return gas
end
"""
    set_Δh!(gas::AbstractGas, Δhspec::Float64, ηp::Float64 = 1.0)

Sets the gas state based on a specified change in enthalpy (Δh) [J/mol],
and a given polytropic efficiency. This represents adding or removing some work
from the gas.
"""
function set_Δh!(gas::AbstractGas, Δhspec::Float64, ηp::Float64 = 1.0)
    P0 = gas.P
    ϕ0 = gas.ϕ
    hf = gas.h + Δhspec
    set_h!(gas, hf)
    gas.P = P0 * exp(ηp / gas.R * (gas.ϕ - ϕ0))
    return gas
end
"""
    set_hP!(gas::AbstractGas, hspec::Float64, P::Float64)

Calculates state of the gas given enthalpy and pressure (h,P)
"""
function set_hP!(gas::AbstractGas, hspec::Float64, P::Float64)
    set_h!(gas, hspec)
    gas.P = P
    return gas
end
"""
    set_TP!(gas::AbstractGas, T::Float64, P::Float64)

Calculates state of the gas given Temperature and pressure (T,P)
in K and Pa respectively.

# Examples
```julia-repl
julia> gas = Gas(); # Create an ideal gas consisting of air at std. conditions
julia> set_TP!(gas, 298.15*2, 101325.0*2)
Ideal Gas at
  T =  596.300 K
  P =  202.650 kPa
 cp =   30.418 J/K/mol
  h =    8.706 kJ/mol
  s =    0.214 kJ/K/mol

with composition:
-----------------------------
 Species        Yᵢ  MW[g/mol]
-----------------------------
     Air     1.000     28.965
-----------------------------
     ΣYᵢ     1.000     28.965
```

"""
function set_TP!(gas::AbstractGas, T::Float64, P::Float64)
    gas.T = T
    gas.P = P
    return gas
end
