abstract type AbstractGas end
"""
    Gas{N}

A type that represents an ideal gas that is calorically perfect 
i.e. ``c_p(T)``, ``h(T)``, ``\\phi(T)`` and ``s(T,P)``.
"""
mutable struct Gas{N, R<:Real} <: AbstractGas
    P::R                   # Pressure [Pa]
    T::R                   # Temperature [K]
    Tarray::MVector{8, R}  # Preallocated temp array
    cp::R                  # Heat capacity [J/mol/K]
    cp_T::R                # Derivative dcp/dT
    h::R                   # Enthalpy [J/mol]
    ϕ::R                   # Entropy complement function
    Y::MVector{N, R}       # Mass fractions
    MW::R                  # Molecular weight [g/mol]
end

"""
    Gas(Y)

Constructs `Gas` with given composition `Y`

"""
function Gas(Y::AbstractVector{R}) where {R<:Real}
  
    gas = Gas(Pstd, Tstd)
    gas.Y = SVector{Nspecies, R}(Y)
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
    Gas{Nspecies, Float64}(
        Pstd,
        Tstd,
        Tarray(Tstd),
        Cp(Tstd, Air),
        (Cp(Tstd + 1.0, Air) - Cp(Tstd - 1.0, Air)) / 2.0, #finite diff dCp/dT
        h(Tstd, Air),
        s(Tstd, Pstd, Air),
        Y,
        Air.MW
    )

end

function Gas(T::R, P::R) where R<:Real
    i = findfirst(x -> x == "Air", spdict.name)
    Air = spdict[i]

    Y = zeros(R, Nspecies)
    Y[i] = one(R)

    Gas{Nspecies, R}(
        P,
        T,
        Tarray(T),  # must return Vector{R}
        Cp(T, Air),  # must return R
        (Cp(T + one(R), Air) - Cp(T - one(R), Air)) / (2one(R)),  # finite diff in R
        h(T, Air),   # must return R
        s(T, P, Air), # must return R
        SVector{Nspecies, R}(Y),
        Air.MW
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

function Base.setproperty!(gas::Gas{N, R1}, sym::Symbol, val::R2) where {N, R1<:Real, R2<:Real}
    if sym === :T
        setfield!(gas, :T, val)
        setfield!(gas, :Tarray, Tarray!(val, gas.Tarray))
        TT = view(gas.Tarray, :)

        A = val < 1000 ? view(spdict.alow, :) : view(spdict.ahigh, :)

        cptemp = zero(R1)
        htemp = zero(R1)
        ϕtemp = zero(R1)
        cp_Ttemp = zero(R1)

        for (Yᵢ, a, m) in zip(gas.Y, A, spdict.MW)
            if Yᵢ != 0
                cptemp   += Yᵢ * Cp(TT, a) / m
                htemp    += Yᵢ * h(TT, a) / m
                ϕtemp    += Yᵢ * 𝜙(TT, a) / m
                cp_Ttemp += Yᵢ * dCpdT(TT, a) / m
            end
        end

        setfield!(gas, :cp, cptemp * 1000)
        setfield!(gas, :h, htemp * 1000)
        setfield!(gas, :ϕ, ϕtemp * 1000)
        setfield!(gas, :cp_T, cp_Ttemp * 1000)

    elseif sym === :P
        setfield!(gas, :P, val)
        TT = gas.Tarray
        A = TT[4] < 1000 ? view(spdict.alow, :) : view(spdict.ahigh, :)
        ϕtemp = zero(R1)

        for (Yᵢ, a, m) in zip(gas.Y, A, spdict.MW)
            if Yᵢ != 0
                ϕtemp += Yᵢ * 𝜙(TT, a) / m
            end
        end

        setfield!(gas, :ϕ, ϕtemp * 1000)

    elseif sym === :h
        set_h!(gas, val)
    elseif sym === :TP
        set_TP!(gas, val[1], val[2])
    else
        setfield!(gas, sym, val)
    end
    return nothing
end

function Base.setproperty!(gas::Gas{N, R1}, sym::Symbol, val::AbstractVector{<:R2}) where {N, R1<:Real, R2<:Real}
    if sym === :Y
        setfield!(gas, :Y, MVector{N, R1}(val))
        setfield!(gas, :MW, MW(gas))
        gas.T = gas.T #Reset T
    elseif sym === :X
        Y = X2Y(val)
        setfield!(gas, :Y, MVector{N, R1}(Y))
        setfield!(gas, :MW, MW(gas))
        gas.T = gas.T #Reset T
    else
        error("Only mass fractions Y/X can be set with a vector. Tried: $sym")
    end
    return nothing
end

function Base.setproperty!(gas::Gas{N, R1}, sym::Symbol, val::AbstractDict{String, R2}) where {N, R1<:Real, R2<:Real}
    names = spdict.name
    Y = zeros(MVector{N, R1})

    if sym === :Y
        for (key, value) in val
            idx = findfirst(==(key), names)
            Y[idx] = value
        end
        gas.Y = Y
        gas.T = gas.T #Reset T
        gas.MW = MW(gas)

    elseif sym === :X
        X = zeros(MVector{N, R1})
        S = zero(R1)
        for (key, value) in val
            idx = findfirst(==(key), names)
            X[idx] = value
            S += value
        end
        X ./= S
        gas.Y = X2Y(X)
        gas.T = gas.T #Reset T
        gas.MW = MW(gas)
    else
        error("Only mass fractions Y/X can be set with a Dict. Tried: $sym")
    end
    return nothing
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
function set_h!(gas::AbstractGas, hspec::R) where R <: Real
    T = gas.T
    dT = T

    itermax = 20
    eps = convert(R, ϵ)  # tolerance with correct type
    for i = 1:itermax # abs(dT) > ϵ
        res = gas.h - hspec # Residual
        res_t = gas.cp  # ∂R/∂T = ∂h/∂T = cp
        dT = -res / res_t # Newton step

        if abs(dT/T) ≤ eps
            break
        end
        #Prevent limit cycles if the iteration count is high
        if i > itermax/2
            dT = dT * i/itermax #Step can no longer be periodic
        end
        T = T + dT
        gas.T = T
    end

    if abs(dT/T) > eps
        error(
            "Error: `set_h!` did not converge:\ngas=",
            print(gas),
            "\n\nabs(dT) = ",
            abs(dT),
            " > ϵ (",
            eps,
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
function set_Δh!(gas::AbstractGas, Δhspec::R, ηp::R = 1.0) where R <: Real
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
function set_hP!(gas::AbstractGas, hspec::R, P::R) where R <: Real
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
function set_TP!(gas::AbstractGas, T::R, P::R) where R <: Real
    gas.T = T
    gas.P = P
    return gas
end
