"""
    Gas1D

Type that represents single component gases.
"""
mutable struct Gas1D <: AbstractGas
    comp_sp::composite_species
    P::Float64 # [Pa]
    T::Float64 # [K]
    Tarray::MVector{8,Float64} # Temperature array to make calcs allocation free

    cp::Float64 #[J/kg/K]
    cp_T::Float64 # derivative dcp/dT
    h::Float64  #[J/kg]
    ϕ::Float64  #[J/kg/K] Entropy complement fn ϕ(T) = ∫ cp(T)/T dT from Tref to T
    #Intentionally not storing s since log is expensive so only calculated when requested

end

"""
    Gas1D()

Constructor that returns a `Gas1D` type representing
Dry Air at standard conditions

!!! warning "Deprecated"
    `Gas1D` is deprecated (see ADR-0002): use [`FrozenGas`](@ref) with the
    pure property functions instead, e.g. `air = FrozenGas(DryAir);
    props(air, T)`. `Gas1D` will be removed in v2.0.

See also [`Gas1D`](@ref).
"""
function Gas1D()
    _legacy_warn(:Gas1D, "FrozenGas(DryAir) with the pure property functions (cp/h/s0/props, temperature)")
    Gas1D(
        DryAir,
        Pstd,
        Tstd,
        Tarray(Tstd),
        Cp(Tstd, DryAir),
        (Cp(Tstd + 1.0, DryAir) - Cp(Tstd - 1.0, DryAir)) / 2.0,
        h(Tstd, DryAir),
        s(Tstd, Pstd, DryAir),
    )
end

"""
    Gas1D(sp::composite_species)

Constructor that returns a `Gas1D` initialized with the given composite
species at standard conditions.

!!! warning "Deprecated"
    `Gas1D` is deprecated (see ADR-0002): use `FrozenGas(sp)` instead.

See also [`Gas1D`](@ref).
"""
function Gas1D(sp::composite_species)
    _legacy_warn(:Gas1D, "FrozenGas(sp) with the pure property functions (cp/h/s0/props, temperature)")
    Gas1D(
        sp,
        Pstd,
        Tstd,
        Tarray(Tstd),
        Cp(Tstd, sp),
        (Cp(Tstd + 1.0, sp) - Cp(Tstd - 1.0, sp)) / 2.0,
        h(Tstd, sp),
        s(Tstd, Pstd, sp),
    )
end

function Base.setproperty!(gas::Gas1D, sym::Symbol, val::Float64)
    ## Setting Temperature
    if sym === :T
        setfield!(gas, :T, val) # first set T
        setfield!(gas, :Tarray, Tarray!(val, getfield(gas, :Tarray))) # update Tarray
        TT = view(getfield(gas, :Tarray), :) # Just convinence
        # Next set the cp, h and s of the gas
        ## Get the right coefficients 
        ## (assumes Tmid is always 1000.0. Check performed in readThermo.jl.):
        sp = getfield(gas, :comp_sp)
        if val < 1000.0 #Value given is the desired temperature
            A = view(sp.alow, :)
        else
            A = view(sp.ahigh, :)
        end

        P = getfield(gas, :P)
        MW = sp.MW / 1000.0 # kg/mol

        setfield!(gas, :cp, Cp(TT, A) / MW)
        setfield!(gas, :h, h(TT, A) / MW)
        setfield!(gas, :ϕ, 𝜙(TT, A) / MW)
        setfield!(gas, :cp_T, dCpdT(TT, A) / MW)

        ## Setting Pressure
    elseif sym === :P
        setfield!(gas, :P, val)
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

# Overload Base.getproperty for convinence
function Base.getproperty(gas::Gas1D, sym::Symbol)
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
        Rgas = Runiv / getfield(gas, :comp_sp).MW * 1000.0
        return getfield(gas, :ϕ) - Rgas * log(getfield(gas, :P) / Pstd)
    elseif sym === :MW
        sp = getfield(gas, :comp_sp)
        return sp.MW
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
    else
        return getfield(gas, sym)
    end
end
