using IdealGasThermo
using ForwardDiff, StaticArrays

function s(x)
    T = x[1]
    P = x[2]
    
    g = Gas(T, P)
    s = g.s
    return s
end

ForwardDiff.gradient(s, [288.15, 101325.0])

gas_ox = Gas(288.15, 101325.0)


function FAR(x)
    T = x[1]
    P = x[2]
    Tf = x[3]
    Tb = x[4]
    
    g = Gas(T, P)
    FAR,_ = IdealGasThermo.gas_burn(g,
        "CH4",
        Tf,
        Tb)
    return FAR
end

ForwardDiff.gradient(FAR, [288.15, 101325.0, 288.15, 1000.0])

function PR(x)
    T = x[1]
    P = x[2]
    PR = x[3]
    eta_p = x[4]
    
    g = Gas(T, P)
    IdealGasThermo.PressureRatio!(g, PR, eta_p)
    return g.T
end

ForwardDiff.gradient(PR, [288.15, 101325.0, 2.0, 1.0])

function FAR_dev(x)
    T = x[1]
    P = x[2]
    Tf = x[3]
    Tburn = x[4]
    ηburn= 1.0
    hvap = 0.0

    gas_ox = Gas(T, P)
    fuel = "CH4"

    #Create variables corresponding to the oxidizer and fuel species and mixtures
    fuel_sps = IdealGasThermo.species_in_spdict(fuel)

    #Extract composite species with oxidizer gas composition
    if gas_ox isa Gas1D
        gas_sps = gas_ox.comp_sp
    else
        if "Air" in keys(gas_ox.Xdict)
            Xin = IdealGasThermo.Xair
        else
            Xin = gas_ox.Xdict
        end
        gas_sps = IdealGasThermo.generate_composite_species(IdealGasThermo.Xidict2Array(Xin))
    end

    #Find the vectors with the fuel mole and mass fractions
    Xfuel = IdealGasThermo.Xidict2Array(Dict([(fuel, 1.0)])) #Mole fraction
    Yfuel = IdealGasThermo.X2Y(Xfuel) #Mass fraction
    gas_fuel = Gas(Tf, gas_ox.P)
    gas_fuel.Y = Yfuel #Create a fuel gas to calculate enthalpy

    #Store enthalpies of oxidizer and fuel at original temperatures
    ho = gas_ox.h
    hf = gas_fuel.h
    #println(ho)

    #Find change in gas composition for a FAR of 1
    nCO2, nN2, nH2O, nO2 = IdealGasThermo.reaction_change_molar_fraction(fuel_sps.name)

    names = ["CO2", "H2O", "N2", "O2"]
    ΔX = [nCO2, nH2O, nN2, nO2]
    Xdict = Dict(zip(names, ΔX))

    Xc = IdealGasThermo.Xidict2Array(Xdict)
    Yc = IdealGasThermo.X2Y(Xc) #mass fraction change in combustion for FAR = 1

    gas_c = Gas(Tburn, gas_ox.P) #Create a "virtual" gas with the changes in combustion, for enthalpy
    #calculations
    gas_c.Y = Yc

    hc = gas_c.h #Enthalpy change for FAR = 1
    #println(hc)

    gas_ox_burnt = deepcopy(gas_ox) #make a copy to avoid modifying the input
    set_TP!(gas_ox_burnt, Tburn, gas_ox.P)
    ha = gas_ox_burnt.h #Enthalpy of original oxidizer at final temperature

    set_TP!(gas_fuel, Tburn, gas_ox.P)
    hff = gas_fuel.h #Enthalpy of fuel at final temperature
    #println(hf)

    #Find FAR corresponding to Tburn
    FAR = (ha - ho) / (hf - ηburn * hc - (1 - ηburn) * hff - abs(hvap)) #solve for FAR 

    #Find product composition and mole and mass fractions
    Xprod_dict = IdealGasThermo.vitiated_mixture(fuel_sps, gas_sps, FAR, ηburn)
    Xprod = IdealGasThermo.Xidict2Array(Xprod_dict)
    Yprod = IdealGasThermo.X2Y(Xprod)

    #Initialize output 
    gas_prod = Gas(Tburn, gas_ox.P)

    #Set the correct composition
    gas_prod.Y = Yprod
    return FAR
end