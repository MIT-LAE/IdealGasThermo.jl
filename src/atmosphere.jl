"""
    standard_atmosphere(z, alt_type = "geometric")
Calculate the gas state at a given desired altitude in the atmosphere using the 1976
US Standard Atmosphere. Model valid between -5 and 86 km in geometric altitude. Accepts 
geometric altitude ("alt_type == "geometric") or geopotential altitude ("alt_type == "geopotential"). 
Model from: `https://ntrs.nasa.gov/citations/19770009539`
"""
function standard_atmosphere(z::R, alt_type = "geometric") where R<:Real

    g = 9.80665 #Acceleration of gravity on Earth surface (z = 0)
    R_e = 6.356766e6 #Radius of Earth

    #Initialize gas (composition does not change)
    T = R(Tstd)
    P = R(Pstd)
    gas = Gas(T, P)
    Rg = gas.R

    #Calculate geopotential altitude
    if alt_type == "geometric"
        H = (R_e * z) / (R_e + z) #Geopotential height

    elseif alt_type == "geopotential"
        H = z
    end

    #Database with 1976 US Standard Atmosphere temperatures and pressures
    H0s = [0, 11e3, 20e3, 32e3, 47e3, 51e3, 71e3] #atm. regions
    λs = [-6.5e-3, 0, 1e-3, 2.8e-3, 0, -2.8e-3, -2e-3] #lapse rates
    T0s = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65] #temperatures at interfaces
    P0s = [101325, 22632.1, 5474.89, 868.019, 110.906, 66.9389, 3.95642] #pressures at interfaces

    #Table lookup to find atmospheric region
    ind = 1
    for i ∈ eachindex(H0s)
        if H > H0s[i]
            ind = i
        end
    end
    #Store base height, temperature, pressure and lapse rate
    H0 = H0s[ind]
    λ = λs[ind]
    T0 = T0s[ind]
    P0 = P0s[ind]

    #Find temperature and pressure using the lapse rate
    #Use closed-form solutions
    if λ != 0
        P = P0 * ((T0 + λ * (H - H0)) / T0)^(-g / (Rg * λ))
        T = T0 + λ * (H - H0)

    elseif λ == 0
        T = T0
        P = P0 * exp(-g * (H - H0) / (Rg * T))
    end

    #set the gas to the correct T and P
    set_TP!(gas, T, P)
    return gas
end
