H2O = species_in_spdict("H2O")
const ε = H2O.MW / DryAir.MW
"""
    saturation_vapor_pressure(T)
Returns the saturation vapor pressure of water in Pa.
See [August-Roche-Magnus formula](https://en.wikipedia.org/wiki/Clausius%E2%80%93Clapeyron_relation#August%E2%80%93Roche%E2%80%93Magnus_approximation).
"""
function saturation_vapor_pressure(T)
    Tc = T - 273.15# convert from K to C
    kern = 17.625 * Tc / (Tc + 243.04)
    return 610.94 * exp(kern)
end  # function saturation_vapor_pressure

"""
    specific_humidity(RH,T,P)
Calculate the specific humidity given the RH, T and P.
"""
function specific_humidity(RH, T, P)
    Psat = saturation_vapor_pressure(T)
    PH2O = Psat * RH
    return ε * PH2O / P
end  # function specific_humidity

function specific_humidity(sp::composite_species)
    comp = sp.composition
    XH2O = comp["H2O"]
    return XH2O * ε / (1 - XH2O)
end
"""
    relative_humidity(Hsp, T, P)
Calculate the relative humidity given specific humidity, T and P.
"""
function relative_humidity(Hsp, T, P)
    PH2O = Hsp * P / ε
    Psat = saturation_vapor_pressure(T)
    return PH2O / Psat
end  # function RH

"""
    generate_humid_air(RH::type, 
    T::type=Tstd, 
    P::type=Pstd) where type<:AbstractFloat

Generates a composite species with the given relative humidity,
temperature, and pressure. Defaults to standard day T, P.
"""
function generate_humid_air(
    RH::type,
    T::type = Tstd,
    P::type = Pstd,
) where {type<:Real}
    q = specific_humidity(RH, T, P)
    Xwater = q / ε
    Xdict = mergewith(+, Xair, Dict("H2O" => Xwater))

    X = zeros(type, Nspecies)
    Xidict2Array!(Xdict, X)

    return generate_composite_species(X, "Wet air with RH = $RH at ($T K; $(P/1000.0) kPa)")
end  # function generate_humid_air
