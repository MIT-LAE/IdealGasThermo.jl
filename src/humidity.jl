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
) where {type<:AbstractFloat}
    q = specific_humidity(RH, T, P)
    Xwater = q / ε
    Xdict = mergewith(+, Xair, Dict("H2O" => Xwater))

    X = zeros(Float64, Nspecies)
    Xidict2Array!(Xdict, X)

    return generate_composite_species(X, "Wet air with RH = $RH at ($T K; $(P/1000.0) kPa)")
end  # function generate_humid_air

"""
    humid_air(; SH = nothing, RH = nothing, T = Tstd, P = Pstd) -> FrozenGas

Humid-air [`FrozenGas`](@ref): dry air (`Xair`) plus water
vapor. Specify the humidity as exactly one of

  - `SH` — specific humidity ω [kg water / kg dry air], or
  - `RH` — relative humidity [-], converted at temperature `T` [K] and
    pressure `P` `Pa` (standard day by default) through the legacy
    [`saturation_vapor_pressure`](@ref)/[`specific_humidity`](@ref) model
    (August–Roche–Magnus).

The composition logic is the legacy [`generate_humid_air`](@ref) one: water
at `SH/ε` moles per mole of dry air (`ε = MW_H2O/MW_air`) is merged into the
dry-air composition and the result renormalized. `humid_air(SH = 0.0)` is
exactly `FrozenGas(DryAir)`.

This is a constructor, not a hot path: it consults the species database and
allocates. Build the gas once; all its property functions are then pure in
`(gas, T)`.

```julia-repl
julia> wet = humid_air(RH = 0.5, T = 303.15, P = 101325.0);

julia> R(wet) > R(FrozenGas(DryAir)) # water is lighter than air
true
```
"""
function humid_air(; SH = nothing, RH = nothing, T = Tstd, P = Pstd)
    count(isnothing, (SH, RH)) == 1 ||
        error("Specify exactly one of SH (specific humidity) or RH (relative humidity)")
    ω = isnothing(SH) ? specific_humidity(RH, T, P) : SH
    ω ≥ 0 || error("Humidity must be non-negative, got specific humidity $ω")
    Xwater = ω / ε # moles of H2O per mole of dry air
    Xdict = mergewith(+, Xair, Dict("H2O" => Xwater))

    X = zeros(Float64, Nspecies)
    Xidict2Array!(Xdict, X) # normalizes

    name = isnothing(SH) ? "Wet air with RH = $RH at ($T K; $(P/1000.0) kPa)" :
           "Wet air with SH = $ω"
    return FrozenGas(generate_composite_species(X, name))
end  # function humid_air
