"""
    Combustor

Precomputed fuel + oxidizer combustion system: the pure, allocation-free
replacement for the Dict-based [`vitiated_species`](@ref) path.

Built **once** from a fuel and an oxidizer (construction may consult the
global species database `spdict`; [`products`](@ref) calls never do).
Construction resolves the fuel, computes the per-mole-of-fuel composition
change for complete combustion (scaled by the burner efficiency `ηburn`),
and stores dense per-species `SVector`/`SMatrix` data over the full species
database so that mixture coefficients can be formed by pure static-array
algebra — no Dicts, no string lookups.

```julia-repl
julia> sys = Combustor("CH4", DryAir);

julia> gas = products(sys, 0.03); # FrozenGas of the burnt mixture
```

Accepted oxidizer forms: a [`composite_species`](@ref) (e.g. `DryAir`), a
database [`species`](@ref) or its name (`"Air"` maps to the dry-air
composition `Xair`, mirroring the legacy path), a mole-fraction
`Dict{String,Float64}`, or a mole-fraction vector ordered as `spdict`.
"""
struct Combustor
    name::String
    ηburn::Float64                               # burner efficiency [-]
    massratio::Float64                           # MW_ox / MW_fuel [-]
    sumΔX::Float64                               # net mole change per mole fuel [-]
    Xin::SVector{Nspecies,Float64}               # oxidizer mole fractions (Σ = 1)
    ΔX::SVector{Nspecies,Float64}                # mole change per mole fuel
    Alow::SMatrix{9,Nspecies,Float64,9 * Nspecies}  # NASA-9 coeffs, T < Tmid
    Ahigh::SMatrix{9,Nspecies,Float64,9 * Nspecies} # NASA-9 coeffs, T ≥ Tmid
    MWvec::SVector{Nspecies,Float64}             # species MW [g/mol]
    Hfvec::SVector{Nspecies,Float64}             # species Hf at 298.15 K [J/mol]
end

"""
    Combustor(fuel, oxidizer=DryAir; ηburn=1.0)

Construct the precomputed combustion system for `fuel`
(an `AbstractString` name or a [`species`](@ref)) burning in `oxidizer`
(see [`Combustor`](@ref) for accepted forms) with burner efficiency
`ηburn` (fraction of fuel burnt; the remainder appears unreacted in the
products). Allocates and consults the species database — do this once,
outside the hot path.
"""
function Combustor(
    fuel::Union{AbstractString,species},
    oxidizer = DryAir;
    ηburn::Float64 = 1.0,
)
    fuelsp = fuel isa species ? fuel : species_in_spdict(fuel)
    Xin, MWox = _X_MW(oxidizer)

    # Per-mole-of-fuel composition change, scaled by ηburn; unburnt fuel
    # passes through (mirrors vitiated_mixture).
    nCO2, nN2, nH2O, nO2 = ηburn .* reaction_change_molar_fraction(fuelsp.name)
    ΔX = zeros(Float64, Nspecies)
    names = spdict.name
    ΔX[findfirst(==(fuelsp.name), names)] += 1.0 - ηburn
    ΔX[findfirst(==("CO2"), names)] += nCO2
    ΔX[findfirst(==("H2O"), names)] += nH2O
    ΔX[findfirst(==("N2"), names)] += nN2
    ΔX[findfirst(==("O2"), names)] += nO2

    Combustor(
        "$(fuelsp.name) + $(oxidizer isa AbstractSpecies ? oxidizer.name : "oxidizer")",
        ηburn,
        MWox / fuelsp.MW,
        sum(ΔX),
        SVector{Nspecies,Float64}(Xin),
        SVector{Nspecies,Float64}(ΔX),
        SMatrix{9,Nspecies,Float64}(reduce(hcat, spdict.alow)),
        SMatrix{9,Nspecies,Float64}(reduce(hcat, spdict.ahigh)),
        SVector{Nspecies,Float64}(spdict.MW),
        SVector{Nspecies,Float64}(spdict.Hf),
    )
end

"""
    products(sys::Combustor, FAR) -> FrozenGas

Combustion-product [`FrozenGas`](@ref) of the system `sys` at fuel–air
(mass) ratio `FAR`. Pure function of `(sys, FAR)`: zero allocations,
no global lookups, smooth in `FAR` and generic over `Real` (ForwardDiff
through `FAR` works).

The product composition is
`X(FAR) = (Xin + molFAR·ΔX) / (1 + molFAR·ΣΔX)` with
`molFAR = FAR·MW_ox/MW_fuel`; equivalent NASA-9 mixture coefficients are
formed by mole-fraction weighting with the entropy of mixing
`-Σ Xᵢ ln Xᵢ` folded into the integration constant (b₂), then mass-scaled
by `1000/MW` — identical (to rounding) to
`FrozenGas(vitiated_species(fuel, oxidizer, FAR))`. Same formation-inclusive
enthalpy datum as every `FrozenGas`.
"""
function products(sys::Combustor, FAR)
    molFAR = FAR * sys.massratio
    X = (sys.Xin + molFAR * sys.ΔX) / (1 + molFAR * sys.sumΔX)

    alow = sys.Alow * X
    ahigh = sys.Ahigh * X
    MW = dot(sys.MWvec, X)
    Hf = dot(sys.Hfvec, X)

    # Entropy of mixing, folded into the integration constant b₂ as
    # in generate_composite_species. Entries with Xᵢ = 0 are identically
    # zero for all FAR (never produced/consumed), so the branch does not
    # break FAR-differentiability.
    Δs_mix = zero(eltype(X))
    @inbounds for i in eachindex(X)
        Xi = X[i]
        if !iszero(Xi)
            Δs_mix += Xi * log(Xi)
        end
    end
    alow = Base.setindex(alow, alow[9] - Δs_mix, 9)
    ahigh = Base.setindex(ahigh, ahigh[9] - Δs_mix, 9)

    scale = 1000 / MW # molar (J/mol) → mass-specific (J/kg)
    FrozenGas(alow * scale, ahigh * scale, MW, 1000 * Runiv / MW, Hf)
end
