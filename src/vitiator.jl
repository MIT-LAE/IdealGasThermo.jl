"""
    Vitiator

Precomputed fuel + oxidizer combustion system: given a fuel–air ratio it
produces the vitiated combustion-product gas via [`products`](@ref). The
pure, allocation-free replacement for the Dict-based [`vitiated_species`](@ref)
path.

Built **once** from a fuel and an oxidizer (construction may consult the
global species database `spdict`; [`products`](@ref) calls never do).
Construction resolves the oxidizer composition and the per-mole-of-fuel
composition change for complete combustion (scaled by the burner efficiency
`ηburn`), and stores **only that reaction stoichiometry** — the NASA-9 mixture
lumping uses the shared module-const basis (`SPALOW`/…), not a per-system copy.

This is the *composition* model of combustion; it is deliberately **not** named
`Combustor`. That noun belongs to the hardware component in a cycle deck
(pressure drop, efficiency map, geometry) so this layer leaves it free.

```julia-repl
julia> sys = Vitiator("CH4", DryAir);

julia> gas = products(sys, 0.03); # FrozenGas of the burnt mixture
```

Accepted oxidizer forms: a [`composite_species`](@ref) (e.g. `DryAir`), a
database [`species`](@ref) or its name (`"Air"` maps to the dry-air
composition `Xair`, mirroring the legacy path), a mole-fraction
`Dict{String,Float64}`, or a mole-fraction vector ordered as `spdict`.
"""
struct Vitiator
    name::String
    ηburn::Float64                               # burner efficiency [-]
    massratio::Float64                           # MW_ox / MW_fuel [-]
    sumΔX::Float64                               # net mole change per mole fuel [-]
    Xin::SVector{Nspecies,Float64}               # oxidizer mole fractions (Σ = 1)
    ΔX::SVector{Nspecies,Float64}                # mole change per mole fuel
end

"""
    Vitiator(fuel, oxidizer=DryAir; ηburn=1.0)

Construct the precomputed combustion system for `fuel`
(an `AbstractString` name or a [`species`](@ref)) burning in `oxidizer`
(see [`Vitiator`](@ref) for accepted forms) with burner efficiency
`ηburn` (fraction of fuel burnt; the remainder appears unreacted in the
products). Allocates and consults the species database — do this once,
outside the hot path. The `AbstractString` form resolves the fuel name and
forwards to the `species` method.
"""
Vitiator(fuel::AbstractString, oxidizer = DryAir; ηburn::Float64 = 1.0) =
    Vitiator(species_in_spdict(fuel), oxidizer; ηburn = ηburn)

function Vitiator(fuel::species, oxidizer = DryAir; ηburn::Float64 = 1.0)
    Xin, MWox = _X_MW(oxidizer)

    # Per-mole-of-fuel composition change, scaled by ηburn; unburnt fuel
    # passes through (mirrors vitiated_mixture).
    nCO2, nN2, nH2O, nO2 = ηburn .* reaction_change_molar_fraction(fuel.name)
    ΔX = zeros(Float64, Nspecies)
    names = spdict.name
    ΔX[findfirst(==(fuel.name), names)] += 1.0 - ηburn
    ΔX[findfirst(==("CO2"), names)] += nCO2
    ΔX[findfirst(==("H2O"), names)] += nH2O
    ΔX[findfirst(==("N2"), names)] += nN2
    ΔX[findfirst(==("O2"), names)] += nO2

    Vitiator(
        "$(fuel.name) + $(oxidizer isa AbstractSpecies ? oxidizer.name : "oxidizer")",
        ηburn,
        MWox / fuel.MW,
        sum(ΔX),
        SVector{Nspecies,Float64}(Xin),
        SVector{Nspecies,Float64}(ΔX),
    )
end

"""
    products(sys::Vitiator, FAR) -> FrozenGas

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
function products(sys::Vitiator, FAR)
    molFAR = FAR * sys.massratio
    X = (sys.Xin + molFAR * sys.ΔX) / (1 + molFAR * sys.sumΔX)
    FrozenGas(X)   # the lumping (basis × X + entropy of mixing) lives in FrozenGas
end

# ── Combustion stoichiometry ──────────────────────────────────────────────────
# These are pure functions of the fuel formula, reachable from the pure-core
# `Vitiator` constructor (above), so they live here rather than in the legacy
# Dict-combustion file (`combustion.jl`), which is scheduled for deletion.

"""
    fuelbreakdown(fuel)

Number of `[C, H, O, N]` atoms in `fuel` — either a chemical-formula string
(e.g. `"CH4"`, `"CH3CH2OH"`) or a database [`species`](@ref) (uses its
`formula`). Used to set up combustion stoichiometry.
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
    reaction_change_molar_fraction(fuel::AbstractString)

Mole-fraction change `[ΔCO2, ΔN2, ΔH2O, ΔO2]` for complete combustion of one
mole of `fuel` (assumed `CᵢHⱼOₖNₗ`):
```
    CᵢHⱼOₖNₗ ⟹ n(CO2)·CO2 + n(H2O)·H2O + n(N2)·N2 − n(O2)·O2
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
    CHON = fuelbreakdown(fuel) # number of C, H, O, and N atoms in fuel
    nCO2 = CHON[1] * 1.0
    nN2 = CHON[4] / 2
    nH2O = CHON[2] / 2
    nO2 = -(CHON[1] + CHON[2] / 4 - CHON[3] / 2) # oxygen is used up
    return [nCO2, nN2, nH2O, nO2]
end
