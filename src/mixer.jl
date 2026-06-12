"""
    Mixer

Precomputed two-stream mixing system: the pure, allocation-free replacement
for the Dict-based composition step of the legacy [`gas_mixing`](@ref) path.

Built **once** from two composition-bearing streams (construction may consult
the global species database `spdict`; [`mixed`](@ref) calls never do).
Construction resolves each stream to mole fractions, converts them to
per-stream mass fractions, and stores dense per-species `SVector`/`SMatrix`
data over the full species database so that mixture coefficients can be
formed by pure static-array algebra — no Dicts, no string lookups.

```julia-repl
julia> sys = Mixer(DryAir, IdealGasThermo.vitiated_species("CH4", "Air", 0.03));

julia> gas = mixed(sys, 0.25); # FrozenGas of the merged mixture
```

Accepted stream forms (same as the `Combustor` oxidizer): a
[`composite_species`](@ref) (e.g. `DryAir`), a database [`species`](@ref) or
its name (`"Air"` maps to the dry-air composition `Xair`, mirroring the
legacy path), a mole-fraction `Dict{String,Float64}`, or a mole-fraction
vector ordered as `spdict`.
"""
struct Mixer
    name::String
    Y1::SVector{Nspecies,Float64}                # stream-1 mass fractions (Σ = 1)
    Y2::SVector{Nspecies,Float64}                # stream-2 mass fractions (Σ = 1)
    Alow::SMatrix{9,Nspecies,Float64,9 * Nspecies}  # NASA-9 coeffs, T < Tmid
    Ahigh::SMatrix{9,Nspecies,Float64,9 * Nspecies} # NASA-9 coeffs, T ≥ Tmid
    MWvec::SVector{Nspecies,Float64}             # species MW [g/mol]
    Hfvec::SVector{Nspecies,Float64}             # species Hf at 298.15 K [J/mol]
end

"""
    Mixer(stream1, stream2)

Construct the precomputed mixing system for `stream1` and `stream2`
(see [`Mixer`](@ref) for accepted forms). Mole fractions are converted to
per-stream mass fractions here, with the species database molecular weights,
so that [`mixed`](@ref) can apply the mass-fraction law of mixtures
directly. Allocates and consults the species database — do this once,
outside the hot path.
"""
function Mixer(stream1, stream2)
    X1, MW1 = _X_MW(stream1)
    X2, MW2 = _X_MW(stream2)
    name(s) = s isa AbstractSpecies ? s.name : (s isa AbstractString ? s : "stream")
    Mixer(
        "$(name(stream1)) + $(name(stream2))",
        SVector{Nspecies,Float64}(X1 .* spdict.MW ./ MW1),
        SVector{Nspecies,Float64}(X2 .* spdict.MW ./ MW2),
        SMatrix{9,Nspecies,Float64}(reduce(hcat, spdict.alow)),
        SMatrix{9,Nspecies,Float64}(reduce(hcat, spdict.ahigh)),
        SVector{Nspecies,Float64}(spdict.MW),
        SVector{Nspecies,Float64}(spdict.Hf),
    )
end

"""
    mixed(sys::Mixer, mratio) -> FrozenGas

Merged [`FrozenGas`](@ref) of the system `sys` at mass ratio
`mratio = mass₂/mass₁`. Pure function of `(sys, mratio)`: zero allocations,
no global lookups, smooth in `mratio` and generic over `Real` (ForwardDiff
through `mratio` works).

The merged composition follows the mass-fraction law of mixtures
`Y(mratio) = (Y₁ + mratio·Y₂)/(1 + mratio)` (as in the legacy
[`gas_mixing`](@ref)), converted to mole fractions
`X = (Y ./ MW)/Σ(Yᵢ/MWᵢ)`; equivalent NASA-9 mixture coefficients are
formed by mole-fraction weighting with the entropy of mixing
`-Σ Xᵢ ln Xᵢ` of the **merged** composition folded into the integration
constant (b₂), then mass-scaled by `1000/MW` — identical (to rounding) to
`FrozenGas(generate_composite_species(X))`. Same formation-inclusive
enthalpy datum as every `FrozenGas`. `mixed(sys, 0)` is stream 1; stream 2
is the `mratio → ∞` limit.
"""
function mixed(sys::Mixer, mratio)
    Y = (sys.Y1 + mratio * sys.Y2) / (1 + mratio)
    W = Y ./ sys.MWvec   # moles per unit mass of mixture
    X = W / sum(W)

    alow = sys.Alow * X
    ahigh = sys.Ahigh * X
    MW = dot(sys.MWvec, X)
    Hf = dot(sys.Hfvec, X)

    # Entropy of mixing of the merged composition, folded into the
    # integration constant b₂ as in generate_composite_species. Entries with
    # Xᵢ = 0 are identically zero for all mratio (absent from both streams),
    # so the branch does not break mratio-differentiability.
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
