using Documenter, IdealGasThermo

push!(LOAD_PATH, "../src")


makedocs(
    repo = Documenter.Remotes.GitHub("MIT-LAE", "IdealGasThermo.jl"),
    sitename = "IdealGasThermo.jl",
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Getting started" => "tutorials/getting_started.md",
            "Combustion and mixing" => "tutorials/combustion.md",
            "Humid air and saturation" => "tutorials/humid_air.md",
        ],
        "Thermodynamic Data" => "readthermo.md",
        "Gas and Species" => [
            "Species and legacy types" => "gas.md",
            "Immutable pure core (FrozenGas)" => "frozengas.md",
        ],
        "Ideal gas thermodynamics" => [
            "Properties" => "idealgasthermo.md",
            "Derivatives" => "derivatives.md",
        ],
        "Combustion" => "combustion.md",
        "Performance benchmarks" => "benchmark.md",
    ],
    format = Documenter.HTML(;
        mathengine = Documenter.KaTeX(
            Dict(
                :delimiters => [
                    Dict(:left => raw"$", :right => raw"$", :display => false),
                    Dict(:left => raw"$$", :right => raw"$$", :display => true),
                    Dict(:left => raw"\[", :right => raw"\]", :display => true),
                ],
                :macros => Dict(
                    "\\RR" => "\\mathbb{R}",
                    "\\genfuel" =>
                        raw"{\mathrm{C}_{x_{\mathrm{C}}}\mathrm{H}_{x_{\mathrm{H}}}\mathrm{O}_{x_{\mathrm{O}}}\mathrm{N}_{x_{\mathrm{N}}}}",
                    raw"\fst" => raw"f_{\mathrm{stoich.}}",
                    raw"\Xi" => raw"X_{i}",
                    raw"\Ru" => raw"R_{\mathrm{univ.}}",
                    raw"\Pstd" => raw"P_{\mathrm{std}}",
                    raw"\Tstd" => raw"T_{\mathrm{std}}",
                    raw"\cphatR" => raw"\frac{\hat{c}_p^\circ(T)}{\Ru}",
                    raw"\hhatRT" => raw"\frac{\hat{h}^\circ (T)}{\Ru T}",
                    raw"\shatR" => raw"\frac{\hat{s}^\circ(T)}{\Ru}",
                    raw"\phihatR" => raw"\frac{\hat{\phi}^\circ(T)}{\Ru}",
                    raw"\cpbarR" => raw"\overline{\frac{\hat{c}_p^\circ(T)}{\Ru}}",
                    raw"\hbarRT" => raw"\overline{\frac{\hat{h}^\circ (T)}{\Ru T}}",
                    raw"\sbarR" => raw"\overline{\frac{\hat{s}^\circ(T)}{\Ru}}",
                    raw"\phibarR" => raw"\overline{\frac{\hat{\phi}^\circ(T)}{\Ru}}",
                ),
            ),
        ),
    ),
)

deploydocs(repo = "github.com/MIT-LAE/IdealGasThermo.jl.git")
