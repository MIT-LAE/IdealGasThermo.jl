using StructArrays
"""
    readThermo(filename)

Reads a NASA 9 polynomial thermo definintion file which can be obtained from
[NASA thermobuild](https://cearun.grc.nasa.gov/ThermoBuild/index_ds.html) 
and returns a dictionary of species.

See [here](https://shepherd.caltech.edu/EDL/PublicResources/sdt/formats/nasa.html) 
for typical data format

Usage:
If a NASA 9 polynomial definition file `thermo.inp` exists then,

```spec = readThermo("thermo.inp")```

will return a dictionary of species.

`readThermo` only considers 2 temperature ranges 
(typically 200-1000 K and 1000-6000 K) but more can be added if needed.
"""
function readThermo(filename)

    f = open(filename)
    input = readlines(f)

    # Get temperatures
    Temps =
        parse.(Float64, filter!(!isempty, split(strip(input[2]), " "))[1:end-1])[1:end-1]
    Tintervals = length(Temps)

    # Count number of species names and lines and 
    # store the species header lines in `sp_lines`
    sp_names, sp_lines = zip(
        [
            [strip(i[1:13]), n] for
            (n, i) in enumerate(input) if startswith(i, r"[a-zA-Z]\w[^END][^thermo]")
        ]...,
    )
    Nspecies = length(sp_names)
    Species = Array{species}(undef, Nspecies)

    a = zeros((2, 9)) #Assume just 2 temp intervals for now
    for (i, sp) in enumerate(sp_names)
        istart = sp_lines[i]
        formula = Iterators.partition(strip(input[istart+1][11:50]), 8)
        comp = split.(join.([i for i in formula]))
        comp = join.(filter!(x -> length(x) > 1, comp)) #Only get part that has species + moles
        formula = join(comp)

        MW = parse(Float64, strip(input[istart+1][53:65]))
        Hf = parse(Float64, strip(input[istart+1][66:80]))
        Tmid = parse(Float64, strip(input[istart+2][12:22]))

        if Tmid != 1000.0
            error("Tmid is not 1000.0")
        end

        coeffs = replace(input[istart+3], "D" => "e") * replace(input[istart+4], "D" => "e")
        a[1, :] =
            parse.(
                Float64,
                [coeffs[(i-1)*16+1:i*16] for i in range(1, 10, step = 1) if i != 8],
            )

        coeffs = replace(input[istart+6], "D" => "e") * replace(input[istart+7], "D" => "e")
        a[2, :] =
            parse.(
                Float64,
                [coeffs[(i-1)*16+1:i*16] for i in range(1, 10, step = 1) if i != 8],
            )

        Species[i] = species(String(sp), Tmid, a[1, :], a[2, :], MW, Hf, formula)
    end
    return StructArray(Species)
end

#Read species data from thermo.inp
const spdict = readThermo(default_thermo_path)
const Nspecies = length(spdict)

species_in_spdict(name::AbstractString) = spdict[findfirst(x -> x == name, spdict.name)]
