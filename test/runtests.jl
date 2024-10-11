using IdealGases
using Test

@testset "IdealGases" verbose = true begin
    include("unit_test_readthermo.jl")
    include("unit_test_mixthermo.jl")
    include("unit_test_composite.jl")
    include("unit_test_vitiated.jl")
    include("unit_test_combustion.jl")
    include("unit_test_humidity.jl")
    include("unit_test_turbo.jl")
    include("unit_test_atmos.jl")
    include("unit_test_utils.jl")
end
