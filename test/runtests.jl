using IdealGasThermo
using Test
using Aqua

@testset "IdealGasThermo" verbose = true begin
    @testset "Aqua" begin
        Aqua.test_all(IdealGasThermo)
    end
    include("unit_test_readthermo.jl")
    include("unit_test_mixthermo.jl")
    include("unit_test_composite.jl")
    include("unit_test_vitiated.jl")
    include("unit_test_combustion.jl")
    include("unit_test_humidity.jl")
    include("unit_test_turbo.jl")
    include("unit_test_atmos.jl")
    include("unit_test_utils.jl")
    include("unit_test_frozengas.jl")
    include("unit_test_products.jl")
    include("unit_test_mixing.jl")
    include("unit_test_humidair.jl")
    include("unit_test_fastgas.jl")
    include("unit_test_gasstate.jl")
    include("unit_test_flow.jl")
    include("unit_test_properties.jl")
    include("unit_test_cea_reference.jl")
end
