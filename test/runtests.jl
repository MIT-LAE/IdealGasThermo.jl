using IdealGasThermo
using Test

const IdealGasThermo = IdealGasThermo

@testset "IdealGasThermo" verbose = true begin
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
    include("unit_test_properties.jl")
end
