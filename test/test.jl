include("../src/Gas.jl")
using .IdealGas

Xair = Dict(
    "N2" => 0.78084,
    "Ar" => 0.009365,
    "Air" => 0.0,
    "H2O" => 0.0,
    "CO2" => 0.000319,
    "O2" => 0.209476,
)

Yair = Dict(
    "O2" => 0.231416,
    "Ar" => 0.012916,
    "Air" => 0.0,
    "H2O" => 0.0,
    "CO2" => 0.000484688,
    "N2" => 0.755184,
)

TT = rand(200.0:600.0, 100)
function _test(TT)
    @views for i in eachindex(TT)
        gas.T = TT[i]
    end
end

function benchmark_Gas(TT::AbstractVector, gas)
    @views for i in eachindex(TT)
        # P = rand(101325.:5*101325,1)
        gas.T = TT[i]
        # gas.P = P[1];
        gas.cp
        gas.ϕ
        gas.h
        gas.cp_T
    end
end

"""
"""
function _test_set_h(h_array, gas)
    for h in h_array
        gas.h = h
    end
end  # function _test_set_h
function _test_set_h2(h_array, gas, T_array)
    T_array .= h2T.(h_array)
    for T in T_array
        gas.T = T
    end
end  # function _test_set_h


"""
"""
function _test2(n)
    for T in rand(200.0:600.0, n)
        gas.T = T
        compress(gas, 2.0, 1.0)
    end
end  # function _test2n

function _test3(n)
    for T in rand(200.0:600.0, n)
        gas.T = T
        compress!(gas, 2.0, dP, dT, 1.0)
    end
end  # function _test2n
