#This file contains benchmarks of some functions, primarily to test the effects of changes to the code
using BenchmarkTools
using IdealGasThermo
TT = rand(200.:600., 100)

function benchmark_Gas(TT::AbstractVector, gas::IdealGasThermo.AbstractGas)
    @views for i in eachindex(TT)
        gas.T = TT[i]; 
        gas.cp
        gas.ϕ
        gas.h
        gas.cp_T
    end
end

gas = Gas()
gas1D = Gas1D()
println("Benchmarking Gas1D")
res1 = @benchmark benchmark_Gas($TT, $gas1D)
display(res1)

println("Benchmarking Gas")
res2 = @benchmark benchmark_Gas($TT, $gas)
display(res2)

function benchmark_enthalpy()
    gas = Gas()
    gas.h = 1e5
    return gas
end

println("Benchmarking change in enthalpy")
res3 = @benchmark benchmark_enthalpy()
display(res3)

function benchmark_combustion()
   gas_ox = Gas()
   fuel = "CH4"

   gas = IdealGasThermo.fuel_combustion(
        gas_ox,
        fuel,
        298.15,
        1e-2
        )
    return gas
end

println("Benchmarking combustion")
res4 = @benchmark benchmark_combustion()
display(res4)
