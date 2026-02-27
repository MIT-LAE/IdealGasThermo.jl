using BenchmarkTools
using IdealGasThermo

const SUITE = BenchmarkGroup()

const TT = rand(200.0:600.0, 100)

function benchmark_gas!(TT::AbstractVector, gas::IdealGasThermo.AbstractGas)
    @views for i in eachindex(TT)
        gas.T = TT[i]
        gas.cp
        gas.ϕ
        gas.h
        gas.cp_T
    end
    return nothing
end

function benchmark_enthalpy()
    gas = Gas()
    gas.h = 1e5
    return gas
end

function benchmark_combustion()
    gas_ox = Gas()
    fuel = "CH4"
    return IdealGasThermo.fuel_combustion(gas_ox, fuel, 298.15, 1e-2)
end

SUITE["gas"]["Gas1D thermo properties"] = @benchmarkable benchmark_gas!($TT, $(Gas1D()))
SUITE["gas"]["Gas thermo properties"] = @benchmarkable benchmark_gas!($TT, $(Gas()))
SUITE["state"]["Set enthalpy"] = @benchmarkable benchmark_enthalpy()
SUITE["combustion"]["fuel_combustion"] = @benchmarkable benchmark_combustion()
