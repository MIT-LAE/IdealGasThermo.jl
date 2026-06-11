abstract type AbstractFlowStations end

mutable struct FlowStation
    gas::AbstractGas
    A::Float64  # Area
    Ma::Float64 # Mach number
    V::Float64  # Velocity, Ma*a = Ma*√(γRT)
end

"""
"""
function get_static(FS::FlowStation)
    Ps = FS.gas.P / isenPR(FS.gas, FS.Ma)
    Ts = FS.gas.T / isenTR(FS.gas, FS.Ma)
    Ps, Ts
end  # function get_static

"""
"""
function isenTR(gas::AbstractGas, Ma)
    1 + 0.5 * (gas.γ - 1) * Ma^2
end  # function isenTR

"""
Isentropic temperatur ratio 
"""
function isenPR(gas::AbstractGas, Ma)
    k = gas.γ / (gas.γ - 1)
    isenTR(gas, Ma)^(k)
end  # function isenPR
