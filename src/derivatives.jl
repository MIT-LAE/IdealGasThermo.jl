using FiniteDifferences

"""
    set_h_with_derivatives!(gas::Gas, h::Float64)

Sets the enthalpy of the gas to `h` and returns dT/dh.
"""
function set_h_with_derivatives!(gas::AbstractGas, h::Float64)
   set_h!(gas, h)
   T_h = 1/gas.cp
   
    return T_h
end

function compress_with_derivatives!(gas::AbstractGas, PR::Float64, ηp::Float64 = 1.0)

    function compress_T(x)
        gas2 = deepcopy(gas)
        PressureRatio!(gas2, x[1], x[2])
        return [gas2.T]
    end
    J = jacobian(forward_fdm(2, 1), compress_T, [PR, ηp])[1]
    
    #Derivatives
    T_PR = J[1,1]
    T_epol = J[1,2]
    P_PR = gas.P
    P_epol = 0.0

    compress!(gas, PR, ηp)
    return T_PR, T_epol, P_PR, P_epol
end

function expand_with_derivatives!(gas::AbstractGas, PR::Float64, ηp::Float64 = 1.0)

    function expand_T(x)
        gas2 = deepcopy(gas)
        PressureRatio!(gas2, x[1], 1/x[2]) #Use inverse efficiency for expansion
        return [gas2.T]
    end
    J = jacobian(forward_fdm(2, 1), expand_T, [PR, ηp])[1]
    
    #Derivatives
    T_PR = J[1,1]
    T_epol = J[1,2]
    P_PR = gas.P
    P_epol = 0.0

    expand!(gas, PR, ηp)
    return T_PR, T_epol, P_PR, P_epol
end

function gas_burn_with_derivatives(
    gas_ox::AbstractGas,
    fuel::String,
    Tf::Float64,
    Tburn::Float64,
    ηburn::Float64 = 1.0,
    hvap::Float64 = 0.0,
)

    function FAR_from_Tf_and_Tburn(x)
        FAR, _ = gas_burn(gas_ox, fuel, x[1], x[2], ηburn, hvap)
    end
    J = jacobian(forward_fdm(2, 1), FAR_from_Tf_and_Tburn, [Tf, Tburn])[1]

    FAR_Tf = J[1, 1]
    FAR_Tburn = J[1, 2]

    FAR, gas_burned = gas_burn(gas_ox, fuel, Tf, Tburn, ηburn, hvap)
    return FAR, gas_burned, FAR_Tf, FAR_Tburn
end

function gas_Mach_with_derivatives!(gas::AbstractGas, M0::Float64, M::Float64, ηp::Float64 = 1.0)
    function find_T_and_P_from_M0_M_and_ηp(x)
        gas2 = deepcopy(gas)
        gas_Mach!(gas2, x[1], x[2], x[3])
        return [gas2.T, gas2.P]
    end

    J = jacobian(forward_fdm(2, 1), find_T_and_P_from_M0_M_and_ηp, [M0, M, ηp])[1]
    #Derivatives
    T_M0 = J[1, 1]
    T_M = J[1, 2]
    T_ηp = J[1, 3]
    P_M0 = J[2, 1]
    P_M = J[2, 2]
    P_ηp = J[2, 3]
    
    gas_Mach!(gas, M0, M, ηp)
    return T_M0, T_M, T_ηp, P_M0, P_M, P_ηp
end
