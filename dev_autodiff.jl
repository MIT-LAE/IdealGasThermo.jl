using IdealGasThermo
using ForwardDiff, StaticArrays

function s(x)
    T = x[1]
    P = x[2]
    
    g = Gas(T, P)
    s = g.s
    return s
end

ForwardDiff.gradient(s, [288.15, 101325.0])

gas_ox = Gas(288.15, 101325.0)


function FAR(x)
    T = x[1]
    P = x[2]
    Tf = x[3]
    Tb = x[4]
    etab = x[5]
    hvap = x[6]
    
    g = Gas(T, P)
    FAR,_ = IdealGasThermo.gas_burn(g,
        "CH4",
        Tf,
        Tb,
        etab,
        hvap)
    return FAR
end

ForwardDiff.gradient(FAR, [288.15, 101325.0, 288.15, 1000.0, 1.0, 0.0])