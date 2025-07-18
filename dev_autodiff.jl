using IdealGasThermo
using ForwardDiff, StaticArrays

function s(x)
    T = x[1]
    P = x[2]
    
    g = Gas(T, P)
    s = g.s
    return s
end



