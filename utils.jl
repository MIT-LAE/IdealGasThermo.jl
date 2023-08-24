
"""
Function to create the required temperature array
"""
function Tarray(T)
   return [T^-2, T^-1, 1.0, T, T^2, T^3, T^4, log(T)]
end


"""
    Tarray!(T, TT)

In place Tarray update that returns
[T^-2, T^-1, 1.0, T, T^2, T^3, T^4, log(T)]
"""
function Tarray!(T, TT)
   TT[1] = T^-2    #T^-2
   TT[2] = TT[1]*T #T^-1
   TT[3] = 1.0     #T^0
   TT[4] = T       #T^1
   TT[5] = T*T     #T^2
   TT[6] = T*TT[5] #T^3
   TT[7] = T*TT[6] #T^4
   TT[8] = log(float(T))
   return TT
end

"""
    thermo_table(gas::Gas; 
    Tstart::Float64=Tstd, Tend::Float64=2000.0, Tinterval::Float64=100.0)

Quickly generate a table of cp, h and s for a gas
"""
function thermo_table(gas::Gas; 
    Tstart::Float64=Tstd, Tend::Float64=2000.0, Tinterval::Float64=100.0)
    
    Trange = range(Tstart, Tend, step=Tinterval)
    thermo_table(gas, Trange)
end

function thermo_table(gas::Gas, Trange::AbstractVector)
    cp_array = zero(Trange)
    h_array  = zero(Trange)
    𝜙_array = zero(Trange)
    s_array  = zero(Trange)
    for (i,T) in enumerate(Trange)
        gas.T = T
        cp_array[i]= gas.cp
        h_array[i] =  gas.h
        𝜙_array[i] =  gas.ϕ
        s_array[i] =  gas.s
    end
   return Trange, cp_array, h_array, 𝜙_array, s_array
end


