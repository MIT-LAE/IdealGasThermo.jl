# Turbomachinery related functions

"""
    PressureRatio(gas::AbstractGas, PR::Float64, ηp::Float64=1.0,)

Generic pressure ratio conversion. See [`compress`](@ref) and [`expand`](@ref)
"""
function PressureRatio(gas::AbstractGas, PR::Float64, ηp::Float64=1.0,)

    T0 = gas.T
    ϕ0 = gas.ϕ
    P0 = gas.P
    R = gas.R
 
    Tfinal = T0 * PR^(R/gas.cp/ηp)
    Pfinal = P0*PR
    dT = Tfinal
    gas.P = Pfinal
    gas.T = Tfinal
    logPR_ηp = log(PR)/ηp
 
    for i in 1:20# abs(dT)>ϵ
       res  = (gas.ϕ - ϕ0)/R - logPR_ηp 
       res_dT = gas.ϕ_T/R
       dT  = - res/res_dT
 
       if abs(dT) ≤ ϵ
          break
       end
 
       Tfinal = Tfinal + dT
       gas.T = Tfinal
       # println("$i: $Tfinal $dT")
    end
 
    if abs(dT) > ϵ
       error("Error: PressureRatio did not converge:\ngas=", print(gas),
        "\n\nabs(dT) = ", abs(dT), " > ϵ (", ϵ, ")")
    end
 
    return gas
 
 end
 
"""
   compress(gas::AbstractGas, PR::Float64, ηp::Float64=1.0,)

Compression with an optional polytropic efficiency.PR should be ≥ 1.0.
See also [`expand`](@ref).
"""
function compress(gas::AbstractGas, PR::Float64, ηp::Float64=1.0,)
   if PR < 1.0
      error("The specified pressure ratio (PR) to compress by needs to be ≥ 1.0.
   Provided PR = $PR. Did you mean to use `expand`?")
   end
   return PressureRatio(gas, PR, ηp)
end

"""
   expand(gas::AbstractGas, PR::Float64, ηp::Float64=1.0,)

Expansion at a given polytropic efficiency. PR should be ≤ 1.0.
See also [`compress`](@ref).
"""
function expand(gas::AbstractGas, PR::Float64, ηp::Float64=1.0,)
   if PR>1.0
      error("The specified pressure ratio (PR) to compress by needs to be ≤ 1.0.
   Provided PR = $PR. Did you mean to use `compress`?")
   end
   return PressureRatio(gas, PR, 1/ηp)
end

"""
   gas_Mach!(gas::AbstractGas, M0::Float64, M::Float64, ηp::Float64 = 1.0)

Calculates the gas state for a change in Mach number with an optional polytropic efficiency.
"""
function gas_Mach!(gas::AbstractGas, M0::Float64, M::Float64, ηp::Float64 = 1.0)
   itmax = 10
   ttol = 0.000001

   uosq = M0^2 * gas.γ * gas.R * gas.T
   h0 = gas.h
   R0 = gas.R 
   ϕ0 = gas.ϕ

   #---- initial guess for temperature, using constant-gamma relation
   t = gas.T * (1.0 + 0.5 * (gas.γ - 1) * M0^2) /
         (1.0 + 0.5 * (gas.γ - 1) * M^2)

   dt = 0.0
   #---- Newton iteration for actual temperature
   for iter = 1:itmax
      set_TP!(gas, t, gas.P) #keep pressure constant for now; will be changed later
      
      #------ usq( t, cp(t))
      usq = M^2 * gas.γ * gas.R * gas.T
      usq_t = M^2 * gas.γ * gas.R
      usq_cp = M^2 * gas.R / (gas.cp - gas.R) * gas.T - usq / (gas.cp - gas.R)

      #------ res( h(t), usq( t, cp(t) ) )
      res = gas.h + 0.5 * usq - h0 - 0.5 * uosq
      res_t = gas.cp + 0.5 * (usq_t + usq_cp * gas.cp_T)
      #
      dt = -res / res_t

      if (abs(dt) < ttol)
         P = gas.P * exp(ηp * (gas.ϕ - ϕ0) / R0) #calculate pressure with polytropic efficiency
         set_TP!(gas, gas.T, P)
         return
      end

      t = t + dt
   end
end

"""
   gas_Deltah(gas_in::AbstractGas, deltah::Float64, epol::Float64 = 1.0)

Calculates the gas state for a given change in specific enthalpy (in compression or expansion) and at 
a given polytropic efficiency. In the case of constant specific heats, this reduces to the classical 
isentropic relations.
"""
function gas_Deltah(gas_in::AbstractGas, deltah::Float64, epol::Float64 = 1.0)
   itmax = 10
   ttol = 0.000001
   gas = deepcopy(gas_in)

   T0 = gas.T
   cp0 = gas.cp
   p0 = gas.P
   h0 = gas.h
   s0 = gas.s

   T = T0 + deltah / cp0

   dT = 0.0
   for iter = 1:itmax
         s = gas.s
         h = gas.h
         h_T = gas.h_T
         cp = gas.cp
         R = gas.R
         res = h - h0 - deltah
         res_T = h_T

         dT = -res / res_T

         if (abs(dT) < ttol)
            p = p0 * exp(epol * (s - s0) / R)
            set_TP!(gas, T, p)
            return gas
         end

         T = T + dT
         set_TP!(gas, T, p0)
   end
      println("gas_Deltah: convergence failed.  dT =", dT)

end # gas_Deltah

"""
   gas_mixing(gas1::AbstractGas, gas2::AbstractGas, mratio::Float64)

Calculates the resulting gas after two gases (gas1 and gas2) are mixed at constant pressure, with a mass ratio
mratio = mass of gas2 / mass gas1.
"""
function gas_mixing(gas1::AbstractGas, gas2::AbstractGas, mratio::Float64)

   #Extract dictionaries with gas molar fractions
   if typeof(gas1) == Gas1D
      X1 = gas1.comp_sp.composition
   else
      if "Air" in keys(gas1.Xdict)
         X1 = Xair
      else
         X1 = gas1.Xdict
      end
  end

   if typeof(gas2) == Gas1D
      X2 = gas2.comp_sp.composition
   else
      if "Air" in keys(gas2.Xdict)
         X2 = Xair
      else
         X2 = gas2.Xdict
      end
   end

   X1arr = Xidict2Array(X1)
   X2arr = Xidict2Array(X2)

   Y1 = X2Y(X1arr) #Vectors with mass fractions
   Y2 = X2Y(X2arr) #Vectors with mass fractions

   Yp = (Y1 + mratio * Y2) / (1 + mratio) #law of mixtures for mass fractions

   #Initialize output 
   gas_prod = Gas(Yp)

   hp = (gas1.h + mratio * gas2.h) / (1 + mratio) #Enthalpy of product by law of mixtures
   set_hP!(gas_prod, hp, gas1.P) #set gas at correct temperature and pressure

   return gas_prod

end