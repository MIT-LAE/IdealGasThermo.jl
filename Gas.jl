using NLsolve
using LinearAlgebra
using BenchmarkTools
using StaticArrays

include("readThermo.jl")

const Runiv = 8.3145 # J/K/mol
const Pstd = 101325.0 # Pa
const Tstd = 298.15 # K


mutable struct Gas
   P::Float64 # [Pa]
   T::Float64 # [K]
   Tarray::MVector{8, Float64} # Temperature array to make calcs allocation free

   cp::Float64 #[J/mol/K]
   h::Float64  #[J/mol]
   s::Float64  #[J/mol/K]
   Y::MVector{length(spdict), Float64} # Mass fraction of species
   MW::Float64 # Molecular weight [g/mol]
end

# Convinence constructors:
function Gas(T, Y)
   Gas(T, Tarray(T), Y)
end
function Gas(Y)
   Gas(Pstd, Tstd, Tarray(Tstd), 0.0, 0.0, 0.0, Y, 28.965)
end

"""
Constructor that returns a `Gas` type representing 
Air at standard conditions
"""
function Gas()
   Air = spdict[findfirst(x->x=="Air", spdict.name)]

   Gas(Pstd, Tstd, Tarray(Tstd),
    Cp(Tstd, Air), 
    h(Tstd, Air),
    s(Tstd, Pstd, Air),
   [0.0, 0.0, 0.0, 0.0, 0.0, 1.0], Air.MW)

end

# Automatically calculates the Tarray if T is set
function Base.setproperty!(gas::Gas, s::Symbol, val)
   if s === :T
      setfield!(gas, :T, val) # first set T
      setfield!(gas, :Tarray, Tarray!(val, getfield(gas, :Tarray))) # update Tarray
      TT = getfield(gas, :Tarray) # Just convinence
      # Next set the cp, h and s of the gas
      ## Get the right coefficients (assumes Tmid is always 1000.0. Check performed in readThermo.jl.):
      if val<1000.0
         A = view(spdict.alow, :)
      else
         A = view(spdict.ahigh, :)
      end   
      ## Initialize temporary vars
      cptemp = 0.0
      htemp  = 0.0
      stemp  = 0.0
      
      P = getfield(gas, :P)
      Y = getfield(gas, :Y)
      # Go through every species where mass fraction is not zero
      @views for (Yᵢ,a) in zip(Y, A)
         if Yᵢ != 0.0
            cptemp = cptemp + Yᵢ * Cp(TT, a)
             htemp = htemp  + Yᵢ * h(TT, a)
             stemp = stemp  + Yᵢ * (𝜙(TT, a) - Runiv*log(P/Pstd))
         end
      end
   
      setfield!(gas, :cp, cptemp)
      setfield!(gas, :h, htemp)
      setfield!(gas, :s, stemp)

   elseif s === :P
      setfield!(gas, :P, val)
      TT = view(getfield(gas, :Tarray), :) # Just convinence
      # Next set s of the gas
      ## Get the right coefficients (assumes Tmid is always 1000.0. Check performed in readThermo.jl.):
      if val<1000.0
         A = view(spdict.alow, :)
      else
         A = view(spdict.ahigh, :)
      end   
      ## Initialize temporary vars
      stemp  = 0.0
      
      P = val
      Y = view(getfield(gas, :Y), :)
      # Go through every species where mass fraction is not zero
      @views for (Yᵢ,a) in zip(Y, A)
         if Yᵢ != 0.0
            stemp = stemp  + Yᵢ * (𝜙(TT, a) - Runiv*log(P/Pstd))
         end
      end

      setfield!(gas, :s, stemp)

   elseif s === :Y # directly set mass fractions Y
      if typeof(val) === Array{Float64, 1}
         # If array directly store in Y
         setfield!(gas, :Y, val) 
      elseif typeof(val) <: Dict
         # If dict provided set each species in the right order
         names = spdict.name
         Y = zeros(MVector{length(names)})
         for (key,value) in val
            index = findfirst(x->x==key, names)
            Y[index] = value
         end
         setfield!(gas, :Y, Y)
      end
      # Update the MW of the gas mixture
      setfield!(gas, :MW, MW(gas))
   end
   # Note: intentionally not including other variables to prevent users from trying to directly set h, s, cp, MW etc.

end

#Read species data from thermo.inp
spd = readThermo("thermo.inp")

Air = Gas()


ϵ = 1e-12

"""
Function to create the required temperature array
"""
function Tarray(T)
   return [T^-2, T^-1, 1.0, T, T^2, T^3, T^4, log(T)]
end

# function Tarray2(T)
#    TT = zeros(Float64, 8)
#    TT[1:7] = [T^i for i in range(-2, stop=4)]
#    TT[8] = log(T)
#    return TT
# end
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
Calculates cp of the given species in J/K/mol
(This is a completely non-allocating operation.)
"""
@views function cp(Tarray, a)
   #  Cp_R = dot(view(a, 1:7), view(Tarray, 1:7))
    Cp_R = dot(a[1:7], Tarray[1:7])
    Cp = Cp_R*Runiv
    return Cp #J/K/mol
end
"""
Calculates cp of a mixture specified by the mass fraction in `gas`
"""
@views function cp(T, g::Gas)
   Cp = 0.0
   g.T = T
   if T<1000.0
      s = :alow
   else
      s = :ahigh
   end
   
   for (key,Yi) in g.Y
      a = getfield(spd[key], s)
      Cp = Cp + Yi * cp(g.Tarray[1:7], a[1:7])
   end
   return Cp
end

function cp(g::Gas)
   cp(g.T, g)
end

"""
Calculates mean molecular weight
"""
function MW(g::gas)
   MW = 0
   for (key,val) in g.Y
      MW = MW + val*spd[key].MW
   end
   return MW/1000
end
"""
Calculates mean molecular weight
"""
function MW(g::species)
   MW = g.MW
   return MW/1000
end

"""
Calculates h of the given **species** in J/mol
Calcualted by:
H0/RT = -a1*T^-2 + a2*T^-1*ln(T) + a3 + a4*T/2 + a5*T^2/3 + a6*T^3/4 + a7*T^4/5 + b1/T
      = -a1*T₁   + a2*T₂*T₈      + a3 + a4*T₄/2 + a5*T₅/3  + a6*T₆/4  + a7*T₇/5  + a₈*T₂
"""
function h(TT, a)
    h_RT  = -a[1]*TT[1] + 
             a[2]*TT[8]*TT[2] + 
             a[3] + 
         0.5*a[4]*TT[4] + 
             a[5]*TT[5]/3.0 + 
        0.25*a[6]*TT[6] + 
        0.20*a[7]*TT[7] + 
             a[8]*TT[2]

    h = h_RT*TT[4]*Runiv
    return h #J/mol
end
"""
Calculates h of a given **mixture** in J/mol
"""
function h(T, g::Gas)
   H = 0.0
   g.T = T
   if T<1000.0
      s = :alow
   else
      s = :ahigh
   end
   
   for (key,Yᵢ) in g.Y
      a = getfield(spd[key], s)
      H = H + Yᵢ * h(g.Tarray, a)
   end
   return H
end
function h(g::Gas)
   h(g.T,g)
end

"""
Calculates the entropy complement function 𝜙=∫(cₚ/T)dT of the given **species** in J/K/mol
This is calculated at standard state. Tref = 298.15 K, Pref = 101325 Pa.

S0/R = -a1*T^-2/2 - a2*T^-1 + a3*ln(T) + a4*T + a5*T^2/2 + a6*T^3/3.0 + a7*T^4/4 + b2 
     = -a1*T₁/2   - a2*T₂   + a3*T₈    + a4*T₄+ a5*T₅/2  + a6*T₆/3.0  + a7*T₇/4  + a₉   
"""
function 𝜙(TT,a)
    so_R = -0.5*a[1] * TT[1] - 
                a[2] * TT[2] + 
                a[3] * TT[8] + 
                a[4] * TT[4] + 
            0.5*a[5] * TT[5] +
                a[6] * TT[6]/3.0 + 
           0.25*a[7] * TT[7] + 
                a[9]

    so = so_R*Runiv
    return so #J/K/mol
end
"""
Calculates the entropy complement function 𝜙=∫(cₚ/T)dT of the given **mixture** in J/K/mol
This is calculated at standard state. Tref = 298.15 K, Pref = 101325 Pa.
"""
function 𝜙(T, g::Gas)
   S = 0.0
   g.T = T
   if T<1000.0
      s = :alow
   else
      s = :ahigh
   end

   for (key,Yᵢ) in g.Y
      a = getfield(spd[key], s)
      S = S + Yᵢ * 𝜙(g.Tarray, a)
   end
   return S
end
function 𝜙(g::Gas)
   𝜙(g.T, g)
end


"""
Returns standard state sᵒ based on the reference point defined at
Tref = 298.15 K
Pref = 101325 Pa

using the entropy complement function and
the entropy change due to pressure.
Δs can then be defined as sᵒ - sᵒ(Tref, Pref) = sᴼ - 𝜙
"""
function s(T, P, gas::Gas)
   Pref = 101325 
   gas.T = T
   sᵒ =  𝜙(gas) - Runiv*log(P/Pref)
   return sᵒ
end



# Specific functions for gas Compression
PR = 10
p2, T2 = 101325, 298.15
ηp = 0.90

""" 
Adiabatic compression given the 
compression pressure ratio (`PR`), the initial pressure (`p`)
and initial temperature (`T`).

Returns `Tfinal` and `pfinal`
"""
function compress(PR, p, T)
   Tfinal = T * PR^(ℜ/cp(T,Air))

   for i in 1:10
      Res  = (𝜙(Tfinal, Air) - 𝜙(T, Air))/ℜ - log(PR)
      Res′ = cp(Tfinal,Air)/ℜ/Tfinal
      dT  = Res/Res′
      Tfinal = Tfinal - dT
      # println(Tfinal)
      if abs(dT) < ϵ
         break
      end
   end

   return Tfinal, p*PR

end
"""
Adiabatic with NL solve
i.e. find x such that F(x)=0
"""
T = 298.15
p = 101325.
PR = 2.0
function f(x)
   s(T,p,Air) - s(x[1],p*PR,Air)
end

"""
Compression with polytropic efficiency
"""
function compress(PR, p, T, ηp )
   Tfinal = T * PR^(ℜ/cp(T,Air)/ηp)

   for i in 1:25
      Res  = (𝜙(Tfinal, Air) - 𝜙(T, Air))/ℜ - log(PR)/ηp
      Res′ = cp(Tfinal,Air)/ℜ/Tfinal
      dT  = Res/Res′
      
      # ω = 1.0
      # if i>10
      #    ω = 0.5
      # end
      # Tfinal = Tfinal - ω*dT
      Tfinal = Tfinal - dT
      # println("$i: $Tfinal $dT")
      if abs(dT) < ϵ
         break
      end
   end

   return Tfinal, p*PR

end

using PyCall
pygui()
using PyPlot
pygui(true)
plt.style.use(["~/prash.mplstyle", "seaborn-colorblind"])

function plot_polycomp()
   p2 = 101.325e3
   T2 = 288.15
   fig, ax = plt.subplots(figsize = (8,5), dpi = 200)
   PR = LinRange(20,50,10)
   for eff in [0.85, 0.90, 0.91, 0.95]
      T = [T for (T,p) in compress.(PR, p2, T2, eff)]
      l, = ax.plot(PR, T); label = "\$ \\eta_p = $(eff*100)\$%"
      ax.text(PR[end], T[end], label, color = l.get_color())
   end
   
   model_point = ax.scatter(32.64, 827.25,label = "MIT CFM56-5B NPSS model",
   color = "k", marker="*", s=80, zorder=10.1)

   GA_model = ax.scatter(28.6, 798.2, label = "GT CFM56-7B model", color = "k",
                      marker = ".", s = 50, zorder = 10.1)
   
   Engs = ["CFM56-7B27","CFM56-5B3", "LEAP-1B28", "PW1133G"]
   PR = [29.0,32.6, 42.0, 38.67]
   annotate_props = Dict("xycoords"=>("data", "axes fraction"), 
                      "arrowprops"=>Dict("arrowstyle"=>"->","connectionstyle"=>"arc3"), 
                      "ha"=>"center", "va" => "center","fontsize" => 8)
                     #  "bbox"=Dict("fc"=>"w", "ec"=>"none")

   # for (pr,eng) in zip(PR, Engs)
   #    ax.annotate(eng, xy=(pr,0), xytext=(pr,0.08),  zorder = 10)
   #    ax.axvline(pr, color = "k", lw = 1.0, alpha = 0.3, ls = "--", zorder = 9)
   # end
   ax.legend()
   ax.set_ylim(691.2326281862709, 1091.441513548923)

   ax.set_xlabel("\$\\pi_{\\mathrm{oo}}\$")
   ax.set_ylabel("T\$_{t3}\$[K]")
   ax.set_title("\$T_{t3}\$ vs \$\\pi_{oo}\$ at SLS conditions (\$T_{amb}\$ = 288.15 K)")
   plt.tight_layout()
end