# Modeling combustion

`IdealGasThermo` allows us to perform simple combustion calculations. 
The following functions are useful utilites to do so.

## Calculating stoichiometric fuel-oxidizer ratio 

For an arbitrary fuel of type ``\genfuel`` we can write a general equation representing the combustion of 1 mole of fuel,

```math
\genfuel + n_{\mathrm{O}_2} \mathrm{O}_2 
\longrightarrow 
n_{\mathrm{CO}_2} \mathrm{CO}_2 + n_{\mathrm{H}_2\mathrm{O}} \mathrm{H}_2\mathrm{O} 
+ n_{\mathrm{N}_2} \mathrm{N}_2

```

where the number of moles of the different species is given by balancing the reaction:
```math
\begin{aligned}
n_{\rm{CO_2}} &= x_{\mathrm{C}} \tag{1}\\
n_{\rm{H_2O}} &= \frac{x_{\mathrm{H}}}{2}\\
n_{\rm{N_2}} &= \frac{x_{\mathrm{N}}}{2}\\
n_{\rm{O_2}} &= x_{\mathrm{C}} + \frac{x_{\mathrm{H}}}{4} - \frac{x_{\mathrm{O}}}{2}.
\end{aligned}
```

The molar fuel-oxygen ratio ``f`` is then ``\displaystyle{\frac{1}{n_{\rm{O_2}}}}``. If
the oxidzier is not pure oxygen then the stoichiometric molar 
*fuel-oxidizer* ratio is given by 
``\displaystyle{f_{\text{stoich.}}=\frac{X_{\rm{O_2}}}{n_{\rm{O_2}} }}``, 
where ``X_{\rm{O_2}}`` is the mole fraction of oxygen in the oxidizer (e.g., ``X_{\rm{O_2}} \approx 0.21`` in dry air). 


This reaction can also be written as
```math
\genfuel
\longrightarrow 
n_{\mathrm{CO}_2} \mathrm{CO}_2 + n_{\mathrm{H}_2\mathrm{O}} \mathrm{H}_2\mathrm{O} 
+ n_{\mathrm{N}_2} \mathrm{N}_2
- n_{\mathrm{O}_2} \mathrm{O}_2 
```
which can be read as "the complete combustion of 1 mole of fuel
``(\genfuel)`` consumes
``n_{\mathrm{O}_2}`` moles of ``\mathrm{O}_2`` and
produces ``n_{\mathrm{CO}_2}`` moles of ``\mathrm{CO}_2``,
``n_{\mathrm{H}_2\mathrm{O}}`` moles of ``\mathrm{H}_2\mathrm{O}``, and
``n_{\mathrm{N}_2}`` moles of ``\mathrm{N}_2``".

```@docs
IdealGasThermo.fuelbreakdown
IdealGasThermo.stoich_molar_fuel_oxy_ratio
IdealGasThermo.stoich_molar_FOR
IdealGasThermo.stoich_FOR
IdealGasThermo.reaction_change_molar_fraction
```

## [Vitiated gas composition](@id vitiated)

If we consider lean combustion (i.e., more oxygen present than required 
to completely react with the fuel) we can write the above equation for some molar fuel-oxygen (or more generally oxidizer) ratio ``(f \leq f_{\mathrm{stoich.}})`` as follows

```math

\begin{aligned}
f \times
\genfuel+ 
1\times \mathrm{O}_2 
&\longrightarrow 

  f n_{\mathrm{CO}_2} \mathrm{CO}_2 
&+& f n_{\mathrm{H}_2\mathrm{O}} \mathrm{H}_2\mathrm{O} 
&+& f n_{\mathrm{N}_2} \mathrm{N}_2 
&+&\left(1 - fn_{\mathrm{O}_2}\right) \mathrm{O}_2 

\\

f \times
\genfuel + 
1\times \mathrm{O}_2 

&\longrightarrow 

  f x_{\mathrm{C}} \mathrm{CO}_2 
&+& f \frac{x_{\mathrm{H}}}{2} \mathrm{H}_2\mathrm{O} 
&+& f \frac{x_{\mathrm{N}}}{2}\mathrm{N}_2 
&+&\left(1 - \frac{f}{f_{\mathrm{stoich.}}}\right)\mathrm{O}_2.

\end{aligned}
```
where ``\fst = 1/n_{\mathrm{O}_2}`` for oxy-combustion. More generally
for some oxidzer that has ``X_{\rm{O_2}}`` moles of oxygen per mole of oxidizer,

```math
\begin{aligned}
f
\genfuel + 
\mathrm{Oxidizer}

\longrightarrow 

  f x_{\mathrm{C}} \mathrm{CO}_2 
+ f \frac{x_{\mathrm{H}}}{2} \mathrm{H}_2\mathrm{O} 
+ f \frac{x_{\mathrm{N}}}{2}\mathrm{N}_2 
&+\mathrm{Oxidizer}\\
& - \left(\frac{f X_{\mathrm{O}_2}}{\fst}\right)\mathrm{O}_2.
\end{aligned}
```



```@docs
IdealGasThermo.vitiated_mixture
IdealGasThermo.vitiated_species
IdealGasThermo.fixed_fuel_vitiated_species
```