const Runiv = 8.3145 # J/K/mol
const Pstd = 101325.0 # Pa
const Tstd = 298.15 # K
const ϵ = sqrt(eps()) #standard tolerance 

# Air composition
function Xair_typed(R)
    Xair = Dict{String, R}(
        "N2" => 0.78084,
        "Ar" => 0.009365,
        "H2O" => 0.0,
        "CO2" => 0.000319,
        "O2" => 0.209476
    )
    return Xair
end

const Xair = Xair_typed(Float64)

const XwetAir = Dict{String, Float64}(

)

const Ytasopt =
    Dict{String, Float64}("N2" => 0.7532, "O2" => 0.2315, "CO2" => 0.0006, "H2O" => 0.0020, "Ar" => 0.0127)
