using ForwardDiff

function FDjacobian(fun, x0)
    delta = 1e-6 #Finite diff perturbation
    f0 = fun(x0)
    jac_FD = zeros(length(f0), length(x0))
    for i in 1:length(x0)
        x = deepcopy(x0)
        x[i] = x[i]*(1+delta)
        fi = fun(x)
        jac_FD[:, i] .= (fi.-f0)./(x[i]-x0[i])
    end
    return jac_FD
end

@testset "thermo properties" begin 

    #Simple equation that returns enthalpy
    function h(x)
        T = x[1]
        P = x[2]
        
        g = Gas(T, P)
        h = g.h
        return h
    end

    x0 = [288.15, 101325.0]
    grad = ForwardDiff.gradient(h, x0)

    #Compare gradient to exact solution
    g = Gas(x0[1], x0[2])
    grad_exact = [g.cp, 0.0]
    for (i,g) in enumerate(grad)
        @test g ≈ grad_exact[i] rtol = 1e-8
    end

    #Elaborate function of thermodynamic parameters
    function thermofun(x)
        T = x[1]
        P = x[2]
        
        g = Gas(T, P)
        s = g.s
        h = g.h
        return s^2/h
    end

    x0 = [288.15, 101325.0]
    grad = ForwardDiff.gradient(thermofun, x0)

    #Compare gradient to finite differences
    grad_FD = FDjacobian(thermofun, x0)
    for (i,g) in enumerate(grad)
        @test g ≈ grad_FD[i] rtol = 1e-4
    end
end

@testset "combustion" begin 
    #Function that returns FAR
    function FARfun(x)
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
            Tb, etab, hvap)
        return FAR
    end

    x0 = [298.15, 101325.0, 298.15, 1000.0, 0.99, 1e5]
    grad = ForwardDiff.gradient(FARfun, x0)

    #Compare gradient to finite differences
    grad_FD = FDjacobian(FARfun, x0)
    for (i,g) in enumerate(grad)
        @test g ≈ grad_FD[i] rtol = 1e-4
    end

    #Function that returns the product temperature
    function Tflame(x)
        T = x[1]
        P = x[2]
        Tf = x[3]
        FAR = x[4]
        etab = x[5]
        hvap = x[6]

        g = Gas(T, P)
        g_prod = IdealGasThermo.fuel_combustion(g,
            "CH4", Tf, FAR,
            etab, hvap) 
        return g_prod.T
    end

    x0 = [298.15, 101325.0, 298.15, 0.01, 0.99, 1e5]
    grad = ForwardDiff.gradient(Tflame, x0)

    #Compare gradient to finite differences
    grad_FD = FDjacobian(Tflame, x0)
    for (i,g) in enumerate(grad)
        @test g ≈ grad_FD[i] rtol = 1e-4
    end
end

@testset "turbo" begin 
    function PRfun(x)
        T = x[1]
        P = x[2]
        PR = x[3]
        eta = x[4]

        gas = Gas(T, P)

        IdealGasThermo.PressureRatio!(gas, PR, eta)
        return [gas.T, gas.P]
    end

    x0 = [298.15, 101325.0, 10.0, 0.9]
    J = ForwardDiff.jacobian(PRfun, x0)

    J_FD = FDjacobian(PRfun, x0)
    for i in 1:2
        for j in 1:length(x0)
            @test J[i, j] ≈ J_FD[i, j] rtol = 1e-4
        end
    end
end
