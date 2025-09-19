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
@testset "autodiff" begin 

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

        function Mach(x)
            T = x[1]
            P = x[2]
            M0 = x[3]
            M = x[4]

            gas = Gas(T, P)
            IdealGasThermo.gas_Mach!(gas, M0, M)
            return [gas.T, gas.P]
        end
        x0 = [298.15, 101325.0, 0.1, 0.8]
        J = ForwardDiff.jacobian(Mach, x0)

        J_FD = FDjacobian(Mach, x0)
        for i in 1:2
            for j in 1:length(x0)
                @test J[i, j] ≈ J_FD[i, j] rtol = 1e-4
            end
        end
    end

    @testset "Gas1D" begin 
        #Elaborate function of thermodynamic parameters
        function thermofun1D(x)
            T = x[1]
            P = x[2]
            
            g = Gas1D(T, P)
            s = g.s
            h = g.h
            return s^2/h
        end
        x0 = [288.15, 101325.0]
        grad = ForwardDiff.gradient(thermofun1D, x0)

        #Compare gradient to finite differences
        grad_FD = FDjacobian(thermofun1D, x0)
        for (i,g) in enumerate(grad)
            @test g ≈ grad_FD[i] rtol = 1e-4
        end
    end

    @testset "atmos" begin 
        #Returns temperature and pressure at a given altitude
        function atmosfun(x)
            z = x[1]
            
            g = IdealGasThermo.standard_atmosphere(z)
            return [g.T, g.P]
        end
        x0 = [1e3]
        J = ForwardDiff.jacobian(atmosfun, x0)

        #Compare gradient to finite differences
        J_FD = FDjacobian(atmosfun, x0)
        for i in 1:2
            for j in 1:length(x0)
                @test J[i, j] ≈ J_FD[i, j] rtol = 1e-4
            end
        end
    end

    @testset "Gas1D" begin 
        #Elaborate function of thermodynamic parameters
        function thermofun1D(x)
            T = x[1]
            P = x[2]
            
            g = Gas1D(T, P)
            s = g.s
            h = g.h
            return s^2/h
        end
        x0 = [288.15, 101325.0]
        grad = ForwardDiff.gradient(thermofun1D, x0)

        #Compare gradient to finite differences
        grad_FD = FDjacobian(thermofun1D, x0)
        for (i,g) in enumerate(grad)
            @test g ≈ grad_FD[i] rtol = 1e-4
        end
    end

    @testset "engine" begin 
        #Basic Brayton thermodynamic cycle to test derivatives
        function JetEngineResiduals(x, p)
            PR_t = x[1]
            mdot_c = x[2]
            R = typeof(PR_t)

            z0 = p[1]
            M0 = p[2]
            PR = R(p[3])
            Tt4 = R(p[4])
            F_N = p[5]
            
            #Intake
            g0 = IdealGasThermo.standard_atmosphere(z0)
            g0 = Gas(R(g0.T), R(g0.P))
            u0 = M0*sqrt(g0.T*g0.gamma*g0.R)
            gt2 = deepcopy(g0) #Copy so as to not modify it
            IdealGasThermo.gas_Mach!(gt2, M0, 0.0) #stagnation properties
            #Compressor
            gt3 = deepcopy(gt2) #Copy so as to not modify it
            IdealGasThermo.compress!(gt3, PR)
            #Burner
            FAR,_ = IdealGasThermo.gas_burn(gt3,
                "CH4",
                R(298.15),
                Tt4, 1.0, 0.0)
            gt4 = IdealGasThermo.fuel_combustion(gt3,
                "CH4", R(298.15), FAR,
                1.0, 0.0)
            #Turbine
            gt5 = deepcopy(gt4) #Make a copy to avoid modifying it
            IdealGasThermo.expand!(gt5, PR_t) #Expand in turbine
            #Nozzle
            g6 = deepcopy(gt5) #Make a copy to avoid modifying it
            PR_n = gt2.P / gt5.P #Nozzle pressure ratio
            IdealGasThermo.expand!(g6, PR_n)

            #Residuals
            u = sqrt(2 * (gt5.h - g6.h))
            Fsp = (1 + FAR)*u - u0
            
            return [F_N - mdot_c*Fsp; gt3.h - gt2.h - (1+FAR)*(gt4.h - gt5.h)] #Thrust and shaft power residuals
        end
        x0 = [0.5, 1.0]
        p = [1e4, 0.8, 10.0, 1700.0, 1e4]
        f(x) = JetEngineResiduals(x, p)
        J = ForwardDiff.jacobian(f, x0)

        #Compare gradient to finite differences
        J_FD = FDjacobian(f, x0)
        for i in 1:2
            for j in 1:length(x0)
                @test J[i, j] ≈ J_FD[i, j] rtol = 1e-4
            end
        end

        #Solve problem for PR_t and mdot_c
        eps = 1e-6
        f0 = [1.0, 1.0]
        x = x0
        while maximum(abs.(f0)) > eps #Use a basic Newton's method
            f0 = f(x)
            J = ForwardDiff.jacobian(f, x)
            dx = -J\f0
            x = x + dx
        end
        x_check = [0.6183259629412734, 10.459646582243932] #Validated against a basic turbojet model
        for (i,xi) in enumerate(x)
            @test xi ≈ x_check[i] rtol = 1e-6
        end
    end
end