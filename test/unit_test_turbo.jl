@testset "turbomachinery" begin
  gas = Gas1D()
  gas.T = Tstd = IdealGases.Tstd
  gas.P = Pstd = IdealGases.Pstd

  IdealGases.compress(gas, 2.0, 1.0)
  @test gas.P == 2*Pstd

  #Test that compress throws right errors
  PR = 0.5
  err_msg = "The specified pressure ratio (PR) to compress by needs to be ≥ 1.0.
      Provided PR = $PR. Did you mean to use `expand`?"
  @test_throws ErrorException(err_msg) IdealGases.compress(gas, PR)
  PR = 1.5
  err_msg = "The specified pressure ratio (PR) to compress by needs to be ≤ 1.0.
      Provided PR = $PR. Did you mean to use `compress`?"
  @test_throws ErrorException(err_msg) IdealGases.expand(gas, PR)

  #Test gas mixing
  set_TP!(gas, Tstd, Pstd)
  gas2 = deepcopy(gas)
  set_TP!(gas2, 3*Tstd, Pstd)
  gas3 = IdealGases.gas_mixing(gas, gas2, 2.0)
  @test gas3.T == 703.6767764998808

  @testset "tasopt comparisons" begin
    gas = Gas()
    gas.Y = IdealGases.Ytasopt
    gas.T = Tstd
    gas.P = Pstd
    @test gas.cp ≈ 1006.5028925107893 rtol = 1e-4
    @test gas.h ≈ -32079.665957469897 rtol = 1e-2

    gas.T = 2000.0
    @test gas.cp ≈ 1253.6089 rtol = 1e-4

    set_h!(gas, 0.0) 
    @test gas.T ≈ 329.99 rtol = 1e-3

    set_TP!(gas, Tstd, Pstd)
    IdealGases.compress(gas, 2.0, 1.0)
    @test gas.T ≈ 363.29 atol=1e-1

    set_TP!(gas, Tstd, Pstd)
    IdealGases.expand(gas, 0.5, 1.0)
    @test gas.T ≈ 244.547 atol=1e-1

    set_TP!(gas, Tstd, Pstd)
    IdealGases.gas_Mach!(gas, 0.0, 1.0, 1.0)
    @test gas.T ≈ 248.41188254462523 rtol = 1e-4

  end

end