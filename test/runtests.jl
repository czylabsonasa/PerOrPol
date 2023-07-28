# run some tests w/ PerOrPol

import Test:
  @test, @testset
import Polynomials: 
  Polynomial,derivative
import PerOrPol: 
  Newton2, Newton3

# Newton step
step(x,f,df)=x-f(x)//df(x)

@testset "Newton2" begin 
  for t in 1:10
    ret=Newton2()
    p=Polynomial(ret.pol[3:-1:1]) # Polynomial expects coeffs in increasing order
    dp=derivative(p)
    x1,x2=ret.orbit
    @test step(x1,p,dp)==x2 && step(x2,p,dp)==x1
  end 
end


@testset "Newton3" begin 
  for t in 1:10
    ret=Newton3()
    p=Polynomial(ret.pol[4:-1:1]) # Polynomial expects coeffs in increasing order
    dp=derivative(p)
    x1,x2,x3=ret.orbit
    @test step(x1,p,dp)==x2 && step(x2,p,dp)==x3 && step(x3,p,dp)==x1
  end 
end

