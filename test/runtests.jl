# run some tests w/ PerOrPol

import Test:
  @test, @testset
import Polynomials: 
  Polynomial,derivative
import PerOrPol: 
  Newton_p2, Newton_p3

# Newton step
step(x,f,df)=x-f(x)//df(x)

@testset "Newton_p2" begin 
  for t in 1:10
    ret=Newton_p2()
    p=Polynomial(ret.pol)
    dp=derivative(p)
    x1,x2=ret.orbit
    @test step(x1,p,dp)==x2 && step(x2,p,dp)==x1
  end 
end


@testset "Newton_p3" begin 
  for t in 1:10
    ret=Newton_p3()
    p=Polynomial(ret.pol)
    dp=derivative(p)
    x1,x2,x3=ret.orbit
    @test step(x1,p,dp)==x2 && step(x2,p,dp)==x3 && step(x3,p,dp)==x1
  end 
end

