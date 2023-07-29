# run some tests w/ PerOrPol

import Test:
  @test, @testset
import Polynomials: 
  Polynomial,derivative
import PerOrPol: 
  Newton_p2, Newton_p3, Newton_p

# Newton step
rstep(x,f,df)=x-f(x)//df(x)

@testset "Newton_p2" begin 
  for t in 1:10
    ret=Newton_p2()
    pol=Polynomial(ret.pol)
    dpol=derivative(pol)
    x1,x2=ret.orbit
    @test rstep(x1,pol,dpol)==x2 && rstep(x2,pol,dpol)==x1
  end 
end


@testset "Newton_p3" begin 
  for t in 1:10
    ret=Newton_p3()
    pol=Polynomial(ret.pol)
    dpol=derivative(pol)
    x1,x2,x3=ret.orbit
    @test rstep(x1,pol,dpol)==x2 && rstep(x2,pol,dpol)==x3 && rstep(x3,pol,dpol)==x1
  end 
end


fstep(x,f,df)=x-f(x)/df(x)
setprecision(BigFloat,128)
@testset "Newton_p(4:8)" begin 
  for p in 4:8
    ret=Newton_p(p,(lhs_type=BigFloat,))
    pol=Polynomial(ret.pol)
    dpol=derivative(pol)
    x=ret.orbit
    x=[x...,x[1]]
    @test all([isapprox(fstep(x[k],pol,dpol),x[k+1]) for k in 1:p])
  end 
end

