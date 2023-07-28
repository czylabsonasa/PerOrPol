# run some tests w/ PerOrPol

using Test
using Polynomials: Polynomial,derivative
import PerOrPol: Newton2

printstyled(stderr,"Newton\n",color=:light_yellow)
# Newton step
step(x,f,df)=x-f(x)//df(x)

@testset begin 
  for t in 1:10
    ret=Newton2()
    p=Polynomial(ret.pol[3:-1:1]) # Polynomial expects coeffs in increasing order
    dp=derivative(p)
    x1,x2=ret.x
    @test step(x1,p,dp)==x2 && step(x2,p,dp)==x1
  end 
end
