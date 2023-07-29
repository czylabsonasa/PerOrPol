### Periodic Orbit Polynomials

* `Newton_p2()`: generates a random quadratic polynomial with length 2 periodic orbit.
* `Newton_p3()`: generates a random cubic polynomial with length 3 periodic orbit.
* `Newton_p(p::Int)`: generates a random degree of `p` polynomial with length `p` periodic orbit.
  * note that, the default approach  - using `Rational{Int}` for solving the linear system - is not appropriate for p≥4, with the default bound settings, bcos of overflows.
  * the `p=4` is amenable with `Rational{BigInt}` (see the `res` dir) and after tweaking the bounds `p=5` can be ok.
  * u can use `Float64/BigFloat`...


* usage:
```julia
import Pkg
Pkg.activate(;temp=true)
Pkg.add(;url="https://github.com/czylabsonasa/PerOrPol")
Pkg.add("Polynomials") # for testing
import PerOrPol: Newton_p2, Newton_p3, Newton_p
import Polynomials: Polynomial, derivative

rstep(x,f,df)=x-f(x)//df(x)
ret=Newton_p2()
pol=Polynomial(ret.pol)
dpol=derivative(pol)
x1,x2=ret.orbit
@assert rstep(x1,pol,dpol)==x2 && rstep(x2,pol,dpol)==x1


ret=Newton_p3()
pol=Polynomial(ret.pol)
dpol=derivative(pol)
x1,x2,x3=ret.orbit
@assert rstep(x1,pol,dpol)==x2 && rstep(x2,pol,dpol)==x3 && rstep(x3,pol,dpol)==x1


fstep(x,f,df)=x-f(x)/df(x)
setprecision(BigFloat,128)
for p in 4:8
  ret=Newton_p(p,(lhs_type=BigFloat,))
  pol=Polynomial(ret.pol)
  dpol=derivative(pol)
  x=ret.orbit
  x=[x...,x[1]]
  @assert all([isapprox(fstep(x[k],pol,dpol),x[k+1]) for k in 1:p])
end 


```
* This is my 1st "package" in julia - still try to learn the stuff. The functions (in different form) 
was used in generating moodle-quizzes for my numerical mathematics courses.
