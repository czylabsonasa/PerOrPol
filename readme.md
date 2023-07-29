### Periodic Orbit Polynomials

* `Newton_p2()`: generates a random quadratic polynomial with length 2 periodic orbit.
* `Newton_p3()`: generates a random cubic polynomial with length 3 periodic orbit.

* usage:
```julia
import Pkg
Pkg.activate(;temp=true)
Pkg.add(;url="https://github.com/czylabsonasa/PerOrPol")
Pkg.add("Polynomials") # for testing
import PerOrPol: Newton_p2, Newton_p3
import Polynomials: Polynomial, derivative

ret=Newton_p2()
p=Polynomial(ret.pol)
dp=derivative(p)
step(x)=x-p(x)//dp(x)
x1,x2=ret.orbit
@assert step(x1)==x2 && step(x2)==x1


ret=Newton_p3()
p=Polynomial(ret.pol)
dp=derivative(p)
step(x)=x-p(x)//dp(x)
x1,x2,x3=ret.orbit
@assert step(x1)==x2 && step(x2)==x3 && step(x3)==x1

```
* This is my 1st "package" in julia - still try to learn the stuff. The functions (in different form) 
was used in generating moodle-quizzes for my numerical mathematics courses.
