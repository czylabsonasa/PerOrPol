* `Newton2()`: generates a random quadratic polynomial with length 2 periodic orbit.
* `Newton3()`: generates a random cubic polynomial with length 3 periodic orbit.

* usage:
```julia
import Pkg
Pkg.activate(;temp=true)
Pkg.add(;url="https://github.com/czylabsonasa/PerOrPol")
Pkg.add("Polynomials") # for testing
import PerOrPol: Newton2, Newton3
import Polynomials: Polynomial, derivative

ret=Newton2()
p=Polynomial(ret.pol[3:-1:1])
dp=derivative(p)
step(x)=x-p(x)//dp(x)
x1,x2=ret.orbit
@assert step(x1)==x2 && step(x2)==x1


ret=Newton3()
p=Polynomial(ret.pol[4:-1:1])
dp=derivative(p)
step(x)=x-p(x)//dp(x)
x1,x2,x3=ret.orbit
@assert step(x1)==x2 && step(x2)==x3 && step(x3)==x1

```
  
