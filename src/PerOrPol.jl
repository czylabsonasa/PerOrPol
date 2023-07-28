"""
### generate periodic orbit polynomials
"""
module PerOrPol
  using StatsBase: sample
  using LinearAlgebra: det

  const _pool_lo=-10
  const _pool_up=10
  const _co_lo=-30
  const _co_up=30

"""  
  it generates a random quadratic poly with Newton period of length 2 
  returns (pol=(a,b,c),x=(x1,x2)), where the required polynomial is ax^2+bx+c, and 
  Newton's method generates: x1->x2->x1...
"""
  function Newton2(
    pool_lo::Int=_pool_lo,
    pool_up::Int=_pool_up,
    co_lo::Int=_co_lo,
    co_up::Int=_co_up
  )
    
    alap=collect(pool_lo:pool_up)
    pool=Rational{Int}[]
    for den in [1,2,4,8]
      pool=vcat(pool,alap//den)
    end
    pool=setdiff(pool, 0) # all coeff and the orbit points are nonzero (ad-hoc condition)

    x1,x2,a,b,c=fill(0//1,5)
    while true
      x1,x2=sample(pool, 2, replace=false)
      LHS=[x1*(x1-2*x2) -x2; x2*(x2-2*x1) -x1]
      iszero(det(LHS)) && continue
      c=rand(pool)
      RHS=[c; c]
      a,b=LHS\RHS
      all(co_lo .≤ [a.num,a.den,b.num,b.den] .≤ co_up) && break
    end
    (pol=(a,b,c),orbity=(x1,x2))
  end
  

end # module PerOrPol
