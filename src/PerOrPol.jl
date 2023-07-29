"""
### generate rational coefficient periodic orbit polynomials
"""
module PerOrPol
  using StatsBase: sample
  using LinearAlgebra: det

  const _pool_lo=-10
  const _pool_up=10
  const _co_lo=-30
  const _co_up=30

"""
    Newton2(;<keyword arguments>)
  it generates a random quadratic poly with Newton-period of length 2 and
  returns a namedtuple `(pol=[c,b,a],x=[x1,x2])`, where the required polynomial is ax^2+bx+c and 
  Newton's method generates: `x1->x2->x1...`.

# Arguments:
* they are lower and upper bounds related to the random selection of the numbers used in the process
"""
  function Newton_p2(
    ;
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
      # based on the condition one can set up a 2x2 linear system w/ unknowns a,b and 
      # [c,c]' rigth hand side - the coeff matrix is:
      LHS=[x1*(x1-2*x2) -x2; x2*(x2-2*x1) -x1]
      iszero(det(LHS)) && continue
      c=rand(pool)
      RHS=[c; c]
      a,b=LHS\RHS
      all(co_lo .≤ [a.num,a.den,b.num,b.den] .≤ co_up) && break
    end
    (pol=[c,b,a],orbit=[x1,x2])
  end

"""  
    Newton3(;<keyword arguments>)

  it generates a random quadratic poly with Newton-period of length 3 and
  returns a namedtuple `(pol=[d,c,b,a],x=[x1,x2,x3])`, where the required polynomial is `ax^3+bx^2+cx+d`, and 
  Newton's method generates the orbit: `x1->x2->x3->x1...`

# Arguments:
* they are lower and upper bounds related to the random selection of the numbers used in the process
"""
  function Newton_p3(
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

    genrow(x,y)=[2*x^3-3*x^2*y, x^2-2*x*y, -y]'
    x1,x2,x3,a,b,c,d=fill(0//1,7)
    while true
      x1,x2,x3=sample(pool, 3, replace=false)
      LHS=vcat(genrow(x1,x2),genrow(x2,x3),genrow(x3,x1))
      iszero(det(LHS)) && continue
      d=rand(pool)
      RHS=[d; d; d]
      a,b,c=LHS\RHS
      all(co_lo .≤ [a.num,a.den,b.num,b.den,c.num,c.den] .≤ co_up) && break
    end
    (pol=[d,c,b,a],orbit=[x1,x2,x3])
  end
  
  

end # module PerOrPol
