"""
### generate periodic orbit polynomials
"""
module PerOrPol
  using StatsBase: sample
  using LinearAlgebra: det

"""
    _def_pars::NamedTuple
* default bounds/types used by the functions
* currently `nums, dens, coeffs, lhs_type` fields
"""  
  _def_pars=(
    nums=(-10,10),
    dens=(1,2,4,8),
    coeffs=(-30,30),
    lhs_type=Rational{Int},
  )

  function adjust_pars(pars)
    nums=get(pars,:nums,_def_pars.nums)
    dens=get(pars,:dens,_def_pars.dens)
    coeffs=get(pars,:coeffs,_def_pars.coeffs)
    lhs_type=get(pars,:lhs_type,_def_pars.lhs_type)

    (nums=nums,dens=dens,coeffs=coeffs,lhs_type=lhs_type)
  end
  

"""
    Newton_p2(pars=_def_pars)
  it generates a random quadratic poly with Newton-period of length 2 and
  returns a namedtuple `(pol=[c,b,a],x=[x1,x2])`, where the required polynomial is ax^2+bx+c and 
  Newton's method generates: `x1->x2->x1...`.

# Arguments:
* through the `pars` namedtuple we can refine the method
  * the `nums` and `dens` is a lower and upper bounds for the numerators and denominators
  * the `coeffs` contains an lower and upper bound for the coeffs of the resulting polynomial
  * lhs_type controls the type of the coefficient matrix in the linear system
"""
  function Newton_p2(pars=_def_pars)
    pars=adjust_pars(pars)
    nums=pars.nums
    dens=pars.dens
    coeffs=pars.coeffs
    lhs_type=pars.lhs_type
    

    
    pre=collect(nums[1]:nums[2])
    pool=Rational{Int}[]
    for den in dens
      pool=vcat(pool,pre//den)
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
      a,b=lhs_type.(LHS)\RHS
      (lhs_type===Float64) || all(coeffs[1] .≤ [a.num,a.den,b.num,b.den] .≤ coeffs[2]) && break
    end
    (pol=lhs_type[c,b,a],orbit=[x1,x2])
  end

"""  
    Newton_p3(pars=_def_pars)

  it generates a random quadratic poly with Newton-period of length 3 and
  returns a namedtuple `(pol=[d,c,b,a],x=[x1,x2,x3])`, where the required polynomial is `ax^3+bx^2+cx+d`, and 
  Newton's method generates the orbit: `x1->x2->x3->x1...`

# Arguments:
* through the `pars` namedtuple we can refine the method
  * the `nums` and `dens` is a lower and upper bounds for the numerators and denominators
  * the `coeffs` contains an lower and upper bound for the coeffs of the resulting polynomial
  * lhs_type controls the type of the coefficient matrix in the linear system
"""
  function Newton_p3(pars=_def_pars)
    pars=adjust_pars(pars)
    
    nums=pars.nums
    dens=pars.dens
    coeffs=pars.coeffs
    lhs_type=pars.lhs_type
    

    
    pre=collect(nums[1]:nums[2])
    pool=Rational{Int}[]
    for den in dens
      pool=vcat(pool,pre//den)
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
      a,b,c=lhs_type.(LHS)\RHS
      (lhs_type===Float64 || all(coeffs[1] .≤ [a.num,a.den,b.num,b.den,c.num,c.den] .≤ coeffs[2])) && break
    end
    (pol=lhs_type[d,c,b,a],orbit=[x1,x2,x3])
  end
  


"""  
    Newton_p(p::Int,pars=_def_pars)

  it generates a random degree of `p` poly with Newton-period of length `p` and
  returns a namedtuple `(pol=a,x=x)`, where the required polynomial is `a[0]+a[1]x+...a[p]x^p`, and 
  Newton's method generates the orbit: `x1->x2->...->xp->x1...`

# Arguments:
* `p` is the length of the orbit.
* through the `pars` namedtuple we can refine the method
  * the `nums` and `dens` is a lower and upper bounds for the numerators and denominators
  * the `coeffs` contains an lower and upper bound for the coeffs of the resulting polynomial
  * lhs_type controls the type of the coefficient matrix in the linear system
"""
  function Newton_p(p::Int, pars=_def_pars)
    if p<2
      error("not implemented for p=$(p)\n")
      return
    end
    pars=adjust_pars(pars)
    # these can be handled w/ the general way too
    if p==2
      return Newton_p2(pars)
    end
    if p==3
      return Newton_p3(pars)
    end

    nums=pars.nums
    dens=pars.dens
    coeffs=pars.coeffs
    lhs_type=pars.lhs_type
    

    
    pre=collect(nums[1]:nums[2])
    pool=Rational{Int}[]
    for den in dens
      pool=vcat(pool,pre//den)
    end
    pool=setdiff(pool, 0) # all coeff and the orbit points are nonzero (ad-hoc condition)

    genrow(x,y)=((xx,yy)=lhs_type.((x,y));[xx^(k-1)*((k-1)*xx-k*yy) for k in 1:p]')
    x=fill(0//1,p)
    a=fill(0//1,p)
    a0=0//1
    while true
      x=sample(pool, p, replace=false)
      LHS=vcat(vcat([genrow(x[k],x[k+1]) for k in 1:p-1]...),genrow(x[p],x[1]))
      tp=promote_type(Float64,lhs_type)
      (abs(det(tp.(LHS)))<1e-9) && continue
      a0=rand(pool)
      RHS=fill(a0,p)
      a=LHS\RHS
      (lhs_type<:AbstractFloat || all(coeffs[1] .≤ vcat([t.num for t in a],[t.den for t in a]) .≤ coeffs[2])) && break
    end
    (pol=vcat(a0,a),orbit=x)
  end

  

end # module PerOrPol
