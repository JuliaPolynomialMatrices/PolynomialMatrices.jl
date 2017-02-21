"""
    triang(p, iterative::Bool=true, dᵤ) -> U

Polynomial Matrix triangularization based on Sylvester matrix [method 3][1].
The method returns unimodular U such that `p U = L`, where `L` is lower triangular.

By default the iterative [version 3'][1] is used. If `iterative` is to `false`,
the [method 3][1] tries to triangularize `p` with a reduction matrix of degree dᵤ.
By default `dᵤ` is set high enough to guarantee triangularization.

Note that not necessarily the hermite form is returned.

# Examples
```julia
julia> s = variable("s")
p = PolyMatrix([s-1 s^2-1; 2 2s+2; 0 3])
U = triang(p)
3x1 Array{Int64,2}:
  Poly(-0.816497 + 0.408248⋅s)  Poly(-0.57735 - 0.57735⋅s)
  Poly(-0.408248)               Poly(0.57735)
```

# References

-  [1]: D. Henrion, M. Sebek "Reliable Numerical Methods for Polynomial Matrix
        Triangularization" IEEE Transactions on Automatic Control, vol. 44,
        no. 3, Mar. 1999.
"""
function triang{T1,M,V,N}(p::PolyMatrix{T1,M,Var{V},N}, iterative::Bool=true, dᵤ::Int=-1)
  U,L = _triang(p, iterative, dᵤ)
  n,m = size(p)
  return PolyMatrix(U, (m,m), V; reverse=true), PolyMatrix(L, (n,m), V; reverse=true)
end

function _triang{T1,M,V,N}(p::PolyMatrix{T1,M,Var{V},N}, iterative::Bool=true, dᵤ::Int=-1)
  # allow user defined dᵤ
  if !iterative && dᵤ < 0
    dᵤ = _mindegree(p)
  elseif dᵤ < 0
    dᵤ = 0
  end
  d   = degree(p)+dᵤ
  n,m = size(p)

  # tolerance for extracting ϵshape
  ϵ = sqrt(n)*sqrt(m)*d*1e-16
  T = float(T1)
  # construct row permuted sylvester matrix
  Rd = zeros(T, n*(d+1), m*(dᵤ+1))
  for (k,v) in coeffs(p)
    for i in 0:n-1, j in 0:dᵤ
      Rd[n*(d+1)-i*(d+1)-k-dᵤ+j, j*m+(1:m)] = v[n-i,:]
    end
  end
  l,q = lq(Rd)
  U   = q.'

  Σb = zeros(Int,m*(dᵤ+1))
  for i in 1:m*(dᵤ+1)
    # TODO if index is zero p is not of full rank.
    # Should we call the method with Identity appended?
    Σb[i] = findfirst(x->abs(x) > ϵ, l[:,i])
    l[1:Σb[i]-1,i] = zeros(T, Σb[i]-1)
  end

  C = [zeros(Int,0) for i in 1:n]
  # calculate index sets
  for k in eachindex(Σb)
    σ    = Σb[k]
    it,r = divrem(σ, d+1)
    i    = r > 0 ? it+1 : it
    push!(C[i], k)
  end

  Σ = zeros(Int,0)
  for i in 1:m*(dᵤ+1)
    if Σb[i] > n
      break
    end
    if !isempty(C[Σb[i]])
      push!(Σ,i)
    end
  end

  if length(Σ) < m
    if iterative
      return _triang(p, iterative, dᵤ+1)
    else
      throw(ErrorException("triang: failed to triangularize"))
    end
  else
    k = [maximum(C[i]) for i in Σ]
    L = zeros((d+1)*n,m)
    for i in 0:d
      for j in 0:n-1
        L[(d-i)*n+j+1,:] = l[(j)*(d+1)+d-i+1,k]
      end
    end
    return U[:,k], L
  end
end

function _mindegree(p::PolyMatrix)
  m = minimum(size(p))
  rowdegs = sort(_rowdegrees(p))
  coldegs = sort(_coldegrees(p))
  return min(sum(rowdegs[end-m+2:end]), sum(coldegs[end-m+2:end]))
end

_rowdegrees(p::PolyMatrix) = [degree(p[i,:]) for i in 1:size(p,1)]
_coldegrees(p::PolyMatrix) = [degree(p[:,i]) for i in 1:size(p,2)]
