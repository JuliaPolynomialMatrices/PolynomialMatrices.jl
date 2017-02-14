function PMcheck(p1::PolyMatrix, p2::PolyMatrix)
  @assert 1==1
end

function +{T1,M1,V,N,T2,M2}(p1::PolyMatrix{T1,M1,Var{V},N}, p2::PolyMatrix{T2,M2,Var{V},N})
  size(p1) == size(p2) || error("incompatible sizes")
  PMcheck(p1,p2)

  # figure out return type
  c1     = coeffs(p1)
  c2     = coeffs(p2)
  k1, v1 = first(c1)
  k2, v2 = first(c2)
  vr     = v1+v2
  M      = typeof(vr)
  r      = PolyMatrix( SortedDict{Int,M,ForwardOrdering}(), size(p1), V)

  cr  = coeffs(r)
  sᵢ  = intersect(keys(c1),keys(c2))
  s₁  = setdiff(keys(c1),sᵢ)
  s₂  = setdiff(keys(c2),sᵢ)
  for k in s₁
    insert!(cr, k, c1[k])
  end
  for k in s₂
    insert!(cr, k, c2[k])
  end
  for k in sᵢ
    insert!(cr, k, c1[k]+ c2[k])
  end
  return r
end

function +{T1,M1,V1,V2,N,T2,M2}(p1::PolyMatrix{T1,M1,Var{V1},N}, p2::PolyMatrix{T2,M2,Var{V2},N})
  warn("p1≈p2: `p1` ($T1,$V1) and `p2` ($T2,$V2) have different variables")
  throw(DomainError())
end

+{T,M,V,N}(p1::PolyMatrix{T,M,Var{V},N}, p2::AbstractArray) = p1 + PolyMatrix(p2, V)
+{T,M,V,N}(p1::AbstractArray, p2::PolyMatrix{T,M,Var{V},N}) = PolyMatrix(p1, V) + p2

function -{T1,M1,V,N}(p::PolyMatrix{T1,M1,Var{V},N})
  # figure out return type
  c     = coeffs(p)
  k1,v1 = first(c)
  vr    = -v1
  M     = typeof(vr)
  r     = PolyMatrix( SortedDict{Int,M,ForwardOrdering}(), size(p), V)
  for (k,v) in c
    insert!(coeffs(r), k, -coeffs(p)[k])
  end
  return r
end

-{T1,M1,V,N,T2,M2}(p1::PolyMatrix{T1,M1,Var{V},N}, p2::PolyMatrix{T2,M2,Var{V},N}) = +(p1,-p2)
function -{T1,M1,V1,V2,N,T2,M2}(p1::PolyMatrix{T1,M1,Var{V1},N}, p2::PolyMatrix{T2,M2,Var{V2},N})
  warn("p1≈p2: `p1` ($T1,$V1) and `p2` ($T2,$V2) have different variables")
  throw(DomainError())
end

-{T1,M1,V,N}(p1::PolyMatrix{T1,M1,Var{V},N}, p2::AbstractArray) = p1 - PolyMatrix(p2, V)
-{T1,M1,V,N}(p1::AbstractArray, p2::PolyMatrix{T1,M1,Var{V},N}) = PolyMatrix(p1, V) - p2

function *{T1,M1,V,N,T2,M2}(p1::PolyMatrix{T1,M1,Var{V},N}, p2::PolyMatrix{T2,M2,Var{V},N})
  size(p1,2) == size(p2,1) || error("incompatible sizes")
  PMcheck(p1,p2)

  # figure out return type
  c1     = coeffs(p1)
  c2     = coeffs(p2)
  k1, v1 = first(c1)
  k2, v2 = first(c2)
  vr     = v1*v2
  M      = typeof(vr)
  r      = PolyMatrix( SortedDict{Int,M,ForwardOrdering}(), size(vr), V)

  # find all new powers k1+k2 and corresponding k1, k2
  klist = Dict{Int,Vector{Tuple}}()
  for k1 in keys(c1)
    for k2 in keys(c2)
      klist[k1+k2] = push!(get(klist,k1+k2, Vector{Tuple}()), tuple(k1,k2))
    end
  end

  # do the calculations
  for k in keys(klist)
    vk = spzeros(size(r)...)
    for v in klist[k]
      vk += c1[v[1]]*c2[v[2]]
    end
    r.coeffs[k] = vk
  end
  return r
end

function *{T1,M1,V1,V2,N,T2,M2}(p1::PolyMatrix{T1,M1,Var{V1},N}, p2::PolyMatrix{T2,M2,Var{V2},N})
  warn("p1≈p2: `p1` ($T1,$V1) and `p2` ($T2,$V2) have different variables")
  throw(DomainError())
end

*{T1,M1,V,N}(p1::PolyMatrix{T1,M1,Var{V},N}, p2::AbstractArray) = p1 * PolyMatrix(p2, V)
*{T1,M1,V,N}(p1::AbstractArray, p2::PolyMatrix{T1,M1,Var{V},N}) = PolyMatrix(p1, V) * p2
