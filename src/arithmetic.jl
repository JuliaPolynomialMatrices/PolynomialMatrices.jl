function PMcheck(p1::PolyMatrix, p2::PolyMatrix)
  @assert 1==1
end

function +{T1,M1,O,N,T2,M2}(p1::PolyMatrix{T1,M1,O,N}, p2::PolyMatrix{T2,M2,O,N})
  size(p1) == size(p2) || error("incompatible sizes")
  PMcheck(p1,p2)

  # figure out return type
  c1     = coeffs(p1)
  c2     = coeffs(p2)
  k1, v1 = first(c1)
  k2, v2 = first(c2)
  vr     = v1+v2
  M      = typeof(vr)
  r      = PolyMatrix( SortedDict{Int,M,O}(), size(p1), p1.var)

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
+(p1::PolyMatrix, p2::AbstractArray) = p1 + PolyMatrix(p2,p1.var)
+(p1::AbstractArray, p2::PolyMatrix) = PolyMatrix(p1,p2.var) + p2

function -{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  # figure out return type
  r  = PolyMatrix( SortedDict{Int,M,O}(), size(p), p.var)
  for (k,v) in coeffs(p)
    insert!(coeffs(r), k, -coeffs(p)[k])
  end
  return r
end
-{T1,M1,O,N,T2,M2}(p1::PolyMatrix{T1,M1,O,N}, p2::PolyMatrix{T2,M2,O,N}) = +(p1,-p2)
-(p1::PolyMatrix, p2::AbstractArray) = p1 - PolyMatrix(p2,p1.var)
-(p1::AbstractArray, p2::PolyMatrix) = PolyMatrix(p1,p2.var) - p2

function *{T1,M1,O,N,T2,M2}(p1::PolyMatrix{T1,M1,O,N}, p2::PolyMatrix{T2,M2,O,N})
  size(p1,2) == size(p2,1) || error("incompatible sizes")
  p1.var == p2.var ||
    error("multiplication of polynomial matrices with different variables not supported")
  PMcheck(p1,p2)

  # figure out return type
  c1     = coeffs(p1)
  c2     = coeffs(p2)
  k1, v1 = first(c1)
  k2, v2 = first(c2)
  vr     = v1*v2
  M      = typeof(vr)
  r      = PolyMatrix( SortedDict{Int,M,O}(), size(vr), p1.var)

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
*(p1::PolyMatrix, p2::AbstractArray) = p1 * PolyMatrix(p2,p1.var)
*(p1::AbstractArray, p2::PolyMatrix) = PolyMatrix(p1,p2.var) * p2
