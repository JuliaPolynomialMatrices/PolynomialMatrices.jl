function PMcheck(p1::PolyMatrix, p2::PolyMatrix)
  @assert 1==1
end

function +{T1,M1,V,N,T2,M2}(p1::PolyMatrix{T1,M1,Val{V},N}, p2::PolyMatrix{T2,M2,Val{V},N})
  size(p1) == size(p2) || error("incompatible sizes")
  PMcheck(p1,p2)

  # figure out return type
  c1     = coeffs(p1)
  c2     = coeffs(p2)
  k1, v1 = first(c1)
  k2, v2 = first(c2)
  vr     = v1+v2
  M      = typeof(vr)
  r      = PolyMatrix( SortedDict{Int,M,ForwardOrdering}(), size(p1), Val{V})

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

function +{T1,M1,V1,V2,N,T2,M2}(p1::PolyMatrix{T1,M1,Val{V1},N}, p2::PolyMatrix{T2,M2,Val{V2},N})
  warn("p1≈p2: `p1` ($T1,$V1) and `p2` ($T2,$V2) have different variables")
  throw(DomainError())
end

+{T,M,V,N}(p1::PolyMatrix{T,M,Val{V},N}, p2::AbstractArray) = p1 + PolyMatrix(p2, Val{V})
+{T,M,V,N}(p1::AbstractArray, p2::PolyMatrix{T,M,Val{V},N}) = PolyMatrix(p1, Val{V}) + p2

function -{T1,M1,V,N}(p::PolyMatrix{T1,M1,Val{V},N})
  # figure out return type
  c     = coeffs(p)
  k1,v1 = first(c)
  vr    = -v1
  M     = typeof(vr)
  r     = PolyMatrix( SortedDict{Int,M,ForwardOrdering}(), size(p), Val{V})
  for (k,v) in c
    insert!(coeffs(r), k, -coeffs(p)[k])
  end
  return r
end

-{T1,M1,V,N,T2,M2}(p1::PolyMatrix{T1,M1,Val{V},N}, p2::PolyMatrix{T2,M2,Val{V},N}) = +(p1,-p2)
function -{T1,M1,V1,V2,N,T2,M2}(p1::PolyMatrix{T1,M1,Val{V1},N}, p2::PolyMatrix{T2,M2,Val{V2},N})
  warn("p1≈p2: `p1` ($T1,$V1) and `p2` ($T2,$V2) have different variables")
  throw(DomainError())
end

-{T1,M1,V,N}(p1::PolyMatrix{T1,M1,Val{V},N}, p2::AbstractArray) = p1 - PolyMatrix(p2, Val{V})
-{T1,M1,V,N}(p1::AbstractArray, p2::PolyMatrix{T1,M1,Val{V},N}) = PolyMatrix(p1, Val{V}) - p2

# heuristic used below was found by benchmarking
# (number of matrix multiplications of _mul is length(c1)*length(c2)
# (number of matrix multiplications of _mulfft is degree(p1)+degree(p2)
function *{T1,M1,V,N,T2,M2}(p1::PolyMatrix{T1,M1,Val{V},N}, p2::PolyMatrix{T2,M2,Val{V},N})
  size(p1,2) == size(p2,1) || error("incompatible sizes")
  PMcheck(p1,p2)

  # figure out return type
  c1     = coeffs(p1)
  c2     = coeffs(p2)
  if length(c1)*length(c2) > 5*(degree(p1)+degree(p2))
    return _mulfft(p1,p2)
  else
    return _mul(p1,p2)
  end
end

function _mul{T1,M1,V,N,T2,M2}(p1::PolyMatrix{T1,M1,Val{V},N}, p2::PolyMatrix{T2,M2,Val{V},N})
  # figure out return type
  c1     = coeffs(p1)
  c2     = coeffs(p2)
  k1, v1 = first(c1)
  k2, v2 = first(c2)
  vr     = v1*v2

  M      = typeof(vr)
  r      = PolyMatrix( SortedDict{Int,M,ForwardOrdering}(), size(vr), Val{V})

  # find all new powers k1+k2 and corresponding k1, k2
  klist = Dict{Int,Vector{Tuple{Int,Int}}}()
  for k1 in keys(c1)
    for k2 in keys(c2)
      klist[k1+k2] = push!(get(klist,k1+k2, Vector{Tuple{Int,Int}}()), tuple(k1,k2))
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

function _mulfft{T1,M1,V,N,T2,M2}(p1::PolyMatrix{T1,M1,Val{V},N}, p2::PolyMatrix{T2,M2,Val{V},N})
  T = promote_type(T1,T2)

  n  = size(p1,1)
  m  = size(p2,2)
  dn = degree(p1)+degree(p2)+1
  # copy all elements into three-dimensional matrices
  A1 = zeros(T, size(p1)..., dn)
  for (k,v) in coeffs(p1)
    A1[:,:,k+1] = v
  end
  A2 = zeros(T, size(p2)..., dn)
  for (k,v) in coeffs(p2)
    A2[:,:,k+1] = v
  end
  # take fft and evaluate multiplication at each interpolation point
  B1 = fft(A1,3)
  B2 = fft(A2,3)

  a = zeros(eltype(B1),n,m,dn)
  @inbounds @simd for k in 1:dn
    a[:,:,k] += B1[:,:,k]*B2[:,:,k]
  end
  # interpolate using fft
  ar = _truncate(T,ifft(a,3))
  return PolyMatrix(ar, Val{V})
end

function *{T1,M1,V1,V2,N,T2,M2}(p1::PolyMatrix{T1,M1,Val{V1},N}, p2::PolyMatrix{T2,M2,Val{V2},N})
  warn("p1≈p2: `p1` ($T1,$V1) and `p2` ($T2,$V2) have different variables")
  throw(DomainError())
end

*{T1,M1,V,N}(p1::PolyMatrix{T1,M1,Val{V},N}, p2::AbstractArray) = p1 * PolyMatrix(p2, Val{V})
*{T1,M1,V,N}(p1::AbstractArray, p2::PolyMatrix{T1,M1,Val{V},N}) = PolyMatrix(p1, Val{V}) * p2

# determinant
function det{T,M,V,N}(p::PolyMatrix{T,M,Val{V},N})
  size(p,1) == size(p,2) || throw(DimensionMismatch("det: PolyMatrix must be square"))
  n  = size(p,1)
  dn = (degree(p))*n+1
  # copy all elements into three-dimensional matrix
  A = zeros(n,n,dn)
  for (k,v) in coeffs(p)
    A[:,:,k+1] = v
  end
  # take fft and evaluate determinant at each interpolation point
  B = fft(A,3)
  a = [det(B[:,:,k]) for k = 1:dn]
  # interpolate using fft
  ar = _truncate(T,ifft(a))
  return Poly(ar, size(ar), Val{V})
end

function _truncate{T<:Real,T2}(::Type{T}, a::AbstractArray{T2})
  real(a)
end

function _truncate{T<:Integer,T2}(::Type{T}, a::AbstractArray{T2})
  round(real(a))
end

# inversion
# return determinant polynomial and adjugate polynomial matrix
function inv{T,M,V,N}(p::PolyMatrix{T,M,Val{V},N})
  size(p,1) == size(p,2) || throw(DimensionMismatch("det: PolyMatrix must be square"))
  n  = size(p,1)
  dn = degree(p)*n+1
  # copy all elements into three-dimensional matrix
  A = zeros(T,n,n,dn)
  for (k,v) in coeffs(p)
    A[:,:,k+1] = v
  end
  # take fft and evaluate determinant at each interpolation point
  B = fft(A,3)
  a = [det(B[:,:,k]) for k = 1:dn]
  ar = _truncate(T,ifft(a))
  rdet = Poly(ar, V)

  v2  = zeros(B)
  for k in 1:dn
    for col in 1:n
      for row in 1:n
        row_idx = _detrange(row,n)
        col_idx = _detrange(col,n)
        v2[row,col,k] = (-1)^(mod(row+col,2))*det(view(view(B, :, row_idx, k), col_idx, :))
      end
    end
  end
  r = _truncate(T, ifft(v2,3))

  return rdet, PolyMatrix(r, size(r), Val{V})
end

function _detrange(i,n)
  a = 1:i-1
  b = i+1:n
  r = vcat(collect(a), collect(b))
  return ifelse(length(r) > 1, r, r[1]:r[1])
end
