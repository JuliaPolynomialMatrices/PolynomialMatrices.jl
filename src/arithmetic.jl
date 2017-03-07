## Basic operations between polynomial matrices
function +{T1,M1,V,N,T2,M2}(p1::PolyMatrix{T1,M1,Val{V},N}, p2::PolyMatrix{T2,M2,Val{V},N})
  if size(p1) ≠ size(p2)
    warn("+(p1,p2): size(p1) ≠ size(p2)")
    throw(DomainError())
  end
  _add(p1,p2)
end

function _add{T1,M1,V,N,T2}(p1::PolyMatrix{T1,M1,Val{V},N}, p2::T2)
  c1 = coeffs(p1)
  c2 = coeffs(p2)
  v1 = first(c1)[end] # for polynomials first(c1) returns index of first element
  v2 = first(c2)[end] # for polynomialmatrices it returns key value pair
  vr = v1+v2
  M  = typeof(vr)
  r  = PolyMatrix(SortedDict{Int,M,ForwardOrdering}(), size(vr), Val{V})

  cr  = coeffs(r)
  sᵢ  = intersect(_keys(c1), _keys(c2))
  s₁  = setdiff(_keys(c1), sᵢ)
  s₂  = setdiff(_keys(c2), sᵢ)
  for k in s₁
    insert!(cr, k, c1[k])
  end
  for k in s₂
    v = T2 <:Poly ? p2[k] : c2[k]
    insert!(cr, k, v)
  end
  for k in sᵢ
    v = T2 <:Poly ? c1[k]+p2[k] : c1[k]+c2[k]
    insert!(cr, k, v)
  end
  return r
end

function +{T1,M1,V1,V2,N,T2,M2}(p1::PolyMatrix{T1,M1,Val{V1},N}, p2::PolyMatrix{T2,M2,Val{V2},N})
  warn("p1+p2: `p1` ($T1,$V1) and `p2` ($T2,$V2) have different variables")
  throw(DomainError())
end

function -{T1,M1,V,N}(p::PolyMatrix{T1,M1,Val{V},N})
  r = PolyMatrix(SortedDict{Int,M1,ForwardOrdering}(), size(p), Val{V})
  for (k,v) in coeffs(p)
    insert!(coeffs(r), k, -coeffs(p)[k])
  end
  return r
end

-{T1,M1,V,N,T2,M2}(p1::PolyMatrix{T1,M1,Val{V},N}, p2::PolyMatrix{T2,M2,Val{V},N}) = +(p1,-p2)

function -{T1,M1,V1,V2,N,T2,M2}(p1::PolyMatrix{T1,M1,Val{V1},N}, p2::PolyMatrix{T2,M2,Val{V2},N})
  warn("p1-p2: `p1` ($T1,$V1) and `p2` ($T2,$V2) have different variables")
  throw(DomainError())
end

# heuristic used below was found by benchmarking
# (number of matrix multiplications of _mul is length(c1)*length(c2)
# (number of matrix multiplications of _mulfft is degree(p1)+degree(p2)
function *{T1,M1,V,T2,M2}(p1::PolyMatrix{T1,M1,Val{V},2}, p2::PolyMatrix{T2,M2,Val{V},2})
  if size(p1,2) ≠ size(p2,1)
    warn("*(p1,p2): size(p1,2) ≠ size(p2,1)")
    throw(DomainError())
  end
  _mul(p1,p2)
end

function *{T1,M1,V,T2,M2}(p1::PolyMatrix{T1,M1,Val{V},2}, p2::PolyMatrix{T2,M2,Val{V},1})
  if size(p1,2) ≠ size(p2,1)
    warn("*(p1,p2): size(p1,2) ≠ size(p2,1)")
    throw(DomainError())
  end
  _mul(p1,p2)
end

function _mul{T1,T2}(p1::T1, p2::T2)
  c1     = coeffs(p1)
  c2     = coeffs(p2)
  if length(c1)*length(c2) > 5*(degree(p1)+degree(p2))
    return _mulfft(p1,p2)
  else
    return _mulconv(p1,p2)
  end
end

function _mulconv{T1,M1,V,N,T2}(p1::PolyMatrix{T1,M1,Val{V},N}, p2::T2)
  # figure out return type
  c1 = coeffs(p1)
  c2 = coeffs(p2)
  v1 = first(c1)[end] # for polynomials first(c1) returns index of first element
  v2 = first(c2)[end] # for polynomialmatrices it returns key value pair
  vr = v1*v2

  M  = typeof(vr)
  r  = PolyMatrix(SortedDict{Int,M,ForwardOrdering}(), size(vr), Val{V})

  # find all new powers k1+k2 and corresponding k1, k2
  klist = Dict{Int,Vector{Tuple{Int,Int}}}()
  for k1 in _keys(c1)
    for k2 in _keys(c2)
      klist[k1+k2] = push!(get(klist,k1+k2, Vector{Tuple{Int,Int}}()), tuple(k1,k2))
    end
  end

  # do the calculations
  for k in keys(klist)
    vk = zeros(vr)
    for v in klist[k]
      vk += T2 <:Poly ? c1[v[1]]*p2[v[2]] : c1[v[1]]*c2[v[2]]
    end
    r.coeffs[k] = vk
  end
  return r
end

_keys{T}(c::T) = keys(c)
_keys{T<:AbstractArray}(c::T) = eachindex(c)-1

function _mulfft{T1,M1,V,N,T2}(p1::PolyMatrix{T1,M1,Val{V},N}, p2::T2)
  T = promote_type(T1, eltype(eltype(T2)))

  c1 = coeffs(p1)
  c2 = coeffs(p2)
  v1 = first(c1)[end] # for polynomials first(c1) returns index of first element
  v2 = first(c2)[end] # for polynomialmatrices it returns key value pair
  vr = v1*v2

  n,m = size(vr)
  dn  = degree(p1)+degree(p2)+1
  # copy all elements into three-dimensional matrices
  A1  = _fftmatrix(p1, T, dn)
  A2  = _fftmatrix(p2, T, dn)
  # take fft and evaluate multiplication at each interpolation point
  B1  = fft(A1, ndims(A1))
  B2  = fft(A2, ndims(A2))

  a = zeros(eltype(B1),n,m,dn)
  @inbounds @simd for k in 1:dn
    a[:,:,k] += T2 <:Poly ? B1[:,:,k]*B2[k] : B1[:,:,k]*B2[:,:,k]
  end
  # interpolate using fft
  ar = _truncate(T,ifft(a,3))
  return PolyMatrix(ar, Val{V})
end

function _fftmatrix{T1,M1,V1,N,T}(p::PolyMatrix{T1,M1,Val{V1},N}, ::Type{T}, dn::Integer)
  A = zeros(T, size(p)..., dn)
  for (k,v) in coeffs(p)
    A[:,:,k+1] = v
  end
  A
end

function _fftmatrix{T1,T}(p::Poly{T1}, ::Type{T}, dn::Integer)
  A = zeros(T,dn)
  c = coeffs(p)
  for i in eachindex(c)
    A[i] = c[i]
  end
  A
end

function *{T1,M1,V1,V2,N,T2,M2}(p1::PolyMatrix{T1,M1,Val{V1},N}, p2::PolyMatrix{T2,M2,Val{V2},N})
  warn("p1*p2: `p1` ($T1,$V1) and `p2` ($T2,$V2) have different variables")
  throw(DomainError())
end

## Basic operations between polynomial matrices and AbstractArrays
+{T,M1,V,N,M2<:AbstractArray}(p1::PolyMatrix{T,M1,Val{V},N}, p2::M2) = p1 + PolyMatrix(p2, Val{V})
+{T,M1,V,N,M2<:AbstractArray}(p1::M2, p2::PolyMatrix{T,M1,Val{V},N}) = PolyMatrix(p1, Val{V}) + p2

-{T,M1,V,N,M2<:AbstractArray}(p1::PolyMatrix{T,M1,Val{V},N}, p2::M2) = p1 - PolyMatrix(p2, Val{V})
-{T,M1,V,N,M2<:AbstractArray}(p1::M2, p2::PolyMatrix{T,M1,Val{V},N}) = PolyMatrix(p1, Val{V}) - p2

*{T,M1,V,N,S}(p1::PolyMatrix{T,M1,Val{V},N}, p2::AbstractArray{S,2}) = p1*PolyMatrix(p2, Val{V})
*{T,M1,V,N,S}(p2::AbstractArray{S,2}, p1::PolyMatrix{T,M1,Val{V},N}) = PolyMatrix(p2, Val{V})*p1

*{T,M1,V,N,S<:Number}(p1::PolyMatrix{T,M1,Val{V},N}, p2::AbstractArray{S,1}) = p1*PolyMatrix(p2, size(p2), Val{V})

## Basic operations between polynomial matrices and polynomials
+{T,M1,V,N,M2<:Poly}(p1::PolyMatrix{T,M1,Val{V},N}, p2::M2) = _add(p1, p2)
+{T,M1,V,N,M2<:Poly}(p1::M2, p2::PolyMatrix{T,M1,Val{V},N}) = _add(p2, p1)

-{T,M1,V,N,M2<:Poly}(p1::PolyMatrix{T,M1,Val{V},N}, p2::M2) = _add(p1, -p2)
-{T,M1,V,N,M2<:Poly}(p1::M2, p2::PolyMatrix{T,M1,Val{V},N}) = _add(-p2, p1)

*{T,M1,V,N,M2<:Poly}(p1::PolyMatrix{T,M1,Val{V},N}, p2::M2) = _mul(p1, p2)
*{T,M1,V,N,M2<:Poly}(p1::M2, p2::PolyMatrix{T,M1,Val{V},N}) = _mul(p2, p1)

## Basic operations between polynomial matrices and Numbers
function _add{T1,M1,V,N,T2<:Number}(p1::PolyMatrix{T1,M1,Val{V},N}, v2::T2)
  T   = promote_type(T1, T2)
  c1  = coeffs(p1)
  v1  = first(c1)[end] # for polynomials first(c1) returns index of first element
  M   = typeof(similar(v1, T))
  r   = PolyMatrix( SortedDict{Int,M,ForwardOrdering}(), size(p1), Val{V})
  cr  = coeffs(r)

  for (k1,v1) in coeffs(p1)
    insert!(cr, k1, v1)
  end
  cr[0] += v2
  return r
end

+{T,M1,V,N,M2<:Number}(p1::PolyMatrix{T,M1,Val{V},N}, p2::M2) = _add(p1, p2)
+{T,M1,V,N,M2<:Number}(p1::M2, p2::PolyMatrix{T,M1,Val{V},N}) = _add(p2, p1)

-{T,M1,V,N,M2<:Number}(p1::PolyMatrix{T,M1,Val{V},N}, p2::M2) = _add(p1, -p2)
-{T,M1,V,N,M2<:Number}(p1::M2, p2::PolyMatrix{T,M1,Val{V},N}) = _add(-p2, p1)

function _mul{T1,M1,V,N,T2<:Number}(p1::PolyMatrix{T1,M1,Val{V},N}, v2::T2)
  T   = promote_type(T1, T2)
  c1  = coeffs(p1)
  v1  = first(c1)[end] # for polynomials first(c1) returns index of first element
  M   = typeof(similar(v1, T))
  r   = PolyMatrix( SortedDict{Int,M,ForwardOrdering}(), size(p1), Val{V})
  cr  = coeffs(r)

  for (k1,v1) in coeffs(p1)
    insert!(cr, k1, v1*v2)
  end
  return r
end

*{T,M1,V,N,M2<:Number}(p1::PolyMatrix{T,M1,Val{V},N}, p2::M2) = _mul(p1, p2)
*{T,M1,V,N,M2<:Number}(p1::M2, p2::PolyMatrix{T,M1,Val{V},N}) = _mul(p2, p1)

/{T,M1,V,N,M2<:Number}(p1::PolyMatrix{T,M1,Val{V},N}, p2::M2) = _mul(p1, 1/p2)

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
  ar = _truncate(T, ifft(a))
  return Poly(ar, V)
end

function _truncate{T<:Real,T2}(::Type{T}, a::AbstractArray{T2})
  r = similar(a, T)
  r = real(a)
end

function _truncate{T<:Integer,T2}(::Type{T}, a::AbstractArray{T2})
  r = similar(a, T)
  for i in eachindex(a)
    r[i] = convert(T, round(real(a[i])))
  end
  r
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
