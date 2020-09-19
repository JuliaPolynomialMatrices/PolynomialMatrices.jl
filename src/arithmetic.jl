zero## Basic operations between polynomial matrices
function +(p1::PolyMatrix{T1,M1,Val{W},N}, p2::PolyMatrix{T2,M2,Val{W},N}) where {T1,M1,W,N,T2,M2}
  if size(p1) ≠ size(p2)
    throw(DimensionMismatch("the two polynomial matrices are of incompatible sizes"))
  end
  _add(p1,p2)
end

function _add(p1::PolyMatrix{T1,M1,Val{W},N},
  p2::PolyMatrix{T2,M2,Val{W},N}) where {T1,M1,W,N,T2,M2}
  c1    = coeffs(p1)
  c2    = coeffs(p2)
  _,v1  = first(c1)
  _,v2  = first(c2)
  vr    = v1+v2
  M     = typeof(vr)

  cr  = SortedDict(0=>zero(vr))
  sᵢ  = intersect(keys(c1), keys(c2))
  s₁  = setdiff(keys(c1), sᵢ)
  s₂  = setdiff(keys(c2), sᵢ)
  for k in s₁
    insert!(cr, k, copy(c1[k]))
  end
  for k in s₂
    insert!(cr, k, copy(c2[k]))
  end
  for k in sᵢ
    insert!(cr, k, c1[k]+c2[k])
  end
  return PolyMatrix(cr, size(vr), Val{W})
end

function _add(p1::PolyMatrix{T1,M1,Val{W},N}, p2::T2) where {T1,M1,W,N,T2<:Polynomial}
  c1    = coeffs(p1)
  c2    = coeffs(p2)
  _,v1  = first(c1) # for polynomials first(c1) returns index of first element
  v2    = first(c2)   # for polynomialmatrices it returns key value pair
  vr    = v1 .+ v2
  M     = typeof(vr)

  cr  = SortedDict(Dict{Int,M}())
  sᵢ  = intersect(_keys(c1), _keys(c2))
  s₁  = setdiff(_keys(c1), sᵢ)
  s₂  = setdiff(_keys(c2), sᵢ)
  for k in s₁
    insert!(cr, k, copy(c1[k]))
  end
  for k in s₂
    insert!(cr, k, copy(p2[k]))
  end
  for k in sᵢ
    insert!(cr, k, c1[k] .+ p2[k])
  end
  return PolyMatrix(cr, size(vr), Val{W})
end

function +(p1::PolyMatrix{T1,M1,Val{W1},N}, p2::PolyMatrix{T2,M2,Val{W2},N}) where {T1,M1,W1,W2,N,T2,M2}
  throw(ArgumentError("the two polynomial matrices have different variables"))
end


function -(p::PolyMatrix{T1,M1,Val{W},N}) where {T1,M1,W,N}
  cr = SortedDict(Dict{Int,M1}())
  for (k,v) in coeffs(p)
    insert!(cr, k, -v)
  end
  return PolyMatrix(cr, size(p), Val{W})
end

-(p1::PolyMatrix{T1,M1,Val{W},N}, p2::PolyMatrix{T2,M2,Val{W},N}) where {T1,M1,W,N,T2,M2} = +(p1,-p2)

function -(p1::PolyMatrix{T1,M1,Val{W1},N}, p2::PolyMatrix{T2,M2,Val{W2},N}) where {T1,M1,W1,W2,N,T2,M2}
  throw(ArgumentError("the two polynomial matrices have different variables"))
end

# heuristic used below was found by benchmarking
function *(p1::PolyMatrix{T1,M1,Val{W},2}, p2::PolyMatrix{T2,M2,Val{W},2}) where {T1,M1,W,T2,M2}
  if size(p1,2) ≠ size(p2,1)
    throw(DimensionMismatch("the two polynomial matrices are of incompatible sizes"))
  end
  _mul(p1,p2)
end

function *(p1::PolyMatrix{T1,M1,Val{W},2}, p2::PolyMatrix{T2,M2,Val{W},1}) where {T1,M1,W,T2,M2}
  if size(p1,2) ≠ size(p2,1)
    throw(DimensionMismatch("the two polynomial matrices are of incompatible sizes"))
  end
  _mul(p1,p2)
end

function _mul(p1::T1, p2::T2) where {T1,T2}
  degree_cutoff = 45
  if degree(p1)+degree(p2) > degree_cutoff
    return _mulfft(p1,p2)
  else
    return _mulconv(p1,p2)
  end
end

function _mulconv(p1::PolyMatrix{T1,M1,Val{W},2},
  p2::PolyMatrix{T2,M2,Val{W},N}) where {T1,M1,W,T2,M2,N}
  # figure out return type
  c1    = coeffs(p1)
  c2    = coeffs(p2)
  _,v1  = first(c1) # for polynomials first(c1) returns index of first element
  _,v2  = first(c2) # for polynomialmatrices it returns key value pair
  vr    = v1*v2
  M     = typeof(vr)
  cr    = SortedDict(Dict{Int,M}())

  # find all new powers k1+k2 and corresponding k1, k2
  klist = Dict{Int,Vector{Tuple{Int,Int}}}()
  for k1 in keys(c1), k2 in keys(c2)
    klist[k1+k2] = push!(get(klist,k1+k2, Vector{Tuple{Int,Int}}()), tuple(k1,k2))
  end

  # do the calculations
  for k in keys(klist)
    vk = zero(vr)
    for v in klist[k]
      vk += c1[v[1]]*c2[v[2]]
    end
    insert!(cr, k, vk)
  end
  return PolyMatrix(cr, size(vr), Val{W})
end

_mulconv(p1::T2, p2::PolyMatrix{T1,M1,Val{W},N}) where {T1,M1,W,N,T2<:Polynomial} = _mulconv(p2,p1)

function _mulconv(p1::PolyMatrix{T1,M1,Val{W},N}, p2::T2) where {T1,M1,W,N,T2<:Polynomial}
  # figure out return type
  c1    = coeffs(p1)
  c2    = coeffs(p2)
  _,v1  = first(c1) # for polynomials first(c1) returns index of first element
  v2    = first(c2) # for polynomialmatrices it returns key value pair
  vr    = v1*v2
  M     = typeof(vr)
  cr    = SortedDict(Dict{Int,M}())

  # find all new powers k1+k2 and corresponding k1, k2
  klist = Dict{Int,Vector{Tuple{Int,Int}}}()
  for k1 in _keys(c1)
    for k2 in _keys(c2)
      klist[k1+k2] = push!(get(klist,k1+k2, Vector{Tuple{Int,Int}}()), tuple(k1,k2))
    end
  end

  # do the calculations
  for k in keys(klist)
    vk = zero(vr)
    for v in klist[k]
      vk += c1[v[1]]*p2[v[2]]
    end
    insert!(cr, k, vk)
  end
  return PolyMatrix(cr, size(vr), Val{W})
end

_keys(c::T) where {T} = keys(c)
_keys(c::T) where {T<:AbstractArray} = eachindex(c) .- 1

function _mulfft(p1::PolyMatrix{T1,M1,Val{W},2},
  p2::PolyMatrix{T2,M2,Val{W},N}) where {T1,M1,W,T2,M2,N}
  T     = promote_type(T1, T2)
  c1    = coeffs(p1)
  c2    = coeffs(p2)
  _,v1  = first(c1) # for polynomials first(c1) returns index of first element
  _,v2  = first(c2) # for polynomialmatrices it returns key value pair
  vr    = v1*v2

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
    a[:,:,k] += B1[:,:,k]*B2[:,:,k]
  end
  # interpolate using fft
  ar = _truncate(T,ifft(a,3))
  return PolyMatrix(ar, Val{W})
end

_mulfft(p1::Polynomial{T2}, p2::PolyMatrix{T1,M1,Val{W},N}) where {T1,M1,W,N,T2} = _mulfft(p2,p1)

function _mulfft(p1::PolyMatrix{T1,M1,Val{W},N}, p2::Polynomial{T2}) where {T1,M1,W,N,T2}
  T     = promote_type(T1, T2)
  c1    = coeffs(p1)
  c2    = coeffs(p2)
  _,v1  = first(c1) # for polynomials first(c1) returns index of first element
  v2    = first(c2) # for polynomialmatrices it returns key value pair
  vr    = v1*v2

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
    a[:,:,k] += B1[:,:,k]*B2[k]
  end
  # interpolate using fft
  ar = _truncate(T,ifft(a,3))
  return PolyMatrix(ar, Val{W})
end

function _fftmatrix(p::PolyMatrix{T1,M1,Val{W1},N}, ::Type{T}, dn::Integer) where {T1,M1,W1,N,T}
  A = zeros(T, size(p)..., dn)
  for (k,v) in coeffs(p)
    A[:,:,k+1] = v
  end
  A
end

function _fftmatrix(p::Polynomial{T1}, ::Type{T}, dn::Integer) where {T1,T}
  A = zeros(T,dn)
  c = coeffs(p)
  for i in eachindex(c)
    A[i] = c[i]
  end
  A
end

function *(p1::PolyMatrix{T1,M1,Val{W1},N}, p2::PolyMatrix{T2,M2,Val{W2},N}) where {T1,M1,W1,W2,N,T2,M2}
  throw(ArgumentError("the two polynomial matrices have different variables"))
end

## Basic operations between polynomial matrices and AbstractArrays
+(p1::PolyMatrix{T,M1,Val{W},N}, p2::M2) where {T,M1,W,N,M2<:AbstractArray} = p1 + PolyMatrix(p2, Val{W})
+(p1::M2, p2::PolyMatrix{T,M1,Val{W},N}) where {T,M1,W,N,M2<:AbstractArray} = PolyMatrix(p1, Val{W}) + p2

-(p1::PolyMatrix{T,M1,Val{W},N}, p2::M2) where {T,M1,W,N,M2<:AbstractArray} = p1 - PolyMatrix(p2, Val{W})
-(p1::M2, p2::PolyMatrix{T,M1,Val{W},N}) where {T,M1,W,N,M2<:AbstractArray} = PolyMatrix(p1, Val{W}) - p2

*(p1::PolyMatrix{T,M1,Val{W},N}, p2::AbstractArray{S,2}) where {T,M1,W,N,S} = p1*PolyMatrix(p2, Val{W})
*(p2::AbstractArray{S,2}, p1::PolyMatrix{T,M1,Val{W},N}) where {T,M1,W,N,S} = PolyMatrix(p2, Val{W})*p1

*(p1::PolyMatrix{T,M1,Val{W},N}, p2::AbstractArray{S,1}) where {T,M1,W,N,S<:Number} = p1*PolyMatrix(p2, size(p2), Val{W})

## Basic operations between polynomial matrices and polynomials
+(p1::PolyMatrix{T,M1,Val{W},N}, p2::M2) where {T,M1,W,N,M2<:Polynomial} = _add(p1, p2)
+(p1::M2, p2::PolyMatrix{T,M1,Val{W},N}) where {T,M1,W,N,M2<:Polynomial} = _add(p2, p1)

-(p1::PolyMatrix{T,M1,Val{W},N}, p2::M2) where {T,M1,W,N,M2<:Polynomial} = _add(p1, -p2)
-(p1::M2, p2::PolyMatrix{T,M1,Val{W},N}) where {T,M1,W,N,M2<:Polynomial} = _add(-p2, p1)

*(p1::PolyMatrix{T,M1,Val{W},N}, p2::M2) where {T,M1,W,N,M2<:Polynomial} = _mul(p1, p2)
*(p1::M2, p2::PolyMatrix{T,M1,Val{W},N}) where {T,M1,W,N,M2<:Polynomial} = _mul(p2, p1)

## Basic operations between polynomial matrices and Numbers
function _add(p1::PolyMatrix{T1,M1,Val{W},N}, v2::T2) where {T1,M1,W,N,T2<:Number}
  T     = promote_type(T1, T2)
  c1    = coeffs(p1)
  _,v1  = first(c1)    # for polynomials first(c1) returns index of first element
  M     = typeof(similar(v1, T))
  cr    = SortedDict(Dict{Int,M}())

  for (k1,v1) in coeffs(p1)
    insert!(cr, k1, copy(v1))
  end
  broadcast!(+,cr[0],cr[0],v2)
  return PolyMatrix(cr, size(p1), Val{W})
end

+(p1::PolyMatrix{T,M1,Val{W},N}, p2::M2) where {T,M1,W,N,M2<:Number} = _add(p1, p2)
+(p1::M2, p2::PolyMatrix{T,M1,Val{W},N}) where {T,M1,W,N,M2<:Number} = _add(p2, p1)

-(p1::PolyMatrix{T,M1,Val{W},N}, p2::M2) where {T,M1,W,N,M2<:Number} = _add(p1, -p2)
-(p1::M2, p2::PolyMatrix{T,M1,Val{W},N}) where {T,M1,W,N,M2<:Number} = _add(-p2, p1)

function _mul(p1::PolyMatrix{T1,M1,Val{W},N}, v2::T2) where {T1,M1,W,N,T2<:Number}
  T     = promote_type(T1, T2)
  c1    = coeffs(p1)
  _,v1  = first(c1) # for polynomials first(c1) returns index of first element
  M     = typeof(similar(v1, T))
  cr    = SortedDict(Dict{Int,M}())

  for (k1,v1) in coeffs(p1)
    insert!(cr, k1, v1*v2)
  end
  return PolyMatrix(cr, size(p1), Val{W})
end

*(p1::PolyMatrix{T,M1,Val{W},N}, p2::M2) where {T,M1,W,N,M2<:Number} = _mul(p1, p2)
*(p1::M2, p2::PolyMatrix{T,M1,Val{W},N}) where {T,M1,W,N,M2<:Number} = _mul(p2, p1)

/(p1::PolyMatrix{T,M1,Val{W},N}, p2::M2) where {T,M1,W,N,M2<:Number} = _mul(p1, 1/p2)

# multiplication with abstractArray
# transpose
# Base.LinAlg.At_mul_B(p1::PolyMatrix{T1,M1,Val{W},N1},
#   p2::AbstractArray{T2,N2}) where {T1,M1,W,N1,T2,N2} = transpose(p1)*p2
#
# Base.LinAlg.A_mul_Bt(p1::PolyMatrix{T1,M1,Val{W},N1},
#   p2::AbstractArray{T2,N2}) where {T1,M1,W,N1,T2,N2} = p1*transpose(p2)
#
# Base.LinAlg.At_mul_Bt(p1::PolyMatrix{T1,M1,Val{W},N1},
#   p2::AbstractArray{T2,N2}) where {T1,M1,W,N1,T2,N2} = transpose(p1)*transpose(p2)
#
# # ctranspose
# Base.LinAlg.Ac_mul_B(p1::PolyMatrix{T1,M1,Val{W},N1},
#   p2::AbstractArray{T2,N2}) where {T1,M1,W,N1,T2,N2} = ctranspose(p1)*p2
#
# Base.LinAlg.A_mul_Bc(p1::PolyMatrix{T1,M1,Val{W},N1},
#   p2::AbstractArray{T2,N2}) where {T1,M1,W,N1,T2,N2} = p1*ctranspose(p2)
#
# Base.LinAlg.Ac_mul_Bc(p1::PolyMatrix{T1,M1,Val{W},N1},
#   p2::AbstractArray{T2,N2}) where {T1,M1,W,N1,T2,N2} = ctranspose(p1)*ctranspose(p2)
#
# # multiplication between Polynomialmatrices
# # transpose
# Base.LinAlg.At_mul_B(p1::PolyMatrix{T1,M1,Val{W},N1},
#   p2::PolyMatrix{T2,M2,Val{W},N2}) where {T1,M1,W,N1,T2,M2,N2} = transpose(p1)*p2
#
# Base.LinAlg.A_mul_Bt(p1::PolyMatrix{T1,M1,Val{W},N1},
#   p2::PolyMatrix{T2,M2,Val{W},N2}) where {T1,M1,W,N1,T2,M2,N2} = p1*transpose(p2)
#
# Base.LinAlg.At_mul_Bt(p1::PolyMatrix{T1,M1,Val{W},N1},
#   p2::PolyMatrix{T2,M2,Val{W},N2}) where {T1,M1,W,N1,T2,M2,N2} = transpose(p1)*transpose(p2)
#
# # ctranspose
# Base.LinAlg.Ac_mul_B(p1::PolyMatrix{T1,M1,Val{W},N1},
#   p2::PolyMatrix{T2,M2,Val{W},N2}) where {T1,M1,W,N1,T2,M2,N2} = ctranspose(p1)*p2
#
# Base.LinAlg.A_mul_Bc(p1::PolyMatrix{T1,M1,Val{W},N1},
#   p2::PolyMatrix{T2,M2,Val{W},N2}) where {T1,M1,W,N1,T2,M2,N2} = p1*ctranspose(p2)
#
# Base.LinAlg.Ac_mul_Bc(p1::PolyMatrix{T1,M1,Val{W},N1},
#   p2::PolyMatrix{T2,M2,Val{W},N2}) where {T1,M1,W,N1,T2,M2,N2} = ctranspose(p1)*ctranspose(p2)

# determinant
function det(p::PolyMatrix{T,M,Val{W},N}) where {T,M,W,N}
  size(p,1) == size(p,2) || throw(DimensionMismatch("the polynomial matrix must be square for a determinant to be defined"))
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
  return Polynomial(ar, W)
end

function _truncate(::Type{T}, a::AbstractArray{T2}) where {T<:Real,T2}
  r = similar(a, T)
  r = real(a)
end

function _truncate(::Type{T}, a::AbstractArray{T2}) where {T<:Integer,T2}
  r = similar(a, T)
  for i in eachindex(a)
    r[i] = convert(T, round(real(a[i])))
  end
  r
end

# inversion
# return determinant polynomial and adjugate polynomial matrix
function inv(p::PolyMatrix{T,M,Val{W},N}) where {T,M,W,N}
  size(p,1) == size(p,2) || throw(DimensionMismatch("the polynomial matrix must be square for a determinant to be defined"))
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
  rdet = Polynomial(ar, W)

  v2  = zero(B)
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

  return rdet, PolyMatrix(r, size(r), Val{W})
end

function _detrange(i,n)
  a = 1:i-1
  b = i+1:n
  r = vcat(collect(a), collect(b))
  return ifelse(length(r) > 1, r, r[1]:r[1])
end
