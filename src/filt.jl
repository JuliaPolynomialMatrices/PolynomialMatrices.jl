function _zerosi(b::PolyMatrix{S}, a::PolyMatrix{G}, T) where {S,G}
  m  = max(degree(a), degree(b))
  si = zeros(promote_type(S, G, T), size(a,1), m)
end

function filt(b::PolyMatrix{T,M1,W,N}, a::PolyMatrix{S,M2,W,N},
  x::AbstractArray{G}, si=_zerosi(b, a, G)) where {T,S,M1,M2,W,N,G}
  filt!(Array{promote_type(T, G, S)}(undef, size(a,1), size(x,2)), b, a, x, si)
end

function filt!(out::AbstractArray{H}, b::PolyMatrix{T,M1,W,N},
  a::PolyMatrix{S,M2,W,N}, x::AbstractArray{G}, si=_zerosi(b, a, G)) where {H,T,S,M1,M2,W,N,G}

  as = degree(a)
  bs = degree(b)
  sz = max(as, bs)
  silen = sz
  if size(si, 2) < silen
    throw(ArgumentError("initial state vector si must have max(degree(a),degree(b))-1 columns"))
  end

  ac = coeffs(a)
  bc = coeffs(b)
  # Filter coefficient normalization
  if !haskey(coeffs(a), 0)
    throw(ArgumentError("First coefficient matrix must be nonzero of a"))
  end
  a0 = ac[0]
  if a0 != Matrix{Float64}(I,size(a0))
    a = inv(a0)*a
    b = inv(a0)*b
  end

  # check for fir case
  if as == 0
    if bs == 0
      # simple scaling
      out = bc[0]*x
    else
      _filt_fir!(out, b, x, si)
    end
    return out
  end

  if bs == 0 && bc[0]==zero(bc[0])
    _filt_ar!(out, a, x, si)
  else
    _filt_iir!(out, b, a, x, si)
  end
  return out
end

function _filt_iir!(out::AbstractArray{T}, b::PolyMatrix{T,M1,W,N},
  a::PolyMatrix{T,M2,W,N}, x::AbstractArray{S}, si::AbstractArray{G}) where {T,S,M1,M2,W,N,G}
  silen = size(si,2)
  bc = coeffs(b)
  ac = coeffs(a)
  v0b = zero(similar(bc[first(keys(bc))]))
  v0a = zero(similar(ac[first(keys(ac))]))
  for i = 0:silen
    get!(bc, i, v0b) # inserts nonexisting entries as zeros
  end
  for i = 0:silen
    get!(ac, i, v0a) # inserts nonexisting entries as zeros
  end

  val = zeros(promote_type(G,S,T), size(a,1), 1)
  @inbounds @simd for i=1:size(x, 2)
    xi  = x[:,i:i]
    val = si[:,1:1] + bc[0]*xi
    for j=1:(silen-1)
       si[:,j] = si[:,j+1:j+1] + bc[j]*xi - ac[j]*val
    end
    si[:,silen] = bc[silen]*xi - ac[silen]*val
    out[:,i] = val
  end

  # remove nonexisting entries again
  truncate!(b)
  truncate!(a)
end

function _filt_fir!(
  out::AbstractMatrix{T}, b::PolyMatrix{T,M1,W,N}, x, si=zeros(T, size(b,1), degree(b))) where {T,M1,W,N}
  silen = size(si,2)
  bc = coeffs(b)
  v0 = zero(similar(bc[first(keys(bc))]))

  # inserts nonexisting entries as zeros
  for i = 0:silen
    get!(bc, i, v0)
  end

  @inbounds @simd for i=1:size(x, 2)
    xi  = x[:,i:i]
    val = si[:,1:1] + bc[0]*xi
    for j=1:(silen-1)
      si[:,j] = si[:,j+1:j+1] + bc[j]*xi
    end
    si[:,silen] = bc[silen]*xi
    out[:,i] = val
  end

  # remove nonexisting entries again
  truncate!(b)
end

function _filt_ar!(
  out::AbstractMatrix{T}, a::PolyMatrix{T,M1,W,N},
  x::AbstractArray{T}, si=zeros(T, size(a,1), degree(a))) where {T,M1,W,N}
  silen = size(si,2)
  ac = coeffs(a)
  val = zeros(T, size(a,1), 1)
  v0 = zero(similar(ac[first(keys(ac))]))

  # inserts nonexisting entries as zeros
  for i = 0:silen
    get!(ac, i, v0)
  end

  @inbounds @simd for i=1:size(x, 2)
    xi = view(x,:,i)
    val = si[:,1] + xi
    for j=1:(silen-1)
       si[:,j] = si[:,j+1] - ac[j]*val
    end
    si[:,silen] = - ac[silen]*val
    out[:,i] = val
  end

  # remove nonexisting entries again
  truncate!(a)
end
