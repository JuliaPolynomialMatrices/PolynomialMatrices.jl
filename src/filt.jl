function _zerosi{S,G}(b::PolyMatrix{S}, a::PolyMatrix{G}, T)
  m  = max(degree(a), degree(b))
  si = zeros(promote_type(S, G, T), size(a,1), m)
end

function filt{T,S,M1,M2,W,N,G}(b::PolyMatrix{T,M1,W,N}, a::PolyMatrix{S,M2,W,N},
  x::AbstractArray{G}, si=_zerosi(b, a, G))
  filt!(Array{promote_type(T, G, S)}(size(a,1), size(x,2)), b, a, x, si)
end

function filt!{H,T,S,M1,M2,W,N,G}(out::AbstractArray{H}, b::PolyMatrix{T,M1,W,N},
  a::PolyMatrix{S,M2,W,N}, x::AbstractArray{G}, si=_zerosi(b, a, G))

  as = degree(a)
  bs = degree(b)
  sz = max(as, bs)
  if as == 0
    if bs == 0
      # simple scaling
      out = (coeffs(a)[0]\coeffs(b)[0])*x
    else
      bc = coeffs(b)
      for i = 0:sz
        get!(bc, i, zeros(similar(bc[first(keys(bc))]))) # inserts nonexisting entries as zeros
      end
      _filt_fir!(out, b, x, si)
    end
    return out
  end
  sz = max(as, bs)
  silen = sz
  if size(si, 2) != silen
      throw(ArgumentError("initial state vector si must have max(length(a),length(b))-1 columns"))
  end

  # Filter coefficient normalization TODO
  a0 = coeffs(a)[0]
  if a0 != eye(a0)
    a = inv(a0)*a
    b = inv(a0)*b
  end
  ac = coeffs(a)
  bc = coeffs(b)
  v0b = zeros(similar(bc[first(keys(bc))]))
  v0a = zeros(similar(ac[first(keys(ac))]))
  for i = 0:sz
    get!(bc, i, v0b) # inserts nonexisting entries as zeros
  end
  for i = 0:sz
    get!(ac, i, v0a) # inserts nonexisting entries as zeros
  end

  _filt_iir!(out, b, a, x, si)
  return out
end

function _filt_iir!{T,S,M1,M2,W,N,G}(out::AbstractArray{T}, b::PolyMatrix{T,M1,W,N},
  a::PolyMatrix{T,M2,W,N}, x::AbstractArray{S}, si::AbstractArray{G})
  silen = size(si,2)
  bc = coeffs(b)
  ac = coeffs(a)
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
end

function _filt_fir!{T,M1,W,N}(
  out::AbstractMatrix{T}, b::PolyMatrix{T,M1,W,N}, x, si=zeros(T, size(b,1), degree(b)))
  silen = size(si,2)
  bc = coeffs(b)
  v0 = zeros(similar(bc[first(keys(bc))]))

  # inserts nonexisting entries as zeros
  for i = 0:degree(b)
    get!(bc, i, v0)
  end

  @inbounds @simd for i=1:size(x, 2)
    xi = view(x,:,i)
    val = si[:,1] + bc[0]*xi
    for j=1:(silen-1)
      si[:,j] = si[:,j+1] + bc[j]*xi
    end
    si[:,silen] = bc[silen]*xi
    out[:,i] = val
  end
end

function _filt_ar!{T,M1,W,N}(
  out::AbstractMatrix{T}, a::PolyMatrix{T,M1,W,N},
  x::AbstractArray{T}, si=zeros(T, size(a,1), degree(a)))
  silen = size(si,2)
  ac = coeffs(a)
  val = zeros(T, size(a,1), 1)
  v0 = zeros(similar(ac[first(keys(ac))]))

  # inserts nonexisting entries as zeros
  for i = 0:degree(a)
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
end
