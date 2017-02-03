size(p::PolyMatrix) = p.dims
function size(p::PolyMatrix, i::Int)
  i > 0 || throw(ArgumentError("size: dimension needs to be positive"))
  return i ≤ length(p.dims) ? p.dims[i] : 1
end

length{T,M,O,N}(p::PolyMatrix{T,M,O,N})      = prod(size(p))
start{T,M,O,N}(p::PolyMatrix{T,M,O,N})       = 1
next{T,M,O,N}(p::PolyMatrix{T,M,O,N}, state) = p[state], state+1
done{T,M,O,N}(p::PolyMatrix{T,M,O,N}, state) = state > length(p)
linearindexing{T<:PolyMatrix}(::Type{T})     = Base.LinearFast()
eltype{T,M,O,N}(p::PolyMatrix{T,M,O,N})      = Poly{T}

# Copying
function copy{T,M}(p::PolyMatrix{T,M})
  r = PolyMatrix( SortedDict(Dict{Int,M}()), size(p), p.var)
  for (k,v) in coeffs(p)
    r.coeffs[k] = copy(v)
  end
  return r
end

# getindex
function getindex{T,M,O,N}(p::PolyMatrix{T,M,O,N}, i::Int)
  @compat @boundscheck checkbounds(p, i)
  r = Poly([v[i] for (k,v) in p.coeffs],p.var)
  return r
end

_PolyMatrix{T,N}(p::Array{Poly{T},N}) = PolyMatrix(p)
_PolyMatrix{T}(p::Poly{T})            = p

# Hacking AbstractArrays getindex to make sure we return PolyMatrix for non-scalar return types
function getindex(A::PolyMatrix, I...)
    r = Base._getindex(linearindexing(A), A, I...)
    _PolyMatrix(r)
end

# setindex!
function setindex!{T,M,O,N,U}(Pm::PolyMatrix{T,M,O,N}, p::Poly{U}, i::Int)
  @compat @boundscheck checkbounds(Pm, i)
  c = coeffs(p)
  Pmc = coeffs(Pm)
  S = Set(eachindex(coeffs(p)))
  for (k, v) in Pmc
    if k+1 ∈ S
      v[i] = c[k+1]
      delete!(S,k+1)
    else
      v[i] = zero(T)
    end
    # NOTE should we delete a key if all elements are made zero?
    # if all(v .== zero(T))
    #   delete!(Pmc, idx-1)
    # end
  end
  # add keys not previously present in Pm
  for k in S
    vk = spzeros(T, Pm.dims...)
    vk[i] = c[k]
    insert!(Pmc, k-1, vk)
  end
end

function insert!{T,M,O,N}(p::PolyMatrix{T,M,O,N}, k::Int, A)
  insert!(coeffs(p), k, A)
end

# Obtain dictionary of coefficient matrices
coeffs(p::PolyMatrix) = p.coeffs

# Maximum degree of the polynomials in a polynomial matrix
degree{T,M,O,N}(p::PolyMatrix{T,M,O,N})  = last(coeffs(p))[1]

function transpose{T,M<:AbstractMatrix,O,N}(p::PolyMatrix{T,M,O,N})
  r = PolyMatrix( SortedDict(Dict{Int,M}()), reverse(p.dims), p.var)
  for (k,v) in p.coeffs
    r.coeffs[k] = transpose(v)
  end
  return r
end

function ctranspose{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  r = PolyMatrix( SortedDict(Dict{Int,M}()), reverse(p.dims), p.var)
  for (k,v) in p.coeffs
    r.coeffs[k] = ctranspose(v)
  end
  return r
end

"""
    isapprox{T,S}(p1::Poly{T}, p2::Poly{S}; reltol::Real = Base.rtoldefault(T,S), abstol::Real = 0, norm::Function = vecnorm)
Truncate `p1` and `p2`, and compare the coefficients with `isapprox`.
The tolerances `reltol` and `abstol` are passed to both `truncate` and `isapprox`.
"""
function isapprox{T1,M1,O,N,T2,M2}(p1::PolyMatrix{T1,M1,O,N}, p2::PolyMatrix{T2,M2,O,N};
  rtol::Real = Base.rtoldefault(T1,T2), atol::Real = 0, norm::Function = vecnorm)
  p1.var == p2.var || throw(DomainError())
  p1t = truncate(p1; rtol=rtol, atol=atol)
  p2t = truncate(p2; rtol=rtol, atol=atol)
  c1 = coeffs(p1t)
  c2 = coeffs(p2t)
  length(c1) == length(p2) || throw(DomainError())
  for ((k1,v1), (k2,v2)) in zip(c1,c2)
    k1 == k2                                          || return false
    isapprox(v1, v2; rtol=rtol, atol=atol, norm=norm) || return false
  end
  return true
end

function isapprox{T,M,O,N,S<:Array}(p1::PolyMatrix{T,M,O,N}, n::S; reltol::Real = Base.rtoldefault(T,S),
  abstol::Real = 0)
  p1t = truncate(p1; reltol = reltol, abstol = abstol)
  degree(p1t) == 0 && isapprox(coeffs(p1), [n]; rtol = reltol, atol = abstol)
end

isapprox{T,M,O,N,S<:Number}(n::S, p1::PolyMatrix{T,M,O,N}; reltol::Real = Base.rtoldefault(T,S),
  abstol::Real = 0) = isapprox(p1, n; reltol = reltol, abstol = abstol)

hash(p::PolyMatrix, h::UInt) = hash(p.var, hash(coeffs(p), h))
isequal(p1::PolyMatrix, p2::PolyMatrix) = hash(p1) == hash(p2)

summary{T,M,O,N}(p::PolyMatrix{T,M,O,N}) =
  string(Base.dims2string(p.dims), " PolyArray{$T,$N}")
