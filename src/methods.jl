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

## Obtain dictionary of coefficient matrices
coeffs(p::PolyMatrix) = p.coeffs

## Maximum degree of the polynomials in a polynomial matrix
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

# TODO is there a definition for matrix norms for polynomial matrices?
# TODO Lₚ norms
# vecnorm
# stacks all coefficient matrices and calls built in vecnorm on resulting tall matrix
function vecnorm{T1,M1,O,N}(p1::PolyMatrix{T1,M1,O,N}, p::Real=2)
  c = coeffs(p1)
  vecnorm(vcat(values(c)...), p)
end

## Comparison
=={T1,M1,O,N,T2,M2}(p1::PolyMatrix{T1,M1,O,N}, p2::PolyMatrix{T2,M2,O,N}) = (p1.var == p2.var && p1.coeffs == p2.coeffs)
==(p1::PolyMatrix, p2::PolyMatrix) = false

function =={T1,M1,O,N,M2<:AbstractArray}(p1::PolyMatrix{T1,M1,O,N}, n::M2)
  c = coeffs(p1)
  has_zero = false
  for (k,v) in c
    if k == 0
      v == n || return false
      has_zero = true
    else
      v == zeros(v) || return false
    end
  end
  return ifelse(has_zero, true, n == zeros(n))
end
=={T1,M1,O,N,M2<:AbstractArray}(n::M2, p1::PolyMatrix{T1,M1,O,N}) = (p1 == n)

hash(p::PolyMatrix, h::UInt) = hash(p.var, hash(coeffs(p), h))
isequal(p1::PolyMatrix, p2::PolyMatrix) = hash(p1) == hash(p2)

function isapprox{T1,M1,O,N,T2,M2}(p1::PolyMatrix{T1,M1,O,N}, p2::PolyMatrix{T2,M2,O,N};
  rtol::Real=Base.rtoldefault(T1,T2), atol::Real=0, norm::Function=vecnorm)
  p1.var == p2.var || throw(DomainError())
  d = norm(p1 - p2)
  if isfinite(d)
    return d <= atol + rtol*max(norm(p1), norm(p2))
  else
    c1 = coeffs(p1)
    c2 = coeffs(p2)
    sᵢ = intersect(keys(c1),keys(c2))
    s₁ = setdiff(keys(c1),sᵢ)
    s₂ = setdiff(keys(c2),sᵢ)
    for k in s₁
      isapprox(c1[k], zeros(c1[k]); rtol=rtol, atol=atol) || return false
    end
    for k in s₂
      isapprox(c2[k], zeros(c2[k]); rtol=rtol, atol=atol) || return false
    end
    for k in sᵢ
      isapprox(c1[k], c2[k]; rtol=rtol, atol=atol)        || return false
    end
    return true
  end
end

function isapprox{T1,M1,O,N,M2<:AbstractArray}(p1::PolyMatrix{T1,M1,O,N}, n::M2;
  rtol::Real=Base.rtoldefault(T1,T2), atol::Real=0, norm::Function=vecnorm)
  p1.var == p2.var || throw(DomainError())
  d = norm(p1 - n)
  if isfinite(d)
    return d <= atol + rtol*max(norm(p1), norm(n))
  else
    c = coeffs(p1)
    has_zero = false
    for (k,v) in c
      if k == 0
        isapprox(v, n; rtol=rtol, atol=atol) || return false
        has_zero = true
      else
        isapprox(v, zeros(v); rtol=rtol, atol=atol) || return false
      end
    end
    return ifelse(has_zero, true, isapprox(n,zeros(n)))
  end
end
isapprox{T1,M1,O,N,M2<:AbstractArray}(n::M2, p1::PolyMatrix{T1,M1,O,N}) = (p1 == n)

summary{T,M,O,N}(p::PolyMatrix{T,M,O,N}) =
  string(Base.dims2string(p.dims), " PolyArray{$T,$N}")
