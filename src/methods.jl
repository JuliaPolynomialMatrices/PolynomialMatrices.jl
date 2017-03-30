size(p::PolyMatrix) = p.dims
size(p::PolyMatrix, i::Integer) = i ≤ length(p.dims) ? p.dims[i] : 1

length{T,M,W,N}(p::PolyMatrix{T,M,Val{W},N})      = prod(size(p))
start{T,M,W,N}(p::PolyMatrix{T,M,Val{W},N})       = 1
next{T,M,W,N}(p::PolyMatrix{T,M,Val{W},N}, state) = p[state], state+1
done{T,M,W,N}(p::PolyMatrix{T,M,Val{W},N}, state) = state > length(p)
eltype{T,M,W,N}(::PolyMatrix{T,M,Val{W},N})       = Poly{T}
vartype{T,M,W,N}(p::PolyMatrix{T,M,Val{W},N})     = W
mattype{T,M,W,N}(p::PolyMatrix{T,M,Val{W},N})     = M
@compat Base.IndexStyle(::Type{<:PolyMatrix})     = IndexLinear()

function similar{T,M,W,N,N2}(p::PolyMatrix{T,M,Val{W},N}, dims::NTuple{N2,Int})
  _,v1 = coeffs(p) |> first
  vr = zeros(similar(v1, T, dims))
  r = PolyMatrix(SortedDict(0=>vr), size(vr), Val{W})
end

function similar{T,M,W,N,S,N2}(p::PolyMatrix{T,M,Val{W},N}, ::Type{S}=T,
  dims::NTuple{N2,Int}=size(p))
  _,v1 = coeffs(p) |> first
  vr = zeros(similar(v1, S, dims))
  r = PolyMatrix(SortedDict(0=>vr), size(vr), Val{W})
end

function similar{T,M,W,N,S,N2}(p::PolyMatrix{T,M,Val{W},N}, ::Type{Poly{S}}=Poly{T},
  dims::NTuple{N2,Int}=size(p))
  _,v1 = coeffs(p) |> first
  vr = zeros(similar(v1, S, dims))
  r = PolyMatrix(SortedDict(0=>vr), size(vr), Val{W})
end

function _matrixofpoly(p::PolyMatrix)
  [p[i,j] for i in 1:size(p,1), j in 1:size(p,2)]
end

function Base.vcat{T,M,W,N}(A::PolyMatrix{T,M,Val{W},N}...)
  mvec  = map(_matrixofpoly, A)
  mpoly = vcat(mvec...)
  return PolyMatrix(mpoly, Val{W})
end

function Base.hcat{T,M,W,N}(A::PolyMatrix{T,M,Val{W},N}...)
  mvec  = map(_matrixofpoly, A)
  mpoly = hcat(mvec...)
  return PolyMatrix(mpoly, Val{W})
end

# works but is not properly since @inferred do not give correct type
function Base.cat{T,M,W,N}(catdims, A::PolyMatrix{T,M,Val{W},N}...)
  mvec  = map(_matrixofpoly, A)
  mpoly = cat(catdims, mvec...)
  return PolyMatrix(mpoly, Val{W})
end

function Base.hvcat{T,M,W,N}(nbc::Integer, A::PolyMatrix{T,M,Val{W},N}...)
  mvec  = map(_matrixofpoly, A)
  mpoly = hvcat(nbc, mvec...)
  return PolyMatrix(mpoly, Val{W})
end

function Base.hvcat{T,M,W,N,N2}(rows::NTuple{N2,Int}, A::PolyMatrix{T,M,Val{W},N}...)
  mvec  = map(_matrixofpoly, A)
  mpoly = hvcat(rows, mvec...)
  return PolyMatrix(mpoly, Val{W})
end

"""
      variable(p::PolyMatrix)
  return variable of `p` as a `Poly` object.
"""
variable{T,M,W,N}(p::PolyMatrix{T,M,Val{W},N}) = variable(T, @compat Symbol(W))

# Copying
function copy{T,M,W,N}(p::PolyMatrix{T,M,Val{W},N})
  cr  = SortedDict(Dict{Int,M}())
  for (k,v) in coeffs(p)
    insert!(cr, k, copy(v))
  end
  PolyMatrix(cr, size(p), Val{W})
end

# getindex
function Base.checkbounds{T,M,W,N}(p::PolyMatrix{T,M,Val{W},N}, I...)
  checkbounds(first(coeffs(p))[2], I...)
end

function getindex{T,M,W,N}(p::PolyMatrix{T,M,Val{W},N}, i::Integer)
  @compat @boundscheck checkbounds(p, i)
  vr = zeros(T, degree(p)+1)
  for (k,v) in coeffs(p)
    vr[k+1] = v[i]
  end
  r = Poly(vr, W)
end

function getindex{T,M,W,N}(p::PolyMatrix{T,M,Val{W},N}, i::Integer, j::Integer)
  vr = zeros(T, degree(p)+1)
  for (k,v) in coeffs(p)
    vr[k+1] = v[i,j]
  end
  r = Poly(vr, W)
end

function getindex{T,M,W,N}(p::PolyMatrix{T,M,Val{W},N}, I...)
  c   = coeffs(p)
  _,v = first(c)
  vr  = getindex(v, I...)
  cr  = SortedDict(Dict{Int,typeof(vr)}())
  for (k,v) in c
    insert!(cr, k, getindex(v, I...))
  end
  PolyMatrix(cr, size(vr), Val{W})
end

_PolyMatrix{T,N}(p::Array{Poly{T},N}) = PolyMatrix(p)
_PolyMatrix{T}(p::Poly{T})            = p

# setindex!
function setindex!{T,M,W,N,U}(Pm::PolyMatrix{T,M,Val{W},N}, p::Poly{U}, i::Integer)
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

function setindex!{T,M,W,U}(Pm::PolyMatrix{T,M,Val{W},2}, p::Poly{U},
    i::Integer, j::Integer)
  @compat @boundscheck checkbounds(Pm, i, j)
  c = coeffs(p)
  Pmc = coeffs(Pm)
  S = Set(eachindex(coeffs(p)))
  for (k, v) in Pmc
    if k+1 ∈ S
      v[i,j] = c[k+1]
      delete!(S,k+1)
    else
      v[i,j] = zero(T)
    end
    # NOTE should we delete a key if all elements are made zero?
    # if all(v .== zero(T))
    #   delete!(Pmc, idx-1)
    # end
  end
  # add keys not previously present in Pm
  for k in S
    vk = spzeros(T, Pm.dims...)
    vk[i,j] = c[k]
    insert!(Pmc, k-1, vk)
  end
end

function setindex!{T,M,W,N,U}(Pm::PolyMatrix{T,M,Val{W},N}, p::Poly{U},
    I...)
  @compat @boundscheck checkbounds(Pm, I...)
  c = coeffs(p)
  Pmc = coeffs(Pm)
  S = Set(eachindex(coeffs(p)))
  for (k, v) in Pmc
    if k+1 ∈ S
      v[I...] = c[k+1]
      delete!(S,k+1)
    else
      v[I...] = zero(T)
    end
    # NOTE should we delete a key if all elements are made zero?
    # if all(v .== zero(T))
    #   delete!(Pmc, idx-1)
    # end
  end
  # add keys not previously present in Pm
  for k in S
    vk = spzeros(T, Pm.dims...)
    vk[I...] = c[k]
    insert!(Pmc, k-1, vk)
  end
end

# setindex for number
function setindex!{T,M,W,N,T2<:Number}(Pm::PolyMatrix{T,M,Val{W},N}, p::T2, i::Integer)
  @compat @boundscheck checkbounds(Pm, i)
  c = coeffs(Pm)
  if haskey(c,0)
    _,v0 = c[0]
    v0[i] = p
    # NOTE should we delete a key if all elements are made zero?
    # if all(v .== zero(T))
    #   delete!(Pmc, idx-1)
    # end
  else
    v0 = spzeros(T, Pm.dims...)
    v0[i] = p
    insert!(c, 0, v0)
  end
end

function setindex!{T,M,W,T2<:Number}(Pm::PolyMatrix{T,M,Val{W},2}, p::T2,
    i::Integer, j::Integer)
  @compat @boundscheck checkbounds(Pm, i, j)
  c = coeffs(Pm)
  if haskey(c,0)
    v0      = c[0]
    v0[i,j] = p
    # NOTE should we delete a key if all elements are made zero?
    # if all(v .== zero(T))
    #   delete!(Pmc, idx-1)
    # end
  else
    v0 = spzeros(T, Pm.dims...)
    v0[i,j] = p
    insert!(c, 0, v0)
  end
end

function setindex!{T,M,W,N,T2<:Number}(Pm::PolyMatrix{T,M,Val{W},N}, p::T2,
    I...)
  @compat @boundscheck checkbounds(Pm, I...)
  c = coeffs(Pm)
  if haskey(c,0)
    _,v0 = c[0]
    v0[I...] = p
    # NOTE should we delete a key if all elements are made zero?
    # if all(v .== zero(T))
    #   delete!(Pmc, idx-1)
    # end
  else
    v0 = spzeros(T, Pm.dims...)
    v0[I...] = p
    insert!(c, 0, v0)
  end
end

function insert!{T,M,W,N}(p::PolyMatrix{T,M,Val{W},N}, k::Int, A)
  insert!(coeffs(p), k, A)
end

## Obtain dictionary of coefficient matrices
coeffs(p::PolyMatrix) = p.coeffs

## Maximum degree of the polynomials in a polynomial matrix
degree{T,M,W,N}(p::PolyMatrix{T,M,Val{W},N})  = last(coeffs(p))[1]

function transpose{T,M<:AbstractMatrix,W,N}(p::PolyMatrix{T,M,Val{W},N})
  cr = SortedDict(Dict{Int,M}())
  for (k,v) in p.coeffs
    insert!(cr, k, transpose(v))
  end
  PolyMatrix(cr, reverse(p.dims), Val{W})
end

function ctranspose{T1,M1,W,N}(p::PolyMatrix{T1,M1,Val{W},N})
  c     = coeffs(p)
  k1,v1 = first(c)
  vr    = ctranspose(v1)
  M     = typeof(vr)
  cr    = SortedDict(Dict{Int,M}())
  for (k,v) in c
    insert!(cr, k, ctranspose(v))
  end
  PolyMatrix(cr, size(p), Val{W})
end

# TODO is there a definition for matrix norms for polynomial matrices?
# TODO Lₚ norms
# vecnorm
# stacks all coefficient matrices and calls built in vecnorm on resulting tall matrix
function vecnorm{T1,M1,W,N}(p₁::PolyMatrix{T1,M1,Val{W},N}, p::Real=2)
  c = coeffs(p₁)
  vecnorm(vcat(values(c)...), p)
end

## Comparison
=={T1,M1,W,N,T2,M2}(p₁::PolyMatrix{T1,M1,Val{W},N}, p₂::PolyMatrix{T2,M2,Val{W},N}) = (p₁.coeffs == p₂.coeffs)
=={T1,M1,W1,W2,N,T2,M2}(p₁::PolyMatrix{T1,M1,Val{W1},N}, p₂::PolyMatrix{T2,M2,Val{W2},N}) = false

function =={T1,M1,W,N,M2<:AbstractArray}(p₁::PolyMatrix{T1,M1,Val{W},N}, n::M2)
  c = coeffs(p₁)
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
=={T1,M1,W,N,M2<:AbstractArray}(n::M2, p₁::PolyMatrix{T1,M1,Val{W},N}) = (p₁ == n)

hash{T1,M1,W,N}(p::PolyMatrix{T1,M1,Val{W},N}, h::UInt) = hash(W, hash(coeffs(p), h))
isequal(p₁::PolyMatrix, p₂::PolyMatrix) = hash(p₁) == hash(p₂)

function isapprox{T1,M1,W,N,T2,M2}(p₁::PolyMatrix{T1,M1,Val{W},N}, p₂::PolyMatrix{T2,M2,Val{W},N};
  rtol::Real=length(p₁)*degree(p₁)*degree(p₂)*Base.rtoldefault(T1,T2),
  atol::Real=0, norm::Function=vecnorm)
  d = norm(p₁ - p₂)
  if isfinite(d)
    return d <= atol + rtol*max(norm(p₁), norm(p₂))
  else
    c₁ = coeffs(p₁)
    c₂ = coeffs(p₂)
    sᵢ = intersect(keys(c₁),keys(c₂))
    s₁ = setdiff(keys(c₁),sᵢ)
    s₂ = setdiff(keys(c₂),sᵢ)
    for k in s₁
      isapprox(c₁[k], zeros(c₁[k]); rtol=rtol, atol=atol) || return false
    end
    for k in s₂
      isapprox(c₂[k], zeros(c₂[k]); rtol=rtol, atol=atol) || return false
    end
    for k in sᵢ
      isapprox(c₁[k], c₂[k]; rtol=rtol, atol=atol)        || return false
    end
    return true
  end
end
function isapprox{T1,M1,W1,W2,N,T2,M2}(p₁::PolyMatrix{T1,M1,Val{W1},N}, p₂::PolyMatrix{T2,M2,Val{W2},N};
  rtol::Real=Base.rtoldefault(T1,T2), atol::Real=0, norm::Function=vecnorm)
  warn("p₁≈p₂: `p₁` ($T1,$W1) and `p₂` ($T2,$W2) have different variables")
  throw(DomainError())
end

function isapprox{T1,M1,W,N,T2}(p₁::PolyMatrix{T1,M1,Val{W},N},
  n::AbstractArray{T2,N}; rtol::Real=Base.rtoldefault(T1,T2), atol::Real=0,
  norm::Function=vecnorm)
  d = norm(p₁ - n)
  if isfinite(d)
    return d <= atol + rtol*max(norm(p₁), norm(n))
  else
    c = coeffs(p₁)
    has_zero = false
    for (k,v) in c
      if k == 0
        isapprox(v, n; rtol=rtol, atol=atol) || return false
        has_zero = true
      else
        isapprox(v, zeros(v); rtol=rtol, atol=atol, norm=norm) || return false
      end
    end
    return ifelse(has_zero, true, isapprox(n,zeros(n); rtol=rtol, atol=atol, norm=norm))
  end
end
isapprox{T1,M1,W,N,T2}(n::AbstractArray{T2,N}, p₁::PolyMatrix{T1,M1,Val{W},N}) = (p₁ ≈ n)

function isapprox{T1,M1,W,N,T2}(p₁::PolyMatrix{T1,M1,Val{W},N},
  n::AbstractArray{Poly{T2},N}; rtol::Real=Base.rtoldefault(T1,T2),
  atol::Real=0, norm::Function=vecnorm)
  return isapprox(p₁, PolyMatrix(n); rtol=rtol, atol=atol, norm=norm)
end

function isapprox{T1,M1,W,N,T2}(n::AbstractArray{Poly{T2},N},
  p₁::PolyMatrix{T1,M1,Val{W},N}; rtol::Real=Base.rtoldefault(T1,T2),
  atol::Real=0, norm::Function=vecnorm)
  return isapprox(PolyMatrix(n), p₁; rtol=rtol, atol=atol, norm=norm)
end


function rank{T,M,W,N}(p::PolyMatrix{T,M,Val{W},N})
  d = degree(p)+1
  # copy all elements into three-dimensional matrix
  A = zeros(size(p)...,d)
  for (k,v) in coeffs(p)
    A[:,:,k+1] = v
  end
  B = fft(A,3)
  a = [rank(B[:,:,k]) for k = 1:d] # rank evaluated at fft frequencies
  return maximum(a)
end

fastrank{T,M,W,N}(p::PolyMatrix{T,M,Val{W},N}) = rank(p(randn(1)...))

summary{T,M,W,N}(p::PolyMatrix{T,M,Val{W},N}) =
  string(Base.dims2string(p.dims), " PolyArray{$T,$N}")
