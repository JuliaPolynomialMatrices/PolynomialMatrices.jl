size(p::PolyMatrix) = p.dims
size(p::PolyMatrix, i::Int) = i ≤ length(p.dims) ? p.dims[i] : 1

length(p::PolyMatrix{T,M,Val{W},N}) where {T,M,W,N}                     = prod(size(p))
iterate(p::PolyMatrix{T,M,Val{W},N}) where {T,M,W,N}                    = (p[1], 1)
iterate(p::PolyMatrix{T,M,Val{W},N}, i::Int) where {T,M,W,N}        = i < lastindex(p) ? (p[i+1], i+1) : nothing
IteratorSize(p::PolyMatrix{T,M,Val{W},N}, i::Int) where {T,M,W,N}   = HasShape{N}()
IteratorEltype(p::PolyMatrix{T,M,Val{W},N}, i::Int) where {T,M,W,N} = HasEltype()
eltype(::PolyMatrix{T,M,Val{W},N}) where {T,M,W,N}                      = Polynomial{T}
vartype(p::PolyMatrix{T,M,Val{W},N}) where {T,M,W,N}                    = W
mattype(p::PolyMatrix{T,M,Val{W},N}) where {T,M,W,N}                    = M
@compat Base.IndexStyle(::Type{<:PolyMatrix})                           = IndexLinear()

function similar(p::PolyMatrix{T,M,Val{W},N}, dims::NTuple{N2,Int}) where {T,M,W,N,N2}
  _,v1 = coeffs(p) |> first
  vr = zero(similar(v1, T, dims))
  r = PolyMatrix(SortedDict(0=>vr), size(vr), Val{W})
end

function similar(p::PolyMatrix{T,M,Val{W},N}, dims::Int) where {T,M,W,N}
  _,v1 = coeffs(p) |> first
  vr = zero(similar(v1, T, dims))
  r = PolyMatrix(SortedDict(0=>vr), size(vr), Val{W})
end

function similar(p::PolyMatrix{T,M,Val{W},N}, ::Type{S}=T,
  dims::NTuple{N2,Int}=size(p)) where {T,M,W,N,S,N2}
  _,v1 = coeffs(p) |> first
  vr = zero(similar(v1, S, dims))
  r = PolyMatrix(SortedDict(0=>vr), size(vr), Val{W})
end

function similar(p::PolyMatrix{T,M,Val{W},N}, ::Type{Polynomial{S}},
  dims::NTuple{N2,Int}=size(p)) where {T,M,W,N,S,N2}
  _,v1 = coeffs(p) |> first
  vr = zero(similar(v1, S, dims))
  r = PolyMatrix(SortedDict(0=>vr), size(vr), Val{W})
end

function _matrixofpoly(p::PolyMatrix)
  [p[i,j] for i in 1:size(p,1), j in 1:size(p,2)]
end

function Base.vcat(A::PolyMatrix{T,M,Val{W},N}...) where {T,M,W,N}
  mvec  = map(_matrixofpoly, A)
  mpoly = vcat(mvec...)
  return PolyMatrix(mpoly, Val{W})
end

function Base.hcat(A::PolyMatrix{T,M,Val{W},N}...) where {T,M,W,N}
  mvec  = map(_matrixofpoly, A)
  mpoly = hcat(mvec...)
  return PolyMatrix(mpoly, Val{W})
end

# works but is not properly since @inferred do not give correct type
function Base.cat(catdims, A::PolyMatrix{T,M,Val{W},N}...) where {T,M,W,N}
  mvec  = map(_matrixofpoly, A)
  mpoly = cat(mvec...,dims=catdims)
  return PolyMatrix(mpoly, Val{W})
end

function Base.hvcat(nbc::Int, A::PolyMatrix{T,M,Val{W},N}...) where {T,M,W,N}
  mvec  = map(_matrixofpoly, A)
  mpoly = hvcat(nbc, mvec...)
  return PolyMatrix(mpoly, Val{W})
end

function Base.hvcat(rows::NTuple{N2,Int}, A::PolyMatrix{T,M,Val{W},N}...) where {T,M,W,N,N2}
  mvec  = map(_matrixofpoly, A)
  mpoly = hvcat(rows, mvec...)
  return PolyMatrix(mpoly, Val{W})
end

"""
      variable(p::PolyMatrix)
  return variable of `p` as a `Polynomial` object.
"""
variable(p::PolyMatrix{T,M,Val{W},N}) where {T,M,W,N} = variable(Polynomial{T}, @compat Symbol(W))

# Copying
function copy(p::PolyMatrix{T,M,Val{W},N}) where {T,M,W,N}
  cr  = SortedDict(Dict{Int,M}())
  for (k,v) in coeffs(p)
    insert!(cr, k, copy(v))
  end
  PolyMatrix(cr, size(p), Val{W})
end

# getindex
function Base.checkbounds(p::PolyMatrix{T,M,Val{W},N}, I...) where {T,M,W,N}
  checkbounds(first(coeffs(p))[2], I...)
end

function getindex(p::PolyMatrix{T,M,Val{W},N}, i::Int) where {T,M,W,N}
  @compat @boundscheck checkbounds(p, i)
  vr = zeros(T, degree(p)+1)
  for (k,v) in coeffs(p)
    vr[k+1] = v[i]
  end
  r = Polynomial(vr, W)
end

function getindex(p::PolyMatrix{T,M,Val{W},N}, i::Int, j::Int) where {T,M,W,N}
  vr = zeros(T, degree(p)+1)
  for (k,v) in coeffs(p)
    vr[k+1] = v[i,j]
  end
  r = Polynomial(vr, W)
end

function getindex(p::PolyMatrix{T,M,Val{W},N}, I::Vararg{Int, N2}) where {T,M,W,N,N2}
  vr = zeros(T, degree(p)+1)
  for (k,v) in coeffs(p)
    vr[k+1] = v[I]
  end
  r = Polynomial(vr, W)
end

function getindex(p::PolyMatrix{T,M,Val{W},N}, I...) where {T,M,W,N}
  c   = coeffs(p)
  _,v = first(c)
  vr  = getindex(v, I...)
  cr  = SortedDict(Dict{Int,typeof(vr)}())
  for (k,v) in c
    insert!(cr, k, getindex(v, I...))
  end
  PolyMatrix(cr, size(vr), Val{W})
end

# setindex!
function setindex!(Pm::PolyMatrix{T,M,Val{W},N}, p::Polynomial{U}, i::Int) where {T,M,W,N,U}
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

function setindex!(Pm::PolyMatrix{T,M,Val{W},2}, p::Polynomial{U},
    i::Int, j::Int) where {T,M,W,U}
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

function setindex!(Pm::PolyMatrix{T,M,Val{W},N}, p::Polynomial{U},
    I...) where {T,M,W,N,U}
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

function setindex!(Pm::PolyMatrix{T,M,Val{W},N}, p::AbstractArray{Polynomial{U}},
    I...) where {T,M,W,N,U}
  @compat @boundscheck checkbounds(Pm, I...)
  c   = coeffs(Pm)
  _,v = first(c)
  for (i,pol) in zip(eachindex(view(v,I...)),p)
    Pm[i] = pol
  end
end

# setindex for number
function setindex!(Pm::PolyMatrix{T,M,Val{W},N}, p::T2, i::Int) where {T,M,W,N,T2<:Number}
  @compat @boundscheck checkbounds(Pm, i)
  c = coeffs(Pm)
  hasconst = false
  for (k,v) in c
    v[i]      = k == 0 ? convert(T, p) : zero(T)
    hasconst |= k == 0
  end
  if !hasconst
    v0 = spzeros(T, Pm.dims...)
    v0[i] = p
    insert!(c, 0, v0)
  end
end

function setindex!(Pm::PolyMatrix{T,M,Val{W},2}, p::T2,
    i::Int, j::Int) where {T,M,W,T2<:Number}
  @compat @boundscheck checkbounds(Pm, i, j)
  c = coeffs(Pm)
  hasconst = false
  for (k,v) in c
    v[i,j]      = k == 0 ? convert(T, p) : zero(T)
    hasconst |= k == 0
  end
  if !hasconst
    v0 = spzeros(T, Pm.dims...)
    v0[i,j] = p
    insert!(c, 0, v0)
  end
end

function setindex!(Pm::PolyMatrix{T,M,Val{W},N}, p::T2, I::Vararg{Int, N2}) where {T,M,W,N,T2<:Number,N2}
  @compat @boundscheck checkbounds(Pm, I...)
  c = coeffs(Pm)
  hasconst = false
  for (k,v) in c
    v[I...] = k == 0 ? convert(T, p) : zero(T)
    hasconst |= k == 0
  end
  if !hasconst
    v0 = spzeros(T, Pm.dims...)
    v0[I...] = p
    insert!(c, 0, v0)
  end
end

function setindex!(Pm::PolyMatrix{T,M,Val{W},N}, p::T2,
    I...) where {T,M,W,N,T2<:Number}
  @compat @boundscheck checkbounds(Pm, I...)
  c = coeffs(Pm)
  hasconst = false
  for (k,v) in c
    v[I...]      = k == 0 ? convert(T, p) : zero(T)
    hasconst |= k == 0
  end
  if !hasconst
    v0 = spzeros(T, Pm.dims...)
    v0[I...] = p
    insert!(c, 0, v0)
  end
end

function setindex!(Pm::PolyMatrix{T,M,Val{W},N}, p::AbstractArray{T2},
    I...) where {T,M,W,N,U,T2<:Number}
  @compat @boundscheck checkbounds(Pm, I...)
  c   = coeffs(Pm)
  _,v = first(c)
  for (i,pol) in zip(eachindex(view(v,I...)),p)
    Pm[i] = pol
  end
end

function firstindex(p::PolyMatrix{T,M,Val{W},N}) where {T,M,W,N}
  c   = coeffs(p)
  _,v = first(c)
  firstindex(v)
end

function lastindex(p::PolyMatrix{T,M,Val{W},N}) where {T,M,W,N}
  c   = coeffs(p)
  _,v = first(c)
  lastindex(v)
end

# insert!
function insert!(p::PolyMatrix{T,M,Val{W},N}, k::Int, A) where {T,M,W,N}
  if size(A) != size(p)
    throw(DomainError((p,k,A),"coefficient matrix to insert does not have the same size as polynomial matrix"))
  end
  insert!(coeffs(p), k, convert(M,A))
end

## Obtain dictionary of coefficient matrices
coeffs(p::PolyMatrix) = p.coeffs

## Maximum degree of the polynomials in a polynomial matrix
degree(p::PolyMatrix{T,M,Val{W},N}) where {T,M,W,N}  = last(coeffs(p))[1]

function transpose(p::PolyMatrix{T,M,Val{W},N}) where {T,M<:AbstractMatrix,W,N}
  cr = SortedDict(Dict{Int,M}())
  for (k,v) in p.coeffs
    insert!(cr, k, transpose(v))
  end
  PolyMatrix(cr, reverse(p.dims), Val{W})
end

function adjoint(p::PolyMatrix{T1,M1,Val{W},N}) where {T1,M1,W,N}
  c     = coeffs(p)
  k1,v1 = first(c)
  vr    = adjoint(v1)
  M     = typeof(vr)
  cr    = SortedDict(Dict{Int,M}())
  for (k,v) in c
    insert!(cr, k, adjoint(v))
  end
  PolyMatrix(cr, size(p), Val{W})
end

# TODO is there a definition for matrix norms for polynomial matrices?
# TODO Lₚ norms
# norm
# stacks all coefficient matrices and calls built in norm on resulting tall matrix
function norm(p₁::PolyMatrix{T1,M1,Val{W},N}, p::Real=2) where {T1,M1,W,N}
  c = coeffs(p₁)
  norm(vcat(values(c)...), p)
end

## Comparison
==(p₁::PolyMatrix{T1,M1,Val{W},N}, p₂::PolyMatrix{T2,M2,Val{W},N}) where {T1,M1,W,N,T2,M2} = (p₁.coeffs == p₂.coeffs)
==(p₁::PolyMatrix{T1,M1,Val{W1},N}, p₂::PolyMatrix{T2,M2,Val{W2},N}) where {T1,M1,W1,W2,N,T2,M2} = false

function ==(p₁::PolyMatrix{T1,M1,Val{W},N}, n::M2) where {T1,M1,W,N,M2<:AbstractArray}
  c = coeffs(p₁)
  has_zero = false
  for (k,v) in c
    if k == 0
      v == n || return false
      has_zero = true
    else
      v == zero(v) || return false
    end
  end
  return ifelse(has_zero, true, n == zero(n))
end
==(n::M2, p₁::PolyMatrix{T1,M1,Val{W},N}) where {T1,M1,W,N,M2<:AbstractArray} = (p₁ == n)

hash(p::PolyMatrix{T1,M1,Val{W},N}, h::UInt) where {T1,M1,W,N} = hash(W, hash(coeffs(p), h))
isequal(p₁::PolyMatrix, p₂::PolyMatrix) = hash(p₁) == hash(p₂)

function isapprox(p₁::PolyMatrix{T1,M1,Val{W},N}, p₂::PolyMatrix{T2,M2,Val{W},N};
  rtol::Real=length(p₁)*degree(p₁)*degree(p₂)*Base.rtoldefault(T1,T2,0),
  atol::Real=0, norm::Function=norm) where {T1,M1,W,N,T2,M2}
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
      isapprox(c₁[k], zero(c₁[k]); rtol=rtol, atol=atol) || return false
    end
    for k in s₂
      isapprox(c₂[k], zero(c₂[k]); rtol=rtol, atol=atol) || return false
    end
    for k in sᵢ
      isapprox(c₁[k], c₂[k]; rtol=rtol, atol=atol)        || return false
    end
    return true
  end
end
function isapprox(p₁::PolyMatrix{T1,M1,Val{W1},N}, p₂::PolyMatrix{T2,M2,Val{W2},N};
  rtol::Real=Base.rtoldefault(T1,T2,0), atol::Real=0, norm::Function=norm) where {T1,M1,W1,W2,N,T2,M2}
  #@warn "p₁≈p₂: `p₁` ($T1,$W1) and `p₂` ($T2,$W2) have different variables"
  throw(DomainError((p₁,p₂),"the two polynomial matrices have different variables."))
end

function isapprox(p₁::PolyMatrix{T1,M1,Val{W},N},
  n::AbstractArray{T2,N}; rtol::Real=Base.rtoldefault(T1,T2,0), atol::Real=0,
  norm::Function=norm) where {T1,M1,W,N,T2}
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
        isapprox(v, zero(v); rtol=rtol, atol=atol, norm=norm) || return false
      end
    end
    return ifelse(has_zero, true, isapprox(n,zeros(n); rtol=rtol, atol=atol, norm=norm))
  end
end
isapprox(n::AbstractArray{T2,N}, p₁::PolyMatrix{T1,M1,Val{W},N}) where {T1,M1,W,N,T2} = (p₁ ≈ n)

function isapprox(p₁::PolyMatrix{T1,M1,Val{W},N},
  n::AbstractArray{Polynomial{T2},N}; rtol::Real=Base.rtoldefault(T1,T2,0),
  atol::Real=0, norm::Function=norm) where {T1,M1,W,N,T2}
  return isapprox(p₁, PolyMatrix(n); rtol=rtol, atol=atol, norm=norm)
end

function isapprox(n::AbstractArray{Polynomial{T2},N},
  p₁::PolyMatrix{T1,M1,Val{W},N}; rtol::Real=Base.rtoldefault(T1,T2,0),
  atol::Real=0, norm::Function=norm) where {T1,M1,W,N,T2}
  return isapprox(PolyMatrix(n), p₁; rtol=rtol, atol=atol, norm=norm)
end


function rank(p::PolyMatrix{T,M,Val{W},N}) where {T,M,W,N}
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

fastrank(p::PolyMatrix{T,M,Val{W},N}) where {T,M,W,N} = rank(p(randn(1)...))

summary(p::PolyMatrix{T,M,Val{W},N}) where {T,M,W,N} =
  string(Base.dims2string(p.dims), " PolyArray{$T,$N}")
