# Parameters:
#   T: type of the polynomials' coefficients
#   M: type of the matrices of coefficients
#   O: ordering type to store the matrices of coefficients
#   N: number of dimensions of polynomial matrix
#
# NOTE: The constructors should make sure that the highest coefficient matrix is nonzero
# NOTE: For dense polynomial matrices, defining `coeffs` as a 3-D array could be much more efficient
immutable PolyMatrix{T,M,V,N} <: AbstractArray{Polynomials.Poly{T},N}
  coeffs::SortedDict{Int,M,ForwardOrdering}
  dims::NTuple{N,Int}

  @compat function (::Type{PolyMatrix}){M,N,V}(
      coeffs::SortedDict{Int,M,ForwardOrdering}, dims::NTuple{N,Int}, ::Type{Val{V}})
    T = eltype(M)
    new{T,M,Val{V},N}(coeffs, dims)
  end
end

# Evaluation of a polynomial matrix at a specific value x
@compat function (p::PolyMatrix{T,M,Val{V},N}){T,M,V,N,S}(x::S)
  degree(p) == 0 && return convert(M, zeros(T,size(p)...))*zero(S)

  c    = p.coeffs
  kvec = collect(keys(c)) #keys assumed sorted
  k    = pop!(kvec)

  # Horner scheme
  R = copy(c[k]) * one(S)
  for k in reverse(kvec)
    R = R*x + c[k]
  end
  R
end

# Outer constructor
function PolyMatrix{M<:AbstractArray}(d::Dict{Int,M}, var::SymbolLike=:x)
  PolyMatrix(d, Val{@compat Symbol(var)})
end

function PolyMatrix{M<:AbstractArray,V}(d::Dict{Int,M}, var::Type{Val{V}})
  if length(d) ≤ 0
    warn("PolyMatrix: lengt(d) == 0")
    throw(DomainError())
  end
  c = SortedDict(d)
  n,m = size(first(c)[2])
  for (k,v) in c
    if size(v) != (n,m)
      warn("PolyMatrix: size of elements not consistent")
      throw(DomainError())
    end
  end
  PolyMatrix(c, (n,m), Val{V})
end

function PolyMatrix{M1<:AbstractArray}(PM::M1)
  eltype(M1) <: Poly   || error("PolyMatrix: Matrix of Polynomials expected, try PolyMatrix(A, dims[, var])")
  length(size(PM)) <= 2 || error("PolyMatrix: higher order arrays not supported at this point")
  T = eltype(eltype(PM))
  M = typeof(similar(PM, T)) # NOTE: Is there a more memory-efficient way to obtain M?
  c = SortedDict(Dict{Int,M}())
  # find the union of all index sets of all polynomials in the matrix PM
  S = Set{Int}()
  for p in PM
    for elem in eachindex(coeffs(p))
      push!(S, elem)
    end
  end
  # initialize all elements to zero
  for idx in S
    insert!(c, idx-1, zeros(T, size(PM)...)) # TODO change to spzeros when Julia drops support for 0.4.7 (there are no sparse vectors)
  end
  # copy all elements
  for pidx in eachindex(PM)
    pc = coeffs(PM[pidx])
    for eidx in eachindex(pc)
      c[eidx-1][pidx] = pc[eidx]
    end
  end
  var = countnz(PM) > 0 ? PM[findfirst(x -> x != zero(x), PM)].var :
                         Poly(T[]).var       # default to Polys default variable
  PolyMatrix(c, size(PM), Val{@compat Symbol(var)})
end

# TODO should we simplify these constructors somehow
# Note: There is now copy in the following constructors
function PolyMatrix{T<:Number}(A::AbstractArray{T}, var::SymbolLike=:x)
  return PolyMatrix(A, size(A), Val{@compat Symbol(var)})
end

function PolyMatrix{T<:Number, V}(A::AbstractArray{T}, ::Type{Val{V}})
  return PolyMatrix(A, size(A), Val{V})
end

# TODO should we support vectors or only matrices
function PolyMatrix{M<:AbstractArray}(A::M, dims::Tuple{Int}, var::SymbolLike=:x)
  PolyMatrix(A, dims, Val{@compat Symbol(var)})
end

function PolyMatrix{M<:AbstractArray,V}(A::M, dims::Tuple{Int}, ::Type{Val{V}})
  ny = dims[1]
  dn = div(size(A,1), ny)
  if rem(size(A,1), ny) != 0
    warn("PolyMatrix: dimensions are not consistent")
    throw(DomainError())
  end
  p0 = dn > 0 ? p0 = A[1:ny] : zeros(eltype(A), dims)
  c  = SortedDict(Dict{Int,typeof(p0)}())
  insert!(c, 0, p0)
  for k = 1:dn-1
    p = A[k*ny+(1:ny)]
    insert!(c, k, p)
  end
  return PolyMatrix(c, dims, Val{V})
end

function PolyMatrix{M<:AbstractArray}(A::M, dims::Tuple{Int,Int}, var::SymbolLike=:x; reverse::Bool=false)
  PolyMatrix(A, dims, Val{@compat Symbol(var)}; reverse=reverse)
end

function PolyMatrix{M<:AbstractArray,V}(A::M, dims::Tuple{Int,Int}, ::Type{Val{V}}; reverse::Bool=false)
  ny = dims[1]
  dn = div(size(A,1), ny)
  if rem(size(A,1), ny) != 0 || size(A,2) != dims[2]
    warn("PolyMatrix: dimensions are not consistent")
    throw(DomainError())
  end
  p0 = dn > 0 ? A[1:ny, :] : zeros(eltype(A),dims)
  c  = SortedDict(Dict{Int,typeof(p0)}())
  for k = 0:dn-1
    idx = reverse ? (dn-k-1)*ny+(1:ny) : k*ny+(1:ny)
    v = A[idx, :]
    insert!(c, k, v)
  end
  return PolyMatrix(c, dims, Val{V})
end

function PolyMatrix{M<:AbstractArray}(A::M, dims::Tuple{Int,Int,Int}, var::SymbolLike=:x)
  PolyMatrix(A, dims, Val{@compat Symbol(var)})
end

function PolyMatrix{M<:AbstractArray,V}(A::M, dims::Tuple{Int,Int,Int}, ::Type{Val{V}})
  if size(A) != dims && dims[3] < 1
    warn("PolyMatrix: dimensions are not consistent")
    throw(DomainError())
  end
  dn = dims[3]
  p0 = A[:, :, 1]
  c  = SortedDict(Dict{Int,typeof(p0)}())
  insert!(c, 0, p0)
  for k = 1:dn-1
    p = A[:, :, k+1]
    if (sum(abs2,p) > zero(eltype(A))) insert!(c, k, p) end
  end
  return PolyMatrix(c, size(A,1,2), Val{V})
end
