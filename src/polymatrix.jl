ordertype{K,D,O}(::Type{SortedDict{K,D,O}}) = O
ordertype{K,D,O}(s::SortedDict{K,D,O}) = O

# Parameters:
#   T: type of the polynomials' coefficients
#   M: type of the matrices of coefficients
#   O: ordering type to store the matrices of coefficients
#   N: number of dimensions of polynomial matrix
#
# NOTE: The constructors should make sure that the highest coefficient matrix is nonzero
# NOTE: For dense polynomial matrices, defining `coeffs` as a 3-D array could be much more efficient
immutable PolyMatrix{T,M,O,N} <: AbstractArray{Polynomials.Poly{T},N}
  coeffs::SortedDict{Int,M,O}
  dims::NTuple{N,Int}
  var::Symbol

  @compat function (::Type{PolyMatrix}){M,O,N}(
      coeffs::SortedDict{Int,M,O}, dims::NTuple{N,Int}, var::Symbol=:x)
    T = eltype(M)
    new{T,M,O,N}(coeffs, dims, var)
  end
end

# Evaluation of a polynomial matrix at a specific value x
@compat function (p::PolyMatrix{T,M,O,N}){T,M,O,N,S}(x::S)
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
  PolyMatrix(c, size(PM), var)
end

function PolyMatrix{M<:AbstractArray,N}(A::M, dims::NTuple{N,Int}, var::Symbol=:x)
  c  = SortedDict(Dict{Int,M}())
  ny = dims[1]
  m  = div(size(A,1), ny)
  if size(A,1) != m*ny || size(A,2) != dims[2]
    warn("PolyMatrix: dimensions are not consistent")
    throw(DomainError())
  end
  for k = 0:m-1
    p = view(A, k*ny+(1:ny), :)
    insert!(c, k, p)
  end
  return PolyMatrix(c, dims, var)
end

function PolyMatrix{T<:Number}(A::AbstractArray{T}, var::Symbol=:x)
  if ndims(A) > 2
    warn("PolyMatrix: higher order arrays not supported at this point")
    throw(DomainError())
  end
  c = SortedDict(Dict{Int,typeof(A)}())
  insert!(c, 0, A)
  return PolyMatrix(c, size(A), var)
end

promote_rule{T1,T2,M2,O,N}(::Type{PolyMatrix{T1,Array{T1,N},O,N}},
  ::Type{PolyMatrix{T2,M2,O,N}}) =
  PolyMatrix{promote_type(T1, T2), promote_type(M1,M2), O, N}

function _convert{T1,N,T2,M1,M2,O}(::Type{PolyMatrix{T1,M1,O,N}},
  p::PolyMatrix{T2,M2,O,N})
  r = PolyMatrix( SortedDict(Dict{Int,AbstractArray{T1,N}}()), size(p), p.var)
  for (k,c) in coeffs(p)
    r.coeffs[k] = map(x->convert(T1,x),c)
  end
  r
end
@generated function convert{T1,N,T2,M1,M2,O}(::Type{PolyMatrix{T1,M1,O,N}},
  p::PolyMatrix{T2,M2,O,N})
  if !(M1 <: AbstractArray{T1,N})
    return :(error("convert: first two parameters of first argument are incompatible"))
  elseif T1 == T2 && M1 == M2
    return :(p)
  else
    return :(_convert(PolyMatrix{T1,M1,O,N}, p::PolyMatrix{T2,M2,O,N}))
  end
end

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

summary{T,M,O,N}(p::PolyMatrix{T,M,O,N}) =
  string(Base.dims2string(p.dims), " PolyArray{$T,$N}")
