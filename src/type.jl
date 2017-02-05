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
