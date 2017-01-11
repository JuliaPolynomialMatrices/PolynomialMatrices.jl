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

  @compat function (::Type{PolyMatrix}){T}(P::AbstractArray{Polynomials.Poly{T}})
    length(size(P)) <= 2 ||
      error("PolyMatrix: higher order arrays not supported at this point")
    maxdegree = maximum(map(degree, P))
    var = countnz(P) > 0 ? P[findfirst(x -> x != zero(x), P)].var :
                           variable(Polynomials.Poly{T}).var
    M = typeof(similar(P, T)) # NOTE: Is there a more memory-efficient way to obtain M?
    coeffs = SortedDict(Dict{Int,M}(i => convert(M,[p[i] for p in P]) for i in 0:maxdegree))
    new{T,M,ordertype(coeffs),ndims(P)}(coeffs, size(P), var)
  end

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
function PolyMatrix{T,N}(A::AbstractArray{T}, dims::NTuple{N,Int}, var::Symbol=:x)
  coeffs = SortedDict(Dict{Int,typeof(A)}())
  ny = dims[1]
  m = div(size(A,1), ny)
  for k = 0:m-1
    p = view(A, k*ny+(1:ny),:)
    insert!(coeffs, k, p)
  end
  return PolyMatrix(coeffs,dims,var)
end
function PolyMatrix{T<:Number}(A::AbstractArray{T}, var::Symbol=:x)
  length(size(A)) <= 2 ||
    error("PolyMatrix: higher order arrays not supported at this point")
  coeffs = SortedDict(Dict{Int,typeof(A)}())
  insert!(coeffs, 0, A)
  return PolyMatrix(coeffs,size(A),var)
end

promote_rule{T1,T2,M1,M2,O,N}(::Type{PolyMatrix{T1,M1,O,N}}, ::Type{PolyMatrix{T2,M2,O,N}}) =
  PolyMatrix{promote_type(T1, T2), promote_type(M1,M2), O, N}

function convert{T1,T2,M1,M2,O,N}(::Type{PolyMatrix{T1,M1,O,N}}, p::PolyMatrix{T2,M2,O,N})
  r = PolyMatrix( SortedDict(Dict{Int,M1}()), p.dims, p.var)
  return r+p
end
convert{T,M,O,N}(::Type{PolyMatrix{T,M,O,N}}, p::PolyMatrix{T,M,O,N}) = p

size(p::PolyMatrix) = p.dims
function size(p::PolyMatrix, i::Int)
  i > 0 || throw(ArgumentError("size: dimension needs to be postive"))
  return i <= length(p.dims) ? p.dims[i] : 1
end

length{T,M,O,N}(p::PolyMatrix{T,M,O,N})      = prod(p.dims)
start{T,M,O,N}(p::PolyMatrix{T,M,O,N})       = 1
next{T,M,O,N}(p::PolyMatrix{T,M,O,N}, state) = p[state], state+1
done{T,M,O,N}(p::PolyMatrix{T,M,O,N}, state) = state > length(p)
Base.linearindexing{T,M,O,N}(::Type{PolyMatrix{T,M,O,N}}) = Base.LinearFast()
#(p::PolyMatrix)       = start(coeffs(p))
#next(p::PolyMatrix, state) = next(coeffs(p), state)
#done(p::PolyMatrix, state) = done(coeffs(p), state)

# Copying
function copy{T,M}(p::PolyMatrix{T,M})
  r = PolyMatrix( SortedDict(Dict{Int,M}()), p.dims, p.var)
  for (k,v) in p.coeffs
    r.coeffs[k] = copy(v)
  end
  return r
end

# Slicing (`getindex`)
function getindex(p::PolyMatrix, row::Int, col::Int)
  1 ≤ row ≤ size(p,1) || error("s[idx,]: idx out of bounds")
  1 ≤ col ≤ size(p,2) || error("s[,idx]: idx out of bounds")
  Poly([v[row, col] for (k,v) in p.coeffs],p.var)
end

function getindex{T,M,O,N}(p::PolyMatrix{T,M,O,N}, idx::Int)
  1 ≤ idx ≤ length(p) || error("p[idx]: idx out of bounds")
  Poly([v[idx] for (k,v) in p.coeffs],p.var)
end

getindex(p::PolyMatrix, ::Colon)           = p
getindex(p::PolyMatrix, ::Colon, ::Colon)  = p

function getindex{T,M,O,N}(p::PolyMatrix{T,M,O,N}, ::Colon, idx::Int)
  1 ≤ idx ≤ length(p) || error("p[idx]: idx out of bounds")
  r = PolyMatrix( SortedDict(Dict{Int,M}()), p.dims, p.var)
  for (k,v) in p.coeffs
    r.coeffs[k] = reshape(v[:,idx], p.dims[1], 1)
  end
  return r
end

function getindex{T,M,O,N}(p::PolyMatrix{T,M,O,N}, idx::Int, ::Colon)
  1 ≤ idx ≤ length(p) || error("p[idx]: idx out of bounds")
  r = PolyMatrix( SortedDict(Dict{Int,M}()), p.dims, p.var)
  for (k,v) in p.coeffs
    r.coeffs[k] = reshape(v[idx,:], 1, p.dims[2])
  end
  return r
end

function insert!{T,M,O,N}(p::PolyMatrix{T,M,O,N}, k::Int, A)
  insert!(coeffs(p), k, A)
end

# Obtain dictionary of coefficient matrices
coeffs(p::PolyMatrix) = p.coeffs

# Maximum degree of the polynomials in a polynomial matrix
degree{T,M,O,N}(p::PolyMatrix{T,M,O,N})  = last(coeffs(p))[1]

function transpose{T,M<:AbstractMatrix,O,N}(p::PolyMatrix{T,M,O,N})
  r = PolyMatrix( SortedDict(Dict{Int,M}()), (p.dims[2], p.dims[1]), p.var)
  for (k,v) in p.coeffs
    r.coeffs[k] = v.'
  end
  return r
end
function transpose{T,M,O}(p::PolyMatrix{T,M,O,1})
  r = PolyMatrix( SortedDict(Dict{Int,AbstractArray{T,2}}()), (1,p.dims[1]), p.var)
  for (k,v) in p.coeffs
    r.coeffs[k] = v.'
  end
  return r
end

function ctranspose{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  r = PolyMatrix( SortedDict(Dict{Int,M}()), (p.dims[2], p.dims[1]), p.var)
  for (k,v) in p.coeffs
    r.coeffs[k] = v'
  end
  return r
end
function ctranspose{T,M,O}(p::PolyMatrix{T,M,O,1})
  r = PolyMatrix( SortedDict(Dict{Int,AbstractArray{T,2}}()), (1,p.dims[1]), p.var)
  for (k,v) in p.coeffs
    r.coeffs[k] = v'
  end
  return r
end

# NOTE: It would be nicer if "; " were replaced by linebreaks
# the "Poly" word were removed, and the entries of each row were
# right justified (as when Julia prints e.g. "[1 100; 100 0]")

summary{T,M,O,N}(p::PolyMatrix{T,M,O,N}) =
 string(Base.dims2string(p.dims), " PolyMatrix{$T}")
