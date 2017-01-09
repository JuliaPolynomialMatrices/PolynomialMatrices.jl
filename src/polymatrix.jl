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
immutable PolyMatrix{T,M,O,N}
  coeffs::SortedDict{Int,M,O}
  dims::NTuple{N,Int}
  var::Symbol

  @compat function (::Type{PolyMatrix}){T}(P::AbstractArray{Polynomials.Poly{T}})
    @assert length(size(P)) <= 2 "PolyMatrix: higher order arrays not supported at this point"
    maxorder = maximum(map(degree, P))
    var = P[findfirst(x -> x != Poly(0) && x != 0, P)].var
    M = typeof(similar(P, T)) # NOTE: Is there a more memory-efficient way to obtain M?
    coeffs = SortedDict(Dict{Int,M}(i => convert(M,[p[i] for p in P]) for i in 0:maxorder))
    new{T,M,ordertype(coeffs),ndims(P)}(coeffs, size(P), var)
  end

  @compat function (::Type{PolyMatrix}){M,O,N}(
      coeffs::SortedDict{Int,M,O}, dims::NTuple{N,Int}, var::Symbol=:x)
    T = eltype(M)
    new{T,M,O,N}(coeffs, dims, var)
  end
end

# Evaluation of a polynomial matrix at a specific value x
# NOTE: Can this be turned into a Horner scheme?
@compat function (p::PolyMatrix{T,M,O,N}){T,M,O,N,S}(x::S)
  length(p) == 0 && return spzeros(T,dims...)

  c  = p.coeffs
  R  = similar(dims -> zeros(promote_type(T,S), dims), indices(c[findfirst(c)]))
  kvec = sort(collect(keys(c)), rev=true)
  show(kvec)
  R += c[kvec[1]]

  length(p) == 1 && return R*x^kvec[1]
  kp = kvec[1]
  for k in kvec[2:end]
    println(kp-k)
    show(R)
    R  = c[k] + R*x^(kp-k)
    kp = k
  end
  return R
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
  @assert length(size(A)) <= 2 "PolyMatrix: higher order arrays not supported at this point"
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

# Total number of entries of a polynomial matrix
length(p::PolyMatrix) = *(p.dims...)

start(p::PolyMatrix)       = 1
next(p::PolyMatrix, state) = p[state], state+1
done(p::PolyMatrix, state) = state > length(p)
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
  @assert 1 ≤ row ≤ size(p,1) "s[idx,]: idx out of bounds"
  @assert 1 ≤ col ≤ size(p,2) "s[,idx]: idx out of bounds"
  Poly([v[row, col] for (k,v) in p.coeffs],p.var)
end

function getindex(p::PolyMatrix, idx::Int)
  @assert 1 ≤ idx ≤ length(p) "p[idx]: idx out of bounds"
  Poly([v[idx] for (k,v) in p.coeffs],p.var)
end

getindex(p::PolyMatrix, ::Colon)           = p
getindex(p::PolyMatrix, ::Colon, ::Colon)  = p

function getindex{T,M,O,N}(p::PolyMatrix{T,M,O,N}, ::Colon, idx::Int)
  @assert 1 ≤ idx ≤ length(p) "p[idx]: idx out of bounds"
  r = PolyMatrix( SortedDict(Dict{Int,M}()), p.dims, p.var)
  for (k,v) in p.coeffs
    r.coeffs[k] = reshape(v[:,idx], p.dims[1], 1)
  end
  return r
end

function getindex{T,M,O,N}(p::PolyMatrix{T,M,O,N}, idx::Int, ::Colon)
  @assert 1 ≤ idx ≤ length(p) "p[idx]: idx out of bounds"
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
# NOTE: Shouldn't this be called "degree" instead?
order(p::PolyMatrix)  = last(coeffs(p))[1]

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

# Printing functions
function show(io::IO, p::PolyMatrix)
  print(io, '[')
  print(io, p)
  print(io, ']')
end

# NOTE: It would be nicer if "; " were replaced by linebreaks
# the "Poly" word were removed, and the entries of each row were
# right justified (as when Julia prints e.g. "[1 100; 100 0]")
function print{T}(io::IO, P::PolyMatrix{T})
  firstrow = true
  for i in 1:size(P,1)
    firstrow ? nothing : print(io, "; ")
    firstrow = false
    for j in 1:size(P,2)
      j == size(P,2) ? print(io, P[i,j]) : print(io, P[i,j], " ")
    end
  end
end

showcompact(io::IO, p::PolyMatrix) = print(io, p)
