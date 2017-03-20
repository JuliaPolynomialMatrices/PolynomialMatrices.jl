module trial

using Compat

immutable MyType{T,M,V} # <: AbstractMatrix{T}
  A::M

  @compat function (::Type{MyType}){M<:AbstractMatrix,V}(A::M, ::Type{Val{V}})
    new{eltype(A),M,V}(A)
  end
end

function MyType{M<:AbstractVector,V}(A::M, ::Type{Val{V}}; b::Bool=false)
  if !b
    return MyType(reshape(A,size(A,1,2)...), Val{V})
  end
end

Base.size(t::MyType) = size(t.A)
Base.length(t::MyType) = length(t.A)
Base.getindex(t::MyType, i::Int) = getindex(t.A, i)

export MyType

end

using trial

MyType(randn(4), Val{:x})

A = randn(4)
reshape(A,size(A,1,2)...)


truncate nonzero coeffs
fix assumption that first exist
rank

using PolynomialMatrices
using DataStructures


using BenchmarkTools

A = randn(1000,1000)
B = randn(1000,1000)
@benchmark $A.'

@benchmark $A*$B.'
@benchmark $A.'*$B

### Description

Some of the implemented and yet to be implemented functionality

- [x] Core PolynomialMatrix
  - [x] basic PolynomialMatrix interface
  - [x] use SortedDict for internal representation
- [ ] Math operations
  - [x] +
  - [x] -
  - [x] *
  - [ ] /
- [x] Indexing:
  - [x] getindex
  - [x] setindex!
  - [x] Iteration
- [x] Constructor
  - [x] constructor from matrix of Polynomials.Poly
  - [x] constructor form dict of AbstractMatrix
  - [x] constructor from tall AbstractMatrix
  - [x] constructor from 3-dimensional AbstractMatrix
- [x] Functions
  - [x] degree, rank, determinant, transpose, inv
  - [x] row and column reduction
  - [ ] gcd
- [ ] transformations
  - [x] triangularization
  - [x] transformation to hermite form
  - [ ] transformation to smith form
