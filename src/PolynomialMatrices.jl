module PolynomialMatrices

using DataStructures
using Polynomials
using Compat

# Import functions for overloading
import Polynomials: coeffs, degree, variable
import Base: start, next, done
import Base: promote_rule, convert, size, length
import Base: +, -, *, /, inv, det
import Base: getindex, setindex!, linearindexing, eltype
import Base: copy
import Base: transpose, ctranspose
import Base: summary
import Base: insert!
import Base: checkbounds
import Base: filt!, filt
import Compat.view
import Base: vecnorm, norm, rank, isapprox, ==, isequal, hash

# Export
export PolyMatrix
export coeffs, degree, variable
export vartype
export col_degree, row_degree
export high_col_deg_matrix, high_row_deg_matrix
export is_col_proper, is_row_proper
export colred, rowred
export ltriang, rtriang, hermite

typealias SymbolLike  Union{Symbol,AbstractString,Char}
typealias ForwardOrdering Base.Order.ForwardOrdering

# Include files
include("type.jl")
include("methods.jl")
include("conversions.jl")
include("arithmetic.jl")
include("reductions.jl")
include("filt.jl")

end # module
