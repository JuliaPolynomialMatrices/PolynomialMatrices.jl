__precompile__(true)

module PolynomialMatrices

using DataStructures
using Polynomials
using Compat
using LinearAlgebra
using DSP
using SparseArrays
using FFTW

# Import functions for overloading
import Polynomials: coeffs, degree, variable
import Base: promote_rule, convert, size, length
import Base: +, -, *, /, inv
import Base: getindex, setindex!, eltype, similar
import Base: copy
import Base: transpose, adjoint
import Base: summary
import Base: insert!
import Base: checkbounds
import Compat.view
import Base: isapprox, ==, isequal, hash
import LinearAlgebra: det, norm, rank
import DSP: filt!, filt

# Export
export PolyMatrix
export coeffs, degree, variable
export vartype, mattype
export col_degree, row_degree
export high_col_deg_matrix, high_row_deg_matrix
export is_col_proper, is_row_proper
export colred, rowred
export ltriang, rtriang, hermite
export gcrd, gcld
export fastrank

const SymbolLike = Union{Symbol,AbstractString,Char}
const ForwardOrdering = Base.Order.ForwardOrdering

# Include files
include("type.jl")
include("methods.jl")
include("conversions.jl")
include("arithmetic.jl")
include("reductions.jl")
include("filt.jl")

end # module
