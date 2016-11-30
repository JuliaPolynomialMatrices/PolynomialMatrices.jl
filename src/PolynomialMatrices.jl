module PolynomialMatrices

using DataStructures
using Polynomials
using Compat

# Import functions for overloading
import Polynomials.coeffs
import Base: start, next, done
import Base: promote_rule, convert, size, length
import Base: +,-,*
import Base: getindex
import Base: transpose, ctranspose
import Base: show, print, showcompact
import Base: insert!

# Export
export PolyMatrix
export order

# Include files
include("polymatrix.jl")
include("arithmetic.jl")


end # module
