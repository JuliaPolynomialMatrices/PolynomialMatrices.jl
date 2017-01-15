using PolynomialMatrices
using Base.Test
using Polynomials

# Include tests heres
# include ("test1.jl")

p1  = Poly([1])
p2  = Poly([2,1,3])
p3  = Poly([2,3,4])
m   =  [p1 p2; p2 p1]
pm1 = PolyMatrix(m)

# test getindex
@test pm1[1].a ≈ p1.a
@test typeof(pm1[1]) == Poly{Int}

t = pm1[1:2,2]
@test coeffs(t[1]) ≈ coeffs(p2)
@test coeffs(t[2]) ≈ coeffs(p1)
@test typeof(pm1[1:2,2]) == PolyMatrix{Int,Vector{Int},Base.Order.ForwardOrdering,1}

t = pm1[1:2,2:2]
@test typeof(t) == PolyMatrix{Int,Matrix{Int},Base.Order.ForwardOrdering,2}

# test setindex!
@test_throws InexactError pm1[1] = Poly([1.5])
pm1[1] = p3
@test coeffs(pm1[1]) ≈ coeffs(p3)
