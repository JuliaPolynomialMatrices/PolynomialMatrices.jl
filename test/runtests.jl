using PolynomialMatrices
using Base.Test
using Polynomials

# Include tests heres
# include ("test1.jl")

# test constructors
p1  = Poly([1])
p2  = Poly([2,1,3])
p3  = Poly([2,3,4])
m   = [p1 p2; p2 p1]
pm1 = PolyMatrix(m)

degreepm2 = 8
ny  = 2
nu  = 2
A   = randn(ny*(degreepm2+1),nu)
B   = eye(Float64,2)

@test_throws DomainError PolyMatrix(A, (5,nu))
@test_throws DomainError PolyMatrix(A, (ny,3))
@test_throws DomainError PolyMatrix(randn(3,3,3))
pm2 = PolyMatrix(A, (ny,nu))
pm3 = PolyMatrix(B)

# test equality and isapprox
p1  = Poly([2,1,3.], :s)
p2  = Poly([2,1,3.1], :s)
p3  = Poly([Inf,1,3.1], :s)
p4  = Poly([Inf,1.,3], :s)
p5  = Poly([NaN64,1.,3], :s)
pm1 = PolyMatrix([p1 p2; p2 p1])
pm2 = PolyMatrix([p1 p2; p2 p2])
pm3 = PolyMatrix([p3 p2; p2 p3])
pm4 = PolyMatrix([p4 p2; p2 p4])
pm5 = PolyMatrix([p5 p2; p2 p4])

@test pm1 != pm2 != pm3 != pm4 != pm5

@test !isapprox(pm3,pm4; rtol=0.001)
@test isapprox(pm3,pm4; rtol=0.1)

B = [2 2; 2 2]
C = [1 1; 1 1]
PolyB = PolyMatrix(B, (2,2))
@test B == PolyB
@test C != PolyB
@test pm2 == copy(pm2)
@test pm1 != pm2
@test pm1 != PolyB

# different variables
B = [0. 1; 1 1]
C = [-0. 1; 1 1]
PolyB = PolyMatrix(B, (2,2), :s)
PolyC = PolyMatrix(B, (2,2), :q)
PolyD = PolyMatrix(C, (2,2), :s)

@test PolyB ≠ PolyC && !isequal(PolyB, PolyC)
@test_throws DomainError PolyB ≈ PolyC
@test PolyB == PolyD && PolyB ≈ PolyD && !isequal(PolyB, PolyD)

# test copy
p1  = Poly([1])
p2  = Poly([2,1,3])
pm1 = PolyMatrix([p1 p2; p2 p1])
pm2 = copy(pm1)

@test isequal(pm1, pm2) && !is(pm1, pm2)

# test getindex
p1  = Poly([1])
p2  = Poly([2,1,3])
p3  = Poly([2,3,4])
m   = [p1 p2; p2 p1]
pm1 = PolyMatrix(m)
degreepm2 = 8
A   = randn(ny*(degreepm2+1),nu)
B   = eye(Float64,2)

pm2 = PolyMatrix(A, (ny,nu))
pm3 = PolyMatrix(B)
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

# test iteration
for elem in pm2
  @test degree(elem) == degreepm2
end

for idx in eachindex(pm3)
  @test degree(pm3[idx]) == 0
end

@test coeffs(pm3[end]) ≈ coeffs(one(Poly{Float64}))

# test transpose and ctranspose
pm4[2] = p3
pm5 = transpose(pm4)
@test coeffs(pm4[2]) ≈ coeffs(pm5[3])

C   = randn(2,2) + randn(2,2)im
pm6 = PolyMatrix(C)
@test coeffs(ctranspose(pm6))[0] ≈ ctranspose(C)

# test filtering
N = 100
x = randn(ny,N)
out = similar(x)

filt(pm2,pm2,x)
filt(pm3,pm2,x)
filt(pm2,pm3,x)
