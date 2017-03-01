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
pm2 = PolyMatrix(A, (ny,nu))
pm3 = PolyMatrix(B)

A = randn(3,2,5)
PolyMatrix(A, :s)
A = randn(3,2)
PolyMatrix(A, :s)
A = randn(3)
PolyMatrix(A, :s)

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

!isapprox(pm3,pm4; rtol=0.001)

B = [2 2; 2 2]
C = [1 1; 1 1]
PolyB = PolyMatrix(B, (2,2), :s)
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
p1  = Poly([1],:s)
p2  = Poly([2,1,3],:s)
p3  = Poly([2,3,4],:s)
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
@test typeof(pm1[1:2,2]) == PolyMatrix{Int,Vector{Int},Val{vartype(pm1)},1}

t = pm1[1:2,2:2]
@test typeof(t) == PolyMatrix{Int,Matrix{Int},Val{vartype(pm1)},2}

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

# test inverse
p1  = Poly([1.0,2.0,3.5],:s)
p2  = Poly([2,1.1,3,4],:s)
p3  = Poly([2,3.1,4,5,7.3],:s)
m   = [p1 p2; p3 p1]
pm1 = PolyMatrix(m)

det1, adj1 = inv(pm1)
t1 = adj1*pm1
@test norm(t1[2,1])/norm(t1[1,1]) < eps(Float64)

# triangularization
s = variable("s")
p = PolyMatrix([s-1 s^2-1; 2*one(s) 2s+2; zero(s) 3*one(s)])
L,U = ltriang(p, false, 1)
L₀ = PolyMatrix([-1.225s+1.225 zero(s); -2.450*one(s) zero(s); -1.225*one(s) 1.732*one(s)])
@test isapprox(p*U, L₀; rtol=1e-3)
@test isapprox(L, L₀; rtol=1e-3)
@test isapprox(p*U, L)

R,U = rtriang(p.', false, 1)
@test isapprox(R.', L₀; rtol=1e-3)
@test isapprox(PolyMatrix(U*p.'), R)

L,U = ltriang(p)
@test isapprox(p*U, L₀; rtol=1e-3)
@test isapprox(L, L₀; rtol=1e-3)
@test isapprox(p*U, L)

R,U = rtriang(p.')
@test isapprox(R.', L₀; rtol=1e-3)
@test isapprox(PolyMatrix(U*p.'), R)

# hermite
s = variable("s")
p = PolyMatrix([-s^3-2s^2+1 -(s+1)^2; (s+2)^2*(s+1) zero(s)])
H₀ = PolyMatrix([s+1 zero(s); (s+2)^2*(s+1) (s+2)^2*(s+1)^2])
U₀ = PolyMatrix([one(s) s+1; -s -s^2-s+1])

H,U = hermite(p)
@test isapprox(p*U, H₀)
@test isapprox(H, H₀)
@test isapprox(U, U₀)

# colred
s = variable("s")

# example 1 from "A Fortran 77 package for column reduction of polynomial matrices" Geurts, A.J. Praagman, C., 1998
p = PolyMatrix([s^4+6s^3+13s^2+12s+4 -s^3-4s^2-5s-2; zero(s) s+2])
U₀ = PolyMatrix([one(s) zero(s); s+2 one(s)])
R₀ = PolyMatrix([zero(s) -(s^3+4s^2+5s+2); s^2+4s+4 s+2])
R,U = colred(p)
@test isapprox(R, R₀)
@test isapprox(U, U₀)
@test isapprox(p*U, R)

R,U = rowred(p.')
@test isapprox(R, R₀.')
@test isapprox(U, U₀.')
@test isapprox(PolyMatrix(U*p.'), R)

# example 2 from "A Fortran 77 package for column reduction of polynomial matrices" Geurts, A.J. Praagman, C., 1998
p = PolyMatrix([s^4 s^2 s^6+1; s^2 one(s) s^4; one(s) zero(s) one(s)])
U₀ = PolyMatrix([one(s) zero(s); s+2 one(s)])
R₀ = PolyMatrix([zero(s) -(s^3+4s^2+5s+2); s^2+4s+4 s+2])
#R,U = colred(p)       # TODO this should work ?!
#@test isapprox(p*U, R)
#@test degree(R) == 0

# R,U = rowred(p.')
#@test isapprox(PolyMatrix(U*p.'), R)

# example 3 from "A Fortran 77 package for column reduction of polynomial matrices" Geurts, A.J. Praagman, C., 1998
ϵ = 0.001
p = PolyMatrix([s^3+s^2 ϵ*s+1 one(s); 2s^2 -one(s) -one(s); 3s^2 one(s) one(s)])
R,U = colred(p)
@test isapprox(p*U, R)

R,U = rowred(p.')
@test isapprox(PolyMatrix(U*p.'), R)

# example 4 from "A Fortran 77 package for column reduction of polynomial matrices" Geurts, A.J. Praagman, C., 1998
ϵ = e-8
p = PolyMatrix([s^3+s^2+2s+1 ϵ*s^2+2s+3 s^2+s+1   s-1;
                s-1          -s+2       2s^2+s-1  2s+1;
                s+3          2s-1       -s^2-2s+1 -s-2;
                one(s)       -one(s)    3s+1       3*one(s)])
#R,U = colred(p) # TODO this should work ?!
# @test isapprox(p*U, R)

#R,U = rowred(p.')
#@test isapprox(PolyMatrix(U*p.'), R)

# left2right matrix fractional descriptions
s = variable("s")
Nₗ = PolyMatrix([s^2 zero(s); -4s s])
Dₗ = PolyMatrix([s^3+2s^2-1 s+1; -5s^2-13s-8 (s+1)*(s+4)])
Nᵣ = PolyMatrix([-s^2 -s; zero(s) -s])
Dᵣ = PolyMatrix([-s^3-2s^2+1 -(s+1)^2; (s+2)^2*(s+1) zero(s)])

rmfd = PolyMatrix(vcat(Dᵣ,Nᵣ))
lmfd = PolyMatrix(hcat(-Nₗ, Dₗ))

# verify that example is correct.
@test vecnorm(lmfd*rmfd) ≈ 0

L,U = ltriang(lmfd)

N = U[3:4,3:4]
D = U[1:2,3:4]

# compare hermite form
Dₕ,U = hermite(D)
Nₕ = N*U

D₀,U = hermite(Dᵣ)
N₀ = Nᵣ*U

@test isapprox(Dₕ,D₀)
@test isapprox(Nₕ,N₀)
