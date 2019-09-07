# setup
p1  = Poly([1,2,3],:s)
p2  = Poly([2,1,3,4],:s)
p3  = Poly([2,3,4,5,7],:s)
m   = [p1 p2; p3 p1]
pm1 = PolyMatrix(m)
pm2 = PolyMatrix(m .+ 1.0)
pm3 = PolyMatrix(hcat(m,m))
pm4 = PolyMatrix(Matrix{Float64}(I,2,2),:x)
n1  = 1
n2  = 0.0
a1  = Matrix{Float64}(I,2,2)
a2  = Matrix{Float64}(I,2,4)
v1  = ones(2)
v2  = ones(4)

# addition
@test pm1+pm1 == PolyMatrix(2*m)
@test pm1+pm2 == PolyMatrix(2*m .+ 1)
@test_throws DimensionMismatch pm1+pm3
@test_throws ArgumentError pm1+pm4
@inferred pm1+pm1
@inferred pm2+pm1

@test pm1+a1 == PolyMatrix(m+a1)
@test a1+pm1 == PolyMatrix(a1+m)
@test_throws DimensionMismatch pm1+a2
@inferred pm1+a1
@inferred a1+pm1

@test pm1+p1 == PolyMatrix(m+fill(p1,2,2))
@test p1+pm1 == PolyMatrix(fill(p1,2,2)+m)
@inferred pm1+p1
@inferred p1+pm1

@test pm1+n1 == PolyMatrix(m+fill(n1,2,2)) # TODO should this be supported
@test n1+pm1 == PolyMatrix(fill(n1,2,2)+m)
@inferred pm1+n1
@inferred n1+pm1

# subtraction
@test pm1-pm1 == 0*PolyMatrix(m)
@test pm1-pm2 ≈ 0*PolyMatrix(m) - fill(1,2,2)
@test_throws DimensionMismatch pm1-pm3
@test_throws ArgumentError pm1-pm4
@inferred pm1-pm1
@inferred pm2-pm1

@test pm1-a1 == PolyMatrix(m-a1)
@test a1-pm1 == PolyMatrix(a1-m)
@test_throws DimensionMismatch pm1-a2
@inferred pm1-a1
@inferred a1-pm1

@test pm1-p1 == PolyMatrix(m-fill(p1,2,2))
@test p1-pm1 == PolyMatrix(fill(p1,2,2)-m)
@inferred pm1-p1
@inferred p1-pm1

@test pm1-n1 == PolyMatrix(m-fill(n1,2,2))
@test n1-pm1 == PolyMatrix(fill(n1,2,2)-m)
@inferred pm1-n1
@inferred n1-pm1

# multiplication
@test pm1*pm1 ≈ PolyMatrix(m*m)
@test pm1*pm2 ≈ PolyMatrix(m*(m+fill(1,2,2)))
@test_throws DimensionMismatch pm3*pm1
@test_throws ArgumentError pm1*pm4
#@inferred pm1*pm1
#@inferred pm2*pm1

@test pm1*a1 ≈ PolyMatrix(m*a1)
@test a1*pm1 ≈ PolyMatrix(a1*m)
@test_throws DimensionMismatch a2*pm1
#@inferred pm1*a1
#@inferred a1+pm1

# fft multiplication
n = 2
d = 25
pmf = PolyMatrix(randn(n,n,d), (n,n,d), :s)
pf  = Poly(randn(d), :s)
@test eltype(pmf*pmf) == Poly{Float64}
@test eltype(pmf*pf)  == Poly{Float64}
@test pf*pmf  == pmf*pf

@test pm1*v1 ≈ PolyMatrix(m*v1)
@test_throws DimensionMismatch pm1*v2

@test (pm1*p1)[1] ≈ m[1]*p1 # Polynomials is problematic
@test (p1*pm1)[1] ≈ m[1]*p1

@test pm1*n1 ≈ PolyMatrix(m)
@test n1*pm1 ≈ PolyMatrix(m)

# division
@test (pm1/(2n1))[1] ≈ m[1]/2

# inverse
#det2, adj2 = inv(pm2) TODO
t2 = adj2*pm2
#@test norm(t2[2,1])/norm(t2[1,1]) < Base.rtoldefault(Float64) TODO

# determinant
@test typeof(det(pm1)) == Poly{Int}
