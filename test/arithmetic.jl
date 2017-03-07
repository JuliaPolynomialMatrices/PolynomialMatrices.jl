# setup
p1  = Poly([1.0,2.0,3.5],:s)
p2  = Poly([2,1.1,3,4],:s)
p3  = Poly([2,3.1,4,5,7.3],:s)
m   = [p1 p2; p3 p1]
pm1 = PolyMatrix(m)
pm2 = PolyMatrix(m+1)
pm3 = PolyMatrix(hcat(m,m))
pm4 = PolyMatrix(eye(2),:x)
n1  = 1
n2  = 0.0
a1  = eye(2,2)
a2  = eye(2,4)
v1  = ones(2)

# addition
@test pm1+pm1 == PolyMatrix(2*m)
@test pm1+pm2 == PolyMatrix(2*m+1)
@test_throws DomainError pm1+pm3
@test_throws DomainError pm1+pm4

@test pm1+a1 == PolyMatrix(m+a1)
@test a1+pm1 == PolyMatrix(a1+m)
@test_throws DomainError pm1+a2

@test pm1+p1 == PolyMatrix(m+fill(p1,2,2))
@test p1+pm1 == PolyMatrix(fill(p1,2,2)+m)

@test pm1+n1 == PolyMatrix(m+fill(n1,2,2))
@test n1+pm1 == PolyMatrix(fill(n1,2,2)+m)

# subtraction
@test pm1-pm1 == 0*PolyMatrix(m)
@test pm1-pm2 == 0*PolyMatrix(m) - fill(1,2,2)
@test_throws DomainError pm1-pm3
@test_throws DomainError pm1-pm4

@test pm1-a1 == PolyMatrix(m-a1)
@test a1-pm1 == PolyMatrix(a1-m)
@test_throws DomainError pm1-a2

@test pm1-p1 == PolyMatrix(m-fill(p1,2,2))
@test p1-pm1 == PolyMatrix(fill(p1,2,2)-m)

@test pm1-n1 == PolyMatrix(m-fill(n1,2,2))
@test n1-pm1 == PolyMatrix(fill(n1,2,2)-m)

# multiplication
@test pm1*pm1 ≈ PolyMatrix(m*m)
@test pm1*pm2 ≈ PolyMatrix(m*(m+fill(1,2,2)))
@test_throws DomainError pm3*pm1
@test_throws DomainError pm1*pm4

@test pm1*a1 ≈ PolyMatrix(m*a1)
@test a1*pm1 ≈ PolyMatrix(a1*m)
@test_throws DomainError a2*pm1

@test pm1*v1 ≈ PolyMatrix(m*v1)

@test (pm1*p1)[1] ≈ m[1]*p1 # Polynomials is problematic
@test (p1*pm1)[1] ≈ m[1]*p1

@test pm1*n1 ≈ PolyMatrix(m)
@test n1*pm1 ≈ PolyMatrix(m)

# division
@test (pm1/(2n1))[1] ≈ m[1]/2

# test inverse
det1, adj1 = inv(pm1)
t1 = adj1*pm1
@test norm(t1[2,1])/norm(t1[1,1]) < Base.rtoldefault(Float64)
