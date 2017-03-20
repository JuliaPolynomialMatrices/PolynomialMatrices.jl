# test constructors
p1  = Poly([1])
p2  = Poly([2,1,3])
p3  = Poly([2,3,4])
m   = [p1 p2; p2 p1]
pm1 = PolyMatrix(m, Val{:x})
#@inferred PolyMatrix(m, Val{:x})

d = Dict(0=>[1 2;2 1], 1=>[0 1;1 0], 2=>[0 3;3 0])
@test pm1 == PolyMatrix(d)
#@inferred PolyMatrix(d, Val{:x})
@test_throws DomainError PolyMatrix(Dict{Int,Matrix{Int}}())

d = Dict(0=>[1 2;2 1], 1=>[0 1;1 0; 0 0])
@test_throws DomainError PolyMatrix(d)

degreepm2 = 8
ny  = 2
nu  = 2
A   = randn(ny*(degreepm2+1),nu)
B   = eye(Float64,2)

@test_throws DomainError PolyMatrix(A, (5,nu))
@test_throws DomainError PolyMatrix(A, (ny,3))
pm2 = PolyMatrix(A, (ny,nu))
pm3 = PolyMatrix(B)
#@inferred PolyMatrix(A, (ny,nu), Val{:x})
#@inferred PolyMatrix(B, Val{:x})

A = randn(3,2,5)
@test PolyMatrix(A, :s) == PolyMatrix(A, Val{:s})
#@inferred PolyMatrix(A, Val{:s})

A = randn(3,2)
@test PolyMatrix(A, :s) == PolyMatrix(A, Val{:s})
#@inferred PolyMatrix(A, Val{:s})

A = randn(3)
@test PolyMatrix(A, :s) == PolyMatrix(A, Val{:s})
#@inferred PolyMatrix(A, Val{:s})

@test pm1(0) == [1 2; 2 1]
