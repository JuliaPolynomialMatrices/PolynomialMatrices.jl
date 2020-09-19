# test equality and isapprox
p1  = Polynomial([2,1,3.], :s)
p2  = Polynomial([2,1,3.1], :s)
p3  = Polynomial([Inf,1,3.1], :s)
p4  = Polynomial([Inf,1.,3], :s)
p5  = Polynomial([NaN64,1.,3], :s)
m1  = [p1 p2; p2 p1]
m2  = [p1 p2; p2 p2]
m3  = [p3 p2; p2 p3]
m4  = [p4 p2; p2 p4]
m5  = [p5 p2; p2 p4]
pm1 = PolyMatrix(m1)
pm2 = PolyMatrix(m2)
pm3 = PolyMatrix(m3)
pm4 = PolyMatrix(m4)
pm5 = PolyMatrix(m5)
pm6 = PolyMatrix(Matrix{Float64}(I,2,2), :s)
pm7 = PolyMatrix(Matrix{Float64}(I,2,2), :z)

@test pm1 != pm2 != pm3 != pm4 != pm5
@test pm6 != pm7

@test !isapprox(pm2,pm3)
@test !isapprox(pm3,pm4; rtol=0.01)
@test isapprox(pm3,pm4; rtol=0.1)
@test !isapprox(pm3,pm4; rtol=0.001)
@test isapprox(Matrix{Float64}(I,2,2),pm6)
@test !isapprox(pm1,pm3)

@test !isapprox(pm3,m4; rtol=0.01)
@test !isapprox(m3,pm4; rtol=0.01)
@test isapprox(pm3,m4; rtol=0.1)
@test isapprox(m3,pm4; rtol=0.1)

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
C = [0. 1; 1 1]
PolyB = PolyMatrix(B, (2,2), :s)
PolyC = PolyMatrix(B, (2,2), :q)
PolyD = PolyMatrix(C, (2,2), :s)

@test variable(PolyB) == Polynomial([zero(eltype(B)), one(eltype(B))], vartype(PolyB))
@test variable(PolyC) == Polynomial([zero(eltype(B)), one(eltype(B))], vartype(PolyC))
@test variable(PolyD) == Polynomial([zero(eltype(C)), one(eltype(C))], vartype(PolyD))

@test PolyB ≠ PolyC && !isequal(PolyB, PolyC)
@test_throws DomainError PolyB ≈ PolyC
@test PolyB == PolyD && PolyB ≈ PolyD && isequal(PolyB, PolyD)

# test copy
p1  = Polynomial([1])
p2  = Polynomial([2,1,3])
pm1 = PolyMatrix([p1 p2; p2 p1])
pm2 = copy(pm1)

@test pm1 == pm2 && pm1 ≢ pm2

# test similar
@test similar(pm1)                      == PolyMatrix(zeros(Int,2,2))
@test similar(pm1, (2))                 == PolyMatrix(zeros(Int,2))
@test similar(pm1, (2,3))               == PolyMatrix(zeros(Int,2,3))
@test similar(pm1, Int)                 == PolyMatrix(zeros(Int,2,2))

sm1                                     = similar(pm1, Float64, (2,3))
@test sm1                               == PolyMatrix(zeros(Float64,2,3))
@test typeof(sm1)                       == typeof(PolyMatrix(zeros(Float64,2,3)))
@inferred similar(pm1)

@test similar(pm1, Polynomial{Float64})       == PolyMatrix(zeros(Float64,2,2))
@test similar(pm1, Float64)             == PolyMatrix(zeros(Float64,2,2))
@inferred similar(pm1, Polynomial{Float64})

@test similar(pm1, Polynomial{Float64}, (2,)) == PolyMatrix(zeros(Float64,2))
@test similar(pm1, Float64, (2,))       == PolyMatrix(zeros(Float64,2))
@inferred similar(pm1, Polynomial{Float64}, (2,))

# test hcat vcat
p1  = Polynomial([1])
p2  = Polynomial([2,1,3])
m1  = [p1 p2; p2 p1]
m2  = [p2 p1; p2 p2]
pm1 = PolyMatrix(m1)
pm2 = PolyMatrix(m2)

@test hcat(pm1,pm2,pm1) == PolyMatrix(hcat(m1,m2,m1))
@inferred hcat(pm1,pm2,pm1)

@test vcat(pm1,pm2,pm1) == PolyMatrix(vcat(m1,m2,m1))
@inferred vcat(pm1,pm2,pm1)

#@test cat(pm1,pm2,pm1;dims=1) == PolyMatrix(vcat(m1,m2,m1))
#@test cat(pm1,pm2,pm1;dims=2) == PolyMatrix(hcat(m1,m2,m1))
#@inferred cat(pm1,pm2,pm1;dims=1)  # TODO cat is not properly implemented since @inferred do not give correct type

hvcat(1, pm1, pm1)
@inferred hvcat(1, pm1, pm1)

@test [pm1 pm2; pm2 pm1] == PolyMatrix([m1 m2; m2 m1])
@inferred hvcat((2,2), pm1,pm2,pm2,pm1)

#@test vcat(pm1, Matrix{Int}(I,2,2)) == PolyMatrix(vcat(m1, Matrix{Int}(I,2,2)))
#@inferred vcat(pm1, Matrix{Int}(I,2,2))

#@test hcat(pm1, Matrix{Float64}(I,2,2)) == PolyMatrix(hcat(m1, Matrix{Float64}(I,2,2)))
#@inferred vcat(pm1, Matrix{Float64}(I,2,2))

# test getindex
p1  = Polynomial([1],:s)
p2  = Polynomial([2,1,3],:s)
p3  = Polynomial([2,3,4],:s)
m   = [p1 p2; p2 p1]
pm1 = PolyMatrix(m)
degreepm2 = 8
ny  = 2
nu  = 2
A   = randn(ny*(degreepm2+1),nu)
B   = Matrix{Float64}(I,2,2)
pm1[1].coeffs
pm2 = PolyMatrix(A, (ny,nu))
pm3 = PolyMatrix(B)
@test pm1[1].coeffs ≈ p1.coeffs
@test typeof(pm1[1]) == Polynomial{Int}
@inferred pm1[1]

t = pm1[1:2,2]
@test coeffs(t[1]) ≈ coeffs(p2)
@test coeffs(t[2]) ≈ coeffs(p1)
@test typeof(pm1[1:2,2]) == PolyMatrix{Int,Vector{Int},Val{vartype(pm1)},1}
@inferred pm1[1:2,2]

t = pm1[1:2,2:2]
@test typeof(t) == PolyMatrix{Int,Matrix{Int},Val{vartype(pm1)},2}
@inferred pm1[1:2,2:2]

# test setindex!
p1  = Polynomial([0, 1])
p2  = Polynomial([2,1,3])
p3  = Polynomial([1,1,2])
pm1 = PolyMatrix([p1 p2; p2 p1])
@test_throws InexactError pm1[1] = Polynomial([1.5])
pm1[1] = p3
@test coeffs(pm1[1]) ≈ coeffs(p3)

pm1[1,2] = p3
@test coeffs(pm1[1,2]) ≈ coeffs(p3)

pm1[1:2,1:2] = [p3 p3; p3 p3]
@test coeffs(pm1[3]) ≈ coeffs(p3)
@test coeffs(pm1[4]) ≈ coeffs(p3)

# insert numbers
pm1[1] = 1
@test pm1[1] ≈ 1

pm1[1,2] = 1
@test pm1[1,2] ≈ 1

pm1[3:4] = [1 1]
@test pm1[3] ≈ 1

pm2 = PolyMatrix([p1 p1; p1 p1])
pm2[1] = 1
@test pm2[1] ≈ 1

pm2 = PolyMatrix([p1 p1; p1 p1])
pm2[1,2] = 1
@test pm2[1,2] ≈ 1

pm2 = PolyMatrix([p1 p1; p1 p1])
pm2[1:2,1:2] = [1 1; 1 1]
@test pm2[3] ≈ 1

# test insert!
pmins = PolyMatrix([p1 p2; p2 p1])
insert!(pmins, 4, ones(2,2))
@test coeffs(pmins)[4] == ones(2,2)
@test_throws DomainError insert!(pmins, 5, ones(3))

# test iteration
pm2 = PolyMatrix(A, (ny,nu))
for elem in pm2
  @test degree(elem) == degreepm2
end

for idx in eachindex(pm3)
  @test degree(pm3[idx]) <= 0   # Polynomial.jl shifted from degree(p-p)=0 to degree(p-p)=-1
end

@test coeffs(pm3[end]) ≈ coeffs(one(Polynomial{Float64}))
@inferred coeffs(pm3[end])

# test transpose and adjoint
pm4[2] = p3
pm5 = transpose(pm4)
@test coeffs(pm4[2]) ≈ coeffs(pm5[3])
@inferred transpose(pm4)

C   = randn(2,2) + randn(2,2)im
pm6 = PolyMatrix(C)
@test coeffs(adjoint(pm6))[0] ≈ adjoint(C)
@inferred adjoint(pm4)

# test rank
p1  = Polynomial([1],:s)
p2  = Polynomial([2,1,3],:s)
p3  = Polynomial([2,3,4],:s)
m1  = [p1 p2; p2 p1]
m2  = [p1 p2; p1 p2]
pm1 = PolyMatrix(m1)
pm2 = PolyMatrix(m2)

@test rank(pm1) == fastrank(pm1) == 2
@test rank(pm2) == fastrank(pm2) == 1
