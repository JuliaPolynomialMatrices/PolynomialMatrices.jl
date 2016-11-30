using PolynomialMatrices
using Polynomials
using BenchmarkTools
using DataStructures

ordertype{K,D,O}(::Type{SortedDict{K,D,O}}) = O

for (p,m) in coeffs(Pm)
  println(m)
end

typeof(Pm1)
typeof(Pm2)

promote_type(Pm1, Pm2)
spzeros(size(speye(2))...)

p1 = Poly([1,2,3])
p2 = Poly([1,4,2])
P1 = [p1 p2; p2 p1]
P2 = fill(Poly([1.0,0]),2,2)

Pm1 = PolyMatrix(P1)
Pm2 = PolyMatrix(P2)
c1 = coeffs(Pm1)
c2 = coeffs(Pm2)

R = Pm1 + Pm2

insert!(Pm1,3,[1 2;3 4])

typeof(R)

order(Pm2)


D = SortedDict{Int,typeof(m),ordertype(Pm1.coeffs)}()
insert!(D,0,m)

first(D)

Pm3 = Pm2*Pm1
typeof(Pm3)

m = randn(10,2)

Pm4 = PolyMatrix(m,(5,2))

promote_type(typeof(C), typeof(m))

C = speye(Float64,2)
Pm5 = PolyMatrix(C,(2,2))

start(Pm5.coeffs)

Pm5.coeffs
typeof(Pm5)
k,v = first(Pm5.coeffs)
v
Pm5[1]

typeof(Pm5+Pm2)


promote_type(typeof(Pm5),typeof(Pm4),typeof(Pm3))

show(Pm3)

show(R)
-Pm1

Pm1 - Pm2

Pm3 = Pm1*Pm1
coeffs = SortedDict(Dict{Int,Array{Float64,2}}())
b = typeof(coeffs)
b[3]
methods(b)
eltype(b)

dims = 2,2
var = :x

Base.Order.ForwardOrdering <: Base.Order.Ordering

length(coeffs)

zero(Array{Float64,2}, (2, 2))

PolyMatrix(coeffs,dims,var)

fill(Poly([]),2,2)

similar(dims -> zeros(Float64, dims), indices(R.coeffs[1]))

similar(Array{Float64,2}(), Float64, (2,1))

spzeros(Float64,dims...)

findfirst(R.coeffs)

R(3)

R

R.coeffs[1]
sort(collect(keys(R.coeffs)), rev=true)
