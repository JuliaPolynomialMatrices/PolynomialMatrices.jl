addprocs(3)

using PolynomialMatrices
using Base.Test
using Polynomials
using BenchmarkTools

n = 2
m = 2
p = 2
d = 10
a = PolyMatrix(randn(n,m,d), (n,m,d), :s)
b = PolyMatrix(randn(m,p,d), (m,p,d), :s)

maxn = 4
maxp = 4
maxd = 10

result = SharedArray{Float64}(maxn,maxp,maxd)
for n in 3:maxn
  m = 5n
  for p in 3:maxp
    for d in 9:maxd
      a = PolyMatrix(randn(n,m,d+7), (n,m,d+7), :s)
      b = PolyMatrix(randn(m,p,d-7), (m,p,d-7), :s)
      r1 = @benchmark PolynomialMatrices._mul($a,$b) seconds=2
      r2 = @benchmark PolynomialMatrices._mulfft($a,$b) seconds=2
      result[n,p,d] = median(r1).time/median(r2).time
    end
  end
end

10*10/(10+10)

result

a = PolyMatrix(randn(n,m,d), (n,m,d), :s)
b = PolyMatrix(randn(m,p,d), (m,p,d), :s)

A = zeros(5,5)
B = randn(5,5)

A2 = spzeros(5,5)
B2 = sprandn(5,5,0.1)

@benchmark sumabs2($A) > 0
@benchmark sumabs2($B) > 0

@benchmark _hasnz($A) > 0
@benchmark _hasnz($B) > 0

@benchmark sumabs2($A2) > 0
@benchmark sumabs2($B2) > 0

@benchmark _hasnz($A2) > 0
@benchmark _hasnz($B2) > 0
