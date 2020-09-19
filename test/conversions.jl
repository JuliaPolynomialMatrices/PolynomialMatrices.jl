p1  = Polynomial([2,1,3], :s)
p2  = Polynomial([2,1,3], :s)
p3  = Polynomial([Inf,1,3.1], :s)
p4  = Polynomial([Inf,1.,3], :s)
p5  = Polynomial([NaN64,1.,3], :s)
pm1 = PolyMatrix([p1 p2; p2 p1])
pm2 = PolyMatrix([p1 p2; p2 p2])
pm3 = PolyMatrix([p3 p2; p2 p3])

@test promote_type(typeof(pm1), typeof(pm3)) == typeof(pm3)

pm1n,pm3n = promote(pm1,pm3)
@test typeof(pm1n) == typeof(pm3)
@test typeof(pm3n) == typeof(pm3)
