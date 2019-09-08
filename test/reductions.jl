# triangularization
s = variable("s")
p = PolyMatrix([s-1 s^2-1; 2*one(s) 2s+2; zero(s) 3*one(s)])
L,U = ltriang(p, false, 1)
L₀ = PolyMatrix([-1.225s+1.225 zero(s); -2.450*one(s) zero(s); -1.225*one(s) 1.732*one(s)])
@test isapprox(p*U, L₀; rtol=1e-3)
@test isapprox(L, L₀; rtol=1e-3)
@test isapprox(p*U, L)

R,U = rtriang(copy(transpose(p)), false, 1)
@test isapprox(copy(transpose(R)), L₀; rtol=1e-3)
@test isapprox(U*transpose(p), R)

L,U = ltriang(p)
@test isapprox(p*U, L₀; rtol=1e-3)
@test isapprox(L, L₀; rtol=1e-3)
@test isapprox(p*U, L)

R,U = rtriang(copy(transpose(p)))
@test isapprox(copy(transpose(R)), L₀; rtol=1e-3)
@test isapprox(U*transpose(p), R)

# hermite
s = variable("s")
p = PolyMatrix([-s^3-2s^2+1 -(s+1)^2; (s+2)^2*(s+1) zero(s)])
H₀ = PolyMatrix([s+1 zero(s); (s+2)^2*(s+1) (s+2)^2*(s+1)^2])
U₀ = PolyMatrix([one(s) s+1; -s -s^2-s+1])

H,U = hermite(p)
@test isapprox(p*U, H₀)
@test isapprox(H, H₀)
@test isapprox(U, U₀)

# gcrd
s   = variable("s")
U   = PolyMatrix([zero(s) zero(s) one(s); zero(s) one(s) zero(s); one(s) s+1 -s^2])
R₀  = PolyMatrix([one(s) s+1; zero(s) s^2])
p₁  = PolyMatrix([s^2 zero(s); zero(s) s^2])
p₂  = PolyMatrix([one(s) s+1])

#R, V₁, V₂ = gcrd(p₁, p₂) TODO

#@test hermite(R)[1] ≈ hermite(R₀)[1] TODO
#@test V₁*R ≈ p₁
#@test V₂*R ≈ p₂

#L, V₁, V₂ = gcld(copy(transpose(p₁)),copy(transpose(p₂))) TODO
#@test hermite(L)[1] ≈ hermite(copy(transpose(R₀)))[1]
#@test L*V₁ ≈ copy(transpose(p₁))
#@test L*V₂ ≈ copy(transpose(p₂))

# colred
s = variable("s")

# example 1 from "A Fortran 77 package for column reduction of polynomial matrices" Geurts, A.J. Praagman, C., 1998
p = PolyMatrix([s^4+6s^3+13s^2+12s+4 -s^3-4s^2-5s-2; zero(s) s+2])
p2 = PolyMatrix([s+4 -s; zero(s) s+2])
U₀ = PolyMatrix([one(s) zero(s); s+2 one(s)])
R₀ = PolyMatrix([zero(s) -(s^3+4s^2+5s+2); s^2+4s+4 s+2])
R,U = colred(p)
@test isapprox(R, R₀)
@test isapprox(U, U₀)
@test isapprox(p*U, R)

R1,R2 = colred(p, p2)
@test isapprox(R, R1)

R,U = rowred(copy(transpose(p)))
@test isapprox(R, copy(transpose(R₀)))
@test isapprox(U, copy(transpose(U₀)))
@test isapprox(U*transpose(p), R)

R1,R2 = rowred(copy(transpose(p)), copy(transpose(p2)))
@test isapprox(R, R1)

# example 2 from "A Fortran 77 package for column reduction of polynomial matrices" Geurts, A.J. Praagman, C., 1998
p = PolyMatrix([s^4 s^2 s^6+1; s^2 one(s) s^4; one(s) zero(s) one(s)])
U₀ = PolyMatrix([one(s) zero(s); s+2 one(s)])
R₀ = PolyMatrix([zero(s) -(s^3+4s^2+5s+2); s^2+4s+4 s+2])
#R,U = colred(p)       # TODO this should work ?!
#@test isapprox(p*U, R)
#@test degree(R) == 0

#R,U = rowred(transpose(p))
#@test isapprox(PolyMatrix(U*p.'), R)

# example 3 from "A Fortran 77 package for column reduction of polynomial matrices" Geurts, A.J. Praagman, C., 1998
ϵ = 0.001
p = PolyMatrix([s^3+s^2 ϵ*s+1 one(s); 2s^2 -one(s) -one(s); 3s^2 one(s) one(s)])
R,U = colred(p)
@test isapprox(p*U, R)

#R,U = rowred(copy(transpose(p))) # TODO does not exit
#@test isapprox(U*transpose(p), R)

# example 4 from "A Fortran 77 package for column reduction of polynomial matrices" Geurts, A.J. Praagman, C., 1998
ϵ = exp(-8)
p = PolyMatrix([s^3+s^2+2s+1 ϵ*s^2+2s+3 s^2+s+1   s-1;
                s-1          -s+2       2s^2+s-1  2s+1;
                s+3          2s-1       -s^2-2s+1 -s-2;
                one(s)       -one(s)    3s+1       3*one(s)])
#R,U = colred(p) # TODO this should work ?! does not exit
# @test isapprox(p*U, R)

#R,U = rowred(copy(transpose(p))) # TODO this should work ?! does not exit
#@test isapprox(PolyMatrix(U*p.'), R)

# left2right matrix fractional descriptions
s = variable("s")
Nₗ = PolyMatrix([s^2 zero(s); -4s s])
Dₗ = PolyMatrix([s^3+2s^2-1 s+1; -5s^2-13s-8 (s+1)*(s+4)])
Nᵣ = PolyMatrix([-s^2 -s; zero(s) -s])
Dᵣ = PolyMatrix([-s^3-2s^2+1 -(s+1)^2; (s+2)^2*(s+1) zero(s)])

rmfd = vcat(Dᵣ,Nᵣ)
lmfd = hcat(-Nₗ, Dₗ)

# verify that example is correct.
@test norm(lmfd*rmfd) ≈ 0

#L,U = ltriang(lmfd)

#N = U[3:4,3:4]
#D = U[1:2,3:4]

# compare hermite form
#Dₕ,U = hermite(D)
#Nₕ = N*U

D₀,U = hermite(Dᵣ)
N₀ = Nᵣ*U

#@test isapprox(Dₕ,D₀)
#@test isapprox(Nₕ,N₀)

#@test norm(lmfd*vcat(Dₕ,Nₕ)) < 1e-12

rmfd2 = vcat(-Dᵣ,Nᵣ)
#R,U = rtriang(rmfd2,false)

#N = U[3:4,1:2]
#D = U[3:4,3:4]
#@test isapprox(N*Dᵣ,D*Nᵣ)

# try to get back rfd from obtained lfd
#lmfd2 = hcat(-N,D)

#@test norm(lmfd2*rmfd) < 1e-14
#L,U = ltriang(lmfd2)

#N = U[3:4,3:4]
#D = U[1:2,3:4]

# compare hermite form
#Dₕ,U = hermite(D)
#Nₕ = N*U

D₀,U = hermite(Dᵣ)
N₀ = Nᵣ*U

#@test isapprox(Dₕ,D₀)
#@test isapprox(Nₕ,N₀)
