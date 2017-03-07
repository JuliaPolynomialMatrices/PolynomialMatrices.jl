# test filtering
A   = randn(ny*(degreepm2+1),nu)
B   = eye(Float64,2)
pm1 = PolyMatrix(A, (ny,nu))
pm2 = PolyMatrix(B)

N = 100
x = randn(ny,N)
out = similar(x)

filt(pm1,pm1,x)
filt(pm2,pm1,x)
filt(pm1,pm2,x)
