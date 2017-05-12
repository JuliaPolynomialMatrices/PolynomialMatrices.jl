# test filtering
ny  = 2
p1  = Poly([2,1,3.], :s)
p2  = Poly([2,1,3.1], :s)
pm1 = PolyMatrix([p1 p2; p2 p1])
pm2 = PolyMatrix([p1 p2; p2 p2])
pm3 = PolyMatrix(eye(2,2), :s)

A0 = [2.0 1.0; 0.5 3.0]
coeffs(pm1)[0] = eye(2)   #
coeffs(pm2)[0] = A0       # set first coefficient matrix to an invertible matrix

N = 100
x = randn(ny,N)
out = similar(x)

filt(pm1,pm1,x)
filt(pm2,pm1,x)
filt(pm1,pm2,x)

filt(pm1,pm3,x)
filt(pm3,pm1,x)
PolynomialMatrices._filt_ar!(zeros(x),pm1,x)
