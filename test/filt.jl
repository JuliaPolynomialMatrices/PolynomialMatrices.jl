# test filtering
ny  = 2
p1  = Poly([2,1,3.], :s)
p2  = Poly([2,1,3.1], :s)
pm1 = PolyMatrix([p1 p2; p2 p1])
pm2 = PolyMatrix([p1 p2; p2 p2])
pm3 = PolyMatrix(eye(2,2), :s)

N = 100
x = randn(ny,N)
out = similar(x)

filt(pm1,pm1,x)
filt(pm2,pm1,x)
filt(pm1,pm2,x)

filt(pm1,pm3,x)
filt(pm3,pm1,x)
PolynomialMatrices._filt_ar!(zeros(x),pm1,x)
