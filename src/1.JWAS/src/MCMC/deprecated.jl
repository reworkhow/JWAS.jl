M   = mme.M[1].genotypes
X   = M
XPX = X'X
lhs = XPX+I*mme.R[1,1]/mme.M[1].G[1,1]
Ch = cholesky(lhs)
iL = inv(Ch.L)
iLhs = inv(Ch)

wArray[1][:] =  wArray[1][:] + X*Mi.α[1]
wArray[2][:] =  wArray[2][:] + X*Mi.α[2]
rhs = X'wArray[1]
sol = iLhs*rhs
p = size(lhs,2)
mysol1 = sol + iL'randn(p)
Mi.α[1] = mysol1

rhs = X'wArray[2]
sol = iLhs*rhs
p = size(lhs,2)
mysol2 = sol + iL'randn(p)
Mi.α[2] = mysol2

wArray[1][:] =  wArray[1][:] - X*mysol1
wArray[2][:] =  wArray[2][:] - X*mysol2
