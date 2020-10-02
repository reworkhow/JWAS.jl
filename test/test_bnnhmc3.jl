using Random,Distributions,Statistics
Random.seed!(1)

n=1500
p=100
nNodes=3

Z0 = rand([0.0,1.0,2.0],n,p)
W0 = randn(p,nNodes)
Z1 = Z0*W0

pigGrowth(z1i)=sqrt(z1i[1]^2/sum(z1i.^2))  #latent trait

g = [pigGrowth(Z1[i,:]) for i=1:n]
genVar = var(g)
h2=0.9
sigma2e = (1-h2)/h2*genVar
y = g + randn(n)*sqrt(sigma2e)


Z0=DataFrame(Z0)
insertcols!(Z0, 1, :ID => collect(1:n))
Z0

CSV.write("C:/Users/ztjsw/.julia/dev/JWAS/src/5.Datasets/data/example/my_x.txt", Z0)

y=DataFrame(ID=collect(1:n),y=y)
CSV.write("C:/Users/ztjsw/.julia/dev/JWAS/src/5.Datasets/data/example/my_y.txt", y)
