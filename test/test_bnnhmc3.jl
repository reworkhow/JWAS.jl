# This file is to simulate data
using Random,Distributions,Statistics,CSV
myseed=3
Random.seed!(myseed)

n=500
p=5000
nNodes=2

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

CSV.write("/Users/tianjing/Box/BNN/simulated_data/my_x.n$n.p$p.nNodes$nNodes.seed$myseed.txt", Z0)

y=DataFrame(ID=collect(1:n),y=y)
CSV.write("/Users/tianjing/Box/BNN/simulated_data/my_y.n$n.p$p.nNodes$nNodes.seed$myseed.txt", y)
