#simulate data
using CSV

Random.seed!(1)
nTrain=100
nTest=20
n=nTrain+nTest
l0=10  #nMarker
l1=2   #latent trait

pigGrowth(z1i)=sqrt(z1i[1]^2/sum(z1i.^2))  #latent trait

X  = rand([0.0,1.0,2.0],n,l0)
W0 = randn(l0,l1)
Z1 = X*W0

g = [pigGrowth(Z1[i,:]) for i=1:n]
genVar = var(g)
h2=0.9
sigma2e = (1-h2)/h2*genVar
y = g + randn(n)*sqrt(sigma2e)
@show sigma2e

Z0=DataFrame(X)
insertcols!(Z0, 1, :ID => collect(1:n))

CSV.write("/Users/tianjing/Box/BNN/simulated_data/x.txt", Z0)

y=DataFrame(ID=collect(1:n),y=y)
CSV.write("/Users/tianjing/Box/BNN/simulated_data/y.txt", y)
