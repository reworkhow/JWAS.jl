include("../src/JWAS.jl")

using XSim
using JWAS
using Distributions

k = 50
h2 = 0.5;

ns    = 10
nd    = 1
no    = 40
n     = ns+ns*nd+ns*nd*no
np    = ns*nd
ind   = 1:n
sires = vcat(fill(0,(ns+ns*nd)),repeat(1:ns,inner=nd*no))
dams  = vcat(fill(0,(ns+ns*nd)),repeat(((ns+1):(ns+ns*nd)),inner=no))
ped   = hcat(ind,sires,dams);

numChr,numLoci,chrLength,mutRate = 10,k,10.0,0.0
mapPos     = linspace(chrLength/(numLoci+1),chrLength-chrLength/(numLoci+1),numLoci)
mapPos     = collect(mapPos)
geneFreq   = fill(0.5,numLoci)

XSim.init(numChr,numLoci,chrLength,geneFreq,mapPos,mutRate);

pedArray=XSim.mkPedArray(ped);
animals = samplePed(pedArray);

M = getOurGenotypes(animals);
M = M[:,vec(var(M,1)).!=0.0];
M = float(M);

colmean = mean(M,1)
p       = vec(mean(M,1))/2;
twopq   = 2*p.*(1-p);

k     = size(M,2)
alpha = randn(k)
a     = M*alpha
scalea= sqrt(var(a))
alpha = alpha/scalea
a     = M*alpha
vara  = var(a)
Ve    = vara*(1-h2)/h2;
println("genetic variance is ",vara)
println("residual variance is ",Ve)

sigma_alpha = mean(alpha'alpha)
y           = a + randn(n)*sqrt(Ve);

using DataFrames
data=DataFrame(y=y);


model_equation    = "y = intercept"
residual_variance = Ve
model             = build_model(model_equation,residual_variance)
add_markers(model,M,vara,center=false);

@time output=runMCMC(model,data,chain_length=5000,methods="GBLUP",printout_frequency=1000);
