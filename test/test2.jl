using DataFrames, JWAS
n=500
p=5000
phenotypes = DataFrame(ID=1:n,y=randn(n));
genotypes = rand([0.0,1.0,2.0],n,p);
using JWAS
model_equation1   ="y = intercept";
R      = 1.0
model1 = build_model(model_equation1,R);

G3 =1.0
add_genotypes(model1,genotypes,G3);

@time out1=runMCMC(model1,phenotypes,methods="BayesC",estimatePi=false,Pi=0.87,chain_length=1000);


using Revise
using JWAS, DataFrames, CSV,JWAS.Datasets
phenofile       = Datasets.dataset("example","phenotypes.txt")
genofile   = Datasets.dataset("example","genotypes.txt")
phenotypes = CSV.read(phenofile,delim = ',',header=true,missingstrings=["NA"],DataFrame)

geno = get_genotypes(genofile,method="RR-BLUP",estimatePi=false)#df=4,G_is_marker_variance=false

model_equations = "x1 = intercept + geno";
model           = build_model(model_equations);


out=runMCMC(model,phenotypes,chain_length=5);
