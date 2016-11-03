include("../src/JWAS.jl")
using DataFrames,JWAS,JWAS.Datasets

phenofile = Datasets.dataset("testMME","data.txt")
data      = readtable(phenofile,separator = ',',header=true)

model_equation  = "nwn = intercept +parity + parity*site + yr + geneticCode + age"

residual_variance = 2.97
model             = build_model(model_equation,residual_variance)

set_covariate(model,"age")

sow_variance      = 0.26
set_random(model,"parity",sow_variance);

outputMCMCsamples(model,"parity","age");
out=runMCMC(model,data,chain_length=50000,output_samples_frequency=100);

out


pedfile   = Datasets.dataset("testMME","pedigree.txt")
phenofile = Datasets.dataset("testMME","simple.txt");
d1 = readtable(phenofile)
ped = get_pedigree(pedfile);
varRes=1.0
model1 = build_model("y = intercept + Age + Animal",varRes)
set_covariate(model1,"Age")
G=2.5
set_random(model1,"Animal",ped,G)

out = solve(model1,d1,solver="GaussSeidel",printout_frequency=40)
out = runMCMC(model1,d1,chain_length=1000);
