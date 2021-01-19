#this is only to test whether the code can run or not
using Revise

using CSV,DataFrames,JWAS,JWAS.Datasets,Random,Distributions
Random.seed!(123)

phenofile  = Datasets.dataset("example","phenotypes.txt")
genofile   = Datasets.dataset("example","genotypes.txt")
phenotypes = CSV.read(phenofile,delim = ',',header=true,missingstrings=["NA"],DataFrame)
geno = get_genotypes(genofile,method="RR-BLUP",estimatePi=false)
model_equations = "x1 = intercept + geno"


# #MH
# fixed_varz=[0.1 0 0; 0 0.2 0; 0 0 0.3]
# model1 = build_model(model_equations,num_latent_traits = 3,nonlinear_function="Neural Network",fixed_varz=fixed_varz);
# out_nn1  = runMCMC(model1,phenotypes,mega_trait=true,chain_length=5);  #mega_trait can be false


# #HMC
fixed_varz=false#[0.1 0 0; 0 0.2 0; 0 0 0.3]
model2 = build_model(model_equations,num_latent_traits=3,nonlinear_function="Neural Network",
                     fixed_varz=fixed_varz,hmc=true)
out_nn2  = runMCMC(model2,phenotypes,mega_trait=true,chain_length=3)  #mega trait must be true






# #Linear
# model = build_model(model_equations)
# out_st  = runMCMC(model,phenotypes,chain_length=10);
