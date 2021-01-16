#this is only to test whether the code can run or not
using Revise

using CSV,DataFrames,JWAS,JWAS.Datasets,Random,Distributions
Random.seed!(123)

phenofile       = Datasets.dataset("example","phenotypes.txt")
genofile   = Datasets.dataset("example","genotypes.txt")
phenotypes = CSV.read(phenofile,delim = ',',header=true,missingstrings=["NA"],DataFrame)
geno = get_genotypes(genofile,method="RR-BLUP",estimatePi=false)
model_equations = "x1 = intercept + geno"


# #MH
# model = build_model(model_equations,num_latent_traits = 3,nonlinear_function="Neural Network");
# out_nn  = runMCMC(model,phenotypes,mega_trait=true,chain_length=50,printout_model_info=false);


# # #HMC
model3 = build_model(model_equations,num_latent_traits=3,nonlinear_function="Neural Network",
                     hmc=true)
out_nn3  = runMCMC(model3,phenotypes,mega_trait=true,chain_length=3)






# #Linear
# model = build_model(model_equations)
# out_st  = runMCMC(model,phenotypes,chain_length=10);
