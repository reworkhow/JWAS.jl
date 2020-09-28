using Revise

using CSV,DataFrames,JWAS,JWAS.Datasets,Random
Random.seed!(123)

phenofile       = Datasets.dataset("example","phenotypes.txt")
genofile   = Datasets.dataset("example","genotypes.txt")
phenotypes = CSV.read(phenofile,delim = ',',header=true,missingstrings=["NA"])




geno = get_genotypes(genofile,method="BayesC")

model_equations = "x1 = intercept + geno"


# #MH
# model = build_model(model_equations,num_latent_traits = nNodes[1],nonlinear_function="Neural Network");
# out_nn  = runMCMC(model,phenotypes,mega_trait=true,chain_length=5,printout_model_info=false);


# # #HMC
nNodes=[2]
L=1
model = build_model(model_equations,num_latent_traits = nNodes[1],L=L,nNodes=nNodes,nonlinear_function="Neural Network")
out_nn  = runMCMC(model,phenotypes,mega_trait=true,chain_length=5)


# #Linear
# model = build_model(model_equations)
# out_st  = runMCMC(model,phenotypes,chain_length=10);
