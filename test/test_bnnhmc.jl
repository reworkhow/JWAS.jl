using Revise

using CSV,DataFrames,JWAS,JWAS.Datasets

phenofile       = Datasets.dataset("example","phenotypes.txt")
genofile   = Datasets.dataset("example","genotypes.txt")
phenotypes = CSV.read(phenofile,delim = ',',header=true,missingstrings=["NA"])


nNodes=[2,3,4]
L=3

geno = get_genotypes(genofile,method="BayesC",estimatePi=true,quality_control=true)

model_equations = "y1 = intercept + geno";
#Neural Network
model = build_model(model_equations,num_latent_traits = nNodes[1],L=L,nNodes=nNodes,nonlinear_function="Neural Network")

out_nn  = runMCMC(model,phenotypes,mega_trait=true,chain_length=10,printout_model_info=false);
