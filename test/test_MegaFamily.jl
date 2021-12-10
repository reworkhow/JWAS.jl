# Step 1: Load packages
using Revise;
using JWAS,DataFrames,CSV,Statistics,JWAS.Datasets;

# Step 2: Read data
phenofile  = dataset("phenotypes.csv")
genofile   = dataset("genotypes.csv")
phenotypes = CSV.read(phenofile,DataFrame,delim = ',',header=true,missingstrings=["NA"])

genotypes  = get_genotypes(genofile,[1 0.5;0.5 1],separator=',',method="BayesC", G_is_marker_variance = true);

# Step 3: Build Model Equations
model_equation  ="y1 = intercept + genotypes"
model = build_model(model_equation,
		    num_hidden_nodes=3,
		    nonlinear_function="tanh");

model = build_model(model_equation, K = 2, Λ=[1 2 3;2 3 4]);
# Step 4: Run Analysis
out=runMCMC(model,phenotypes,chain_length = 10,
Ymatrix_megaFamily = DataFrame(randn(100,3), :auto));

model.yobs
model.nModels
model.K
model.X
model.ySparse
model.MCMCinfo.mega_trait
model.R
model.sol
model.M[1].G
