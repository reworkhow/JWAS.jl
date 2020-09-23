using Revise

using DataFrames,CSV,Random,JWAS,DelimitedFiles,Statistics
Random.seed!(123)

nNodes=[2,3,4]
L=3
chainLength = 10_000
burnin=5_000

############ READ DATA ##########
#read phenotype
data_path =  "C:/Users/ztjsw/Box/BNN/data/soy/"
phenofile  = data_path*"soy_pheno.csv"   ##may change
phenotypes = DataFrame!(CSV.File(phenofile))
phenotypes[!,:ID]=string.(phenotypes[!,:ID]);

#read genotype
genofile = data_path*"soy_geno_clean.txt"      ##may change

# read trainID for one specific CV
nObs=size(phenotypes,1)
ID = randperm(nObs)
trainID = ID[1:4000]
testID=ID[4001:end]

############ JWAS ##########
geno = get_genotypes(genofile,method="BayesC",estimatePi=true,quality_control=true);
model_equations = "HT = intercept + geno";

#Neural Network
model = build_model(model_equations,num_latent_traits = nNodes[1],L=L,nNodes=nNodes,nonlinear_function="Neural Network")
out_nn  = runMCMC(model,phenotypes[trainID,:],mega_trait=true,chain_length=chainLength,burnin=burnin);

out_nn_ebv = out_nn["EBV_NonLinear"]
out_nn_ebv[!,:ID] = string.(out_nn_ebv[!,:ID]);
results_nn = innerjoin(out_nn_ebv, phenotypes, on = :ID);
accuruacy=cor(results_nn[testID,:EBV],results_nn[testID,:HT])
