
include("../src/JWAS.jl")
using DataFrames,JWAS,JWAS.Datasets


phenofile = Datasets.dataset("testMME","simple.txt")
genofile  = Datasets.dataset("testMME","genotype.txt")
pedfile   = Datasets.dataset("testMME","pedigree.txt");
phenotype = readtable(phenofile,separator = ',',header=true);
pedigree = get_pedigree(pedfile);
residual_variance = 1.0
genetic_variance  = 2.5
genetic_variance_by_marker    = 1.5
genetic_variance_by_polygenic = 1.0;
model = build_model("y = intercept + Age + Animal",residual_variance)
set_covariate(model,"Age")
set_random(model,"Animal",pedigree,genetic_variance_by_polygenic)
add_markers(model,genofile,genetic_variance_by_marker,separator=',');

output=runMCMC(model,phenotype,chain_length=5000,methods="GBLUP",printout_frequency=100);
