#include("../src/JWAS.jl")

using DataFrames,JWAS,JWAS.Datasets

phenofile = Datasets.dataset("testMME","simple.txt")
genofile  = Datasets.dataset("testMME","genotype.txt")
pedfile   = Datasets.dataset("testMME","pedigree.txt");
phenotype = readtable(phenofile,separator = ',',header=true);
pedigree = get_pedigree(pedfile);
residual_variance = 1.0
genetic_variance  = 2.5
model = build_model("y = intercept + Age + Animal",residual_variance)
set_covariate(model,"Age")
set_random(model,"Animal",pedigree,genetic_variance)


d1 = adjust_phenotypes(model,phenotype)


add_markers(model,genofile,genetic_variance,separator=',',center=false);

pe=LOOCV(model,d1,genetic_variance,residual_variance)
