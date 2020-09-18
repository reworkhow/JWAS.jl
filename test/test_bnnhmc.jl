using Revise

using CSV,DataFrames,JWAS,JWAS.Datasets

phenofile       = Datasets.dataset("example","phenotypes.txt")
genofile   = Datasets.dataset("example","genotypes.txt")
phenotypes = CSV.read(phenofile,delim = ',',header=true,missingstrings=["NA"])
