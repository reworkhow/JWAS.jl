include("../src/JWAS.jl")
using DataFrames,JWAS,JWAS.Datasets

using DataFrames,JWAS,JWAS.Datasets
phenofile = Datasets.dataset("testMT","phenotype.txt")
genofile  = Datasets.dataset("testMT","genotype.txt")
pedfile   = Datasets.dataset("testMT","pedigree.txt");

data = readtable(phenofile,separator = ',',header=true);

R      = [10.0 2.0
           2.0 1.0]

G      = [20.0 1.0
           1.0 2.0]

model_equations = "BW = intercept + age + sex
                   CW = intercept + age + sex";

model1 = build_model(model_equations,R)

set_covariate(model1,"age");

add_markers(model1,genofile,G,separator=',',header=true);
Pi=Dict([1.0; 1.0]=>0.7,[1.0;0.0]=>0.1,[0.0,1.0]=>0.1,[0.0; 0.0]=>0.1)
out = runMCMC(model1,data,Pi=Pi,chain_length=5000,methods="BayesC",
              estimatePi=true,output_samples_frequency=5)
