# Step 1: Load packages
using JWAS,DataFrames,CSV,Statistics,JWAS.Datasets,DelimitedFiles

# Step 2: Read data
phenofile  = dataset("phenotypes.csv")
pedfile    = dataset("pedigree.csv")
genofile   = dataset("genotypes.csv")
phenotypes = CSV.read(phenofile,DataFrame,delim = ',',header=true,missingstrings=["NA"])
pedigree   = get_pedigree(pedfile,separator=",",header=true);
genotypes  = get_genotypes(genofile,separator=',',method="BayesC");
first(phenotypes,5)
# Step 3: Build Model Equations
model_equation  ="y1 = intercept + x1 + x2 + x2*x3 + ID + dam + genotypes
                  y2 = intercept + x1 + x2 + ID + genotypes
                  y3 = intercept + x1 + ID + genotypes";
model = build_model(model_equation);
# Step 4: Set Factors or Covariates
set_covariate(model,"x1");
# Step 5: Set Random or Fixed Effects
set_random(model,"x2");
set_random(model,"ID dam",pedigree);
# Step 6: Run Analysis
# If `causal_structure` is provided, e.g., causal_structure = [0.0 0.0 0.0;1.0 0.0 0.0;1.0 0.0 0.0] for
# trait 2 -> trait 1 and trait 3 -> trait 1, phenotypic causal networks will be incorporated using structure equation models.
my_structure = [0.0 0.0 0.0
                1.0 0.0 0.0
                1.0 0.0 0.0]
out=runMCMC(model,phenotypes,causal_structure=my_structure);
