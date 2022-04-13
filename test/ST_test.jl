# Step 1: Load packages
using JWAS,DataFrames,CSV,Statistics,JWAS.Datasets

# Step 2: Read data
phenofile  = dataset("phenotypes.csv")
pedfile    = dataset("pedigree.csv")
genofile   = dataset("genotypes.csv")
phenotypes = CSV.read(phenofile,DataFrame,delim = ',',header=true,missingstrings=["NA"])
pedigree   = get_pedigree(pedfile,separator=",",header=true);
genotypes  = get_genotypes(genofile,separator=',',method="BayesC");
first(phenotypes,5)

# Step 3: Build Model Equations

model_equation  ="y1 = intercept + x1 + x2 + x2*x3 + ID + dam + genotypes"
model = build_model(model_equation);

# Step 4: Set Factors or Covariates
set_covariate(model,"x1");

# Step 5: Set Random or Fixed Effects
set_random(model,"x2");
set_random(model,"ID dam",pedigree);

# Step 6: Run Analysis
out=runMCMC(model,phenotypes);

# Step 7: Check Accuruacy
results    = innerjoin(out["EBV_y1"], phenotypes, on = :ID) ;
accuruacy  = cor(results[!,:EBV],results[!,:bv1])
