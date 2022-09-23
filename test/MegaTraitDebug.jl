using Revise
# Step 1: Load packages
using JWAS,DataFrames,CSV,Statistics,JWAS.Datasets,Random,DelimitedFiles
Random.seed!(123)
# Step 2: Read data
phenofile  = dataset("phenotypes.csv")
pedfile    = dataset("pedigree.csv")
genofile   = dataset("genotypes.csv")
phenotypes = CSV.read(phenofile,DataFrame,delim = ',',header=true,missingstrings=["NA"])
pedigree   = get_pedigree(pedfile,separator=",",header=true);

genofile, header  = readdlm(genofile,',',header=true)
rowID = genofile[:,1]
genofile =  genofile[:,2:end]
genotypes  = get_genotypes(genofile,header=header,rowID=rowID,method="BayesC");

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
out=runMCMC(model,phenotypes,chain_length=5000,mega_trait=true);

mean(out["pi_genotypes"][out["pi_genotypes"][!,:Trait] .== "y3",:Estimate])

# Step 7: Check Accuruacy
results    = innerjoin(out["EBV_y1"], phenotypes, on = :ID)
accuruacy  = cor(results[!,:EBV],results[!,:bv1])

results    = innerjoin(out["EBV_y2"], phenotypes, on = :ID)
accuruacy  = cor(results[!,:EBV],results[!,:bv2])

results    = innerjoin(out["EBV_y3"], phenotypes, on = :ID)
accuruacy  = cor(results[!,:EBV],results[!,:bv3])





# Step 3: Build Model Equations
model_equation  ="y1 = intercept + genotypes";
model = build_model(model_equation,
                    num_hidden_nodes=3,
                    nonlinear_function="tanh");
# Step 6: Run Analysis
out=runMCMC(model,phenotypes,chain_length=5000,mega_trait=true); #by default mega_trait=true
# Step 7: Check Accuruacy
results    = innerjoin(out["EBV_NonLinear"], phenotypes, on = :ID)
accuruacy  = cor(results[!,:EBV],results[!,:bv1])
