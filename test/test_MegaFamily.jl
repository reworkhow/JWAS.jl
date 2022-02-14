# Step 1: Load packages
using Revise;
using JWAS,DataFrames,CSV,Statistics,JWAS.Datasets;

# Step 2: Read data
path = "/Users/apple/Desktop/JiayiQu/UCD_PhD/SparseFactorModel/MegaLMM_data_analysis/Wheat/BayesC_Simulation/SimpleSim"
phenotypes  = CSV.read(path*"/Y_2Factors_4Traits_rep1.csv", DataFrame)
genofile   = path*"/X2_rep1.csv"
genotypes  = get_genotypes(genofile,separator=',',method="BayesC", estimatePi = true);

# Step 3: Build Model Equations
model_equation  ="y = intercept + genotypes"

model = build_model(model_equation, K = 2, Λ = [1 1 0 0; 0 0 1 1], obsTraitNames = ["focal","secondary1","secondary2","secondary3"]);
model.Lamb.is_estimate
# Step 4: Run Analysis
out=runMCMC(model,phenotypes,chain_length =1000);


out
out["heritability"]
#out["pi_Lambda"]

out["EBV_f1"][!,:EBV]
FΛ = [out["EBV_f1"][!,:EBV] out["EBV_f2"][!,:EBV]] * [ 1 1 0 0; 1 0 1 1.]
cor(phenotypes[!,:focal], FΛ[:,1])
