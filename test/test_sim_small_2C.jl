using Revise;
using JWAS,DataFrames,CSV,Statistics,JWAS.Datasets,DelimitedFiles,HTTP;
path = "/Users/apple/Library/CloudStorage/Box-Box/Jiayi_Computer/UCD_PhD/FAANG prediction/Codes/SimData/"
pheno_file = path * "pheno_sim_h2=0.5_small_annoMat_rmepi.csv"
phenotypes = CSV.read(pheno_file,DataFrame,delim = ',',header=true,missingstrings=["NA"])

geno1 = CSV.read(path * "geno1_small.csv",DataFrame,delim = ',',header=true,missingstrings=["NA"])
geno2 = CSV.read(path * "geno2_small.csv",DataFrame,delim = ',',header=true,missingstrings=["NA"])

rowID1 = string.(geno1[:,:ID])
M1 =  Matrix{Float64}(geno1[:,2:end])
header1 = permutedims(names(geno1))

rowID2 = string.(geno2[:,:ID])
M2 =  Matrix{Float64}(geno2[:,2:end])
header2 = permutedims(names(geno2))

genotype1 = get_genotypes(M1,header=header1,rowID=rowID1,method="BayesC", quality_control = false);
genotype2 = get_genotypes(M2,header=header2,rowID=rowID2,method="BayesC", quality_control = false);

model_equation  ="y = intercept + genotype1 + genotype2"
model = build_model(model_equation);
out=runMCMC(model,phenotypes,chain_length = 10000,seed = 123, burnin = 2000, output_folder = "MultiClass_smallGeno_sim2_rmepi");

out["pi_genotype1"]
out["pi_genotype2"]
