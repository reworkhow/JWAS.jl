using Revise
using JWAS,DataFrames,CSV,Statistics,JWAS.Datasets,DelimitedFiles;

phenofile  = dataset("phenotypes.csv")
pedfile    = dataset("pedigree.csv")
genofile   = dataset("genotypes.csv")


phenotypes = CSV.read(phenofile,DataFrame,delim = ',',header=true,missingstrings=["NA"])
pedigree   = get_pedigree(pedfile,separator=",",header=true);

genofile2, header  = readdlm(genofile,',',header=true)
rowID = genofile2[:,1]
genofile2 =  genofile2[:,2:end]
genotypes  = get_genotypes(genofile2,header=header,rowID=rowID,method="BayesC", quality_control = false, annMat = [ones(1000) ones(1000)]);
#annMat = [ones(1000) ones(1000)]
# Pi = 0.5
genotypes.annMat

model_equation  ="y1 = intercept + x1 + x2 + x2*x3 + ID + dam + genotypes"
model = build_model(model_equation);

set_covariate(model,"x1");
set_random(model,"x2");
set_random(model,"ID dam", pedigree);

out=runMCMC(model,phenotypes,chain_length = 5000,seed = 123);
out["pi_genotypes"]

results    = innerjoin(out["EBV_y1"], phenotypes, on = :ID)
accuruacy  = cor(results[!,:EBV],results[!,:bv1])

out=runMCMC(model,phenotypes,chain_length = 5000,seed = 234);
out["pi_genotypes"]
results    = innerjoin(out["EBV_y1"], phenotypes, on = :ID)
accuruacy  = cor(results[!,:EBV],results[!,:bv1])
