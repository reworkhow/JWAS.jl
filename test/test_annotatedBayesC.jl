using Revise
using JWAS,DataFrames,CSV,Statistics,JWAS.Datasets,DelimitedFiles;

pheno_file = dataset("pheno_annoMat.csv")
phenotypes = CSV.read(pheno_file,DataFrame,delim = ',',header=true,missingstrings=["NA"])
geno_file  = dataset("geno_annoMat.csv")

nMarker = 100
anno_mat = zeros(nMarker,2)
nMarkerc = Int(nMarker/2)
for i in 1:2
    anno_mat[((i-1)*nMarkerc+ 1):(i*nMarkerc),i] .= 1.
end
anno_mat = [ones(nMarker) anno_mat]
genotypes = get_genotypes(geno_file,separator=',',method="BayesC", quality_control = false, annMat = anno_mat)

model_equation  ="y = intercept + genotypes"
model = build_model(model_equation);
out=runMCMC(model,phenotypes,chain_length = 5000,seed = 123);

results = innerjoin(out["EBV_y"],phenotypes,on=:ID)
accuracy = cor(results[!,:EBV], results[!,:BV])
