using JWAS,DataFrames,CSV,Statistics,JWAS.Datasets,DelimitedFiles,HTTP;
pheno_file = "/Users/apple/Library/CloudStorage/Box-Box/Jiayi_Computer/UCD_PhD/FAANG prediction/Codes/pheno_sim_h2=0.5.csv"
phenotypes = CSV.read(pheno_file,DataFrame,delim = ',',header=true,missingstrings=["NA"])

function getdata(file_name)  # function to load data from github folder
    http_obj = HTTP.get("https://raw.githubusercontent.com/zhaotianjing/bio_protocol/main/data/$file_name.txt")
    data = CSV.read(http_obj.body, DataFrame,header=true,missingstrings=["."])
end
genotypes = getdata("genotypes_5k");
genotypes

cd("/Users/apple/Library/CloudStorage/Box-Box/Jiayi_Computer/UCD_PhD/FAANG prediction/")

rowID = genotypes[:,:ID]
M =  Matrix(genotypes[:,2:end])
header = permutedims(names(genotypes))

anno_mat = zeros(5000,5)
for i = 1:5
    anno_mat[((i-1)*1000 + 1):(i*1000),i] .= 1.
end
anno_mat

geno  = get_genotypes(M,header=header,rowID=rowID,method="BayesC", quality_control = false, annMat = anno_mat);
geno_wo = get_genotypes(M,header=header,rowID=rowID,method="BayesC", quality_control = false);

model_equation  ="y = intercept + geno"
model = build_model(model_equation);
out=runMCMC(model,phenotypes,chain_length = 20000,seed = 123, burnin = 5000);
out["pi_geno"]


results    = innerjoin(out["EBV_y"], phenotypes, on = :ID)
accuruacy  = cor(results[!,:EBV],results[!,:BV])
