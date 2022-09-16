using Revise
using JWAS,DataFrames,CSV,Statistics,JWAS.Datasets,DelimitedFiles,HTTP;
cd("/Users/apple/Library/CloudStorage/Box-Box/Jiayi_Computer/UCD_PhD/FAANG prediction/")
path = "/Users/apple/Library/CloudStorage/Box-Box/Jiayi_Computer/UCD_PhD/FAANG prediction/Codes/SimData/"
pheno_file = path * "pheno_sim_h2=0.5_small_annoMat_rmepi.csv"
phenotypes = CSV.read(pheno_file,DataFrame,delim = ',',header=true,missingstrings=["NA"])

geno = CSV.read(path * "geno_small.csv",DataFrame,delim = ',',header=true,missingstrings=["NA"])
rowID = string.(geno[:,:ID])
M =  Matrix{Float64}(geno[:,2:end])
header = permutedims(names(geno))

nMarker = 100
anno_mat = zeros(nMarker,2)
nMarkerc = Int(nMarker/2)
for i in 1:2
    anno_mat[((i-1)*nMarkerc+ 1):(i*nMarkerc),i] .= 1.
end
anno_mat = [ones(nMarker) anno_mat]
anno_mat
#=
anno_mat[1:66,1] .= 1.
anno_mat[33:100,2] .= 1.
anno_mat
=#

genotype = get_genotypes(M,header=header,rowID=rowID,method="BayesC", quality_control = false, annMat = anno_mat);

model_equation  ="y = intercept + genotype"
model = build_model(model_equation);
out=runMCMC(model,phenotypes,chain_length = 10000,seed = 123, burnin = 2000, output_folder = "AnnMat_smallGeno_sim2_adj_pi_10K");

out["pi_genotype"]

QTL_info = CSV.read(path * "QTL_info_small_annoMat_rmepi.csv", DataFrame)

mean(out["pi_genotype"][1:50,:Estimate])
mean(out["pi_genotype"][51:100,:Estimate])

model.M[1].annCoef

QTL_info_cat1 = QTL_info[!,:position][1:19]
QTL_info_cat12 = QTL_info[!,:position][20:34]
QTL_info_cat2 = QTL_info[!,:position][35:47]
SNP_names = names(geno)[2:end]

QTL_index = collect(1:nMarker)[in(string.(QTL_info[!,:position])).(SNP_names)]
mean(out["pi_genotype"][QTL_index,:Estimate])
nonQTL_index = deleteat!(collect(1:nMarker), QTL_index)
mean(out["pi_genotype"][nonQTL_index,:Estimate])

QTL_index1 = collect(1:nMarker)[in(string.(QTL_info_cat1)).(SNP_names)]
mean(out["pi_genotype"][QTL_index1,:Estimate])

QTL_index12 = collect(1:nMarker)[in(string.(QTL_info_cat12)).(SNP_names)]
mean(out["pi_genotype"][QTL_index12,:Estimate])

QTL_index2 = collect(1:nMarker)[in(string.(QTL_info_cat2)).(SNP_names)]
mean(out["pi_genotype"][QTL_index2,:Estimate])
