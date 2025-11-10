# Step 1: Load packages
cd("/Users/tianjingzhao/Library/CloudStorage/Dropbox/mt_nnmm") # so intermediate files are not saved in the project folder


# ] activate #go to global env
# ] dev /Users/tianjingzhao/Library/CloudStorage/Dropbox/JWAS.jl  #use local JWAS.jl
using Revise #in global env
using JWAS
using Random,DataFrames,CSV,Statistics,JWAS.Datasets, DelimitedFiles
using LinearAlgebra
########################################################
#simulate data
########################################################

# Random.seed!(123)
# genofile   = dataset("genotypes0.csv")
# geno=get_genotypes(genofile)
# M=geno.genotypes
# nind,nsnp=size(M)
# nomics=5
# a=randn(nsnp,nomics)
# omics_g=M*a
# var(omics_g, dims=1)
# omics_g=(omics_g./std(omics_g, dims=1))*sqrt(0.8)
# var(omics_g, dims=1)
# omics = omics_g + randn(nind,nomics).*sqrt(0.2)
# var(omics, dims=1)
# omics_id=["o" * string(i) for i in 1:nomics]
# omics_df=DataFrame(omics,Symbol.(omics_id))
# ind_id=geno.obsID
# insertcols!(omics_df, 1, :ID => ind_id)
# CSV.write("omics.csv",omics_df) 

# phenofile   = dataset("phenotypes.csv")
# phenotypes = CSV.read(phenofile,DataFrame)
# phenotypes=phenotypes[:,[:ID,:y1,:y2]]
# CSV.write("nnmm_phenotypes.csv",phenotypes) 


# musk_matrix=rand([0,1],nsnp,nomics) #nsnp-by-nomics
# writedlm("musk_matrix.txt",musk_matrix,',')


# omics = CSV.read("/Users/tianjingzhao/Library/CloudStorage/Dropbox/mt_nnmm/omics.csv",DataFrame)
# omics.o1 = convert(Vector{Union{Float64, String}}, omics.o1)
# omics[1:3,:o1] .= "NA"
# omics.o2 = convert(Vector{Union{Float64, String}}, omics.o2)
# omics[:,:o2] .= "NA"

# omics.x1 = randn(size(omics,1))
# omics.x2 = randn(size(omics,1))
# omics.x3 = rand(["m","f"],size(omics,1))
# CSV.write("/Users/tianjingzhao/Library/CloudStorage/Dropbox/mt_nnmm/omics_w_na.csv",omics)

# omics = CSV.read("/Users/tianjingzhao/Library/CloudStorage/Dropbox/mt_nnmm/omics_w_na.csv",DataFrame)
# omics.o2 = convert(Vector{Union{Float64, String3}}, omics.o2)
# omics[[1:3;6:100],:o2] = randn(98)
# omics
# CSV.write("/Users/tianjingzhao/Library/CloudStorage/Dropbox/mt_nnmm/omics_w_na2.csv",omics)




# y = CSV.read("/Users/tianjingzhao/Library/CloudStorage/Dropbox/mt_nnmm/nnmm_phenotypes_y1.csv",DataFrame)
# y.x4=randn(size(y,1))
# y.x5=rand([1,2,3],size(y,1))
# y
# CSV.write("/Users/tianjingzhao/Library/CloudStorage/Dropbox/mt_nnmm/nnmm_phenotypes_y1.csv",y)


# y = CSV.read("/Users/tianjingzhao/Library/CloudStorage/Dropbox/mt_nnmm/nnmm_phenotypes_y1.csv",DataFrame)
# #change type of y1 to Union{Missing,Float64}
# y[!,:y1] = convert(Vector{Union{Float64,String}}, y[!,:y1])
# y
# y[99,:y1] = "NA"
# y[100,:y1] = "NA"
# y
# CSV.write("/Users/tianjingzhao/Library/CloudStorage/Dropbox/mt_nnmm/nnmm_phenotypes_y1_w_na.csv",y)

########################################################
#test NNMM
########################################################
g_path   = ["/Users/tianjingzhao/Library/CloudStorage/Dropbox/mt_nnmm/genotypes0.csv"]
# g_path   = repeat(["/Users/tianjingzhao/Library/CloudStorage/Dropbox/mt_nnmm/genotypes0.csv"], 5)

# o_path="/Users/tianjingzhao/Library/CloudStorage/Dropbox/mt_nnmm/omics_w_na.csv"
o_path="/Users/tianjingzhao/Library/CloudStorage/Dropbox/mt_nnmm/omics_w_na2.csv"


# y_path="/Users/tianjingzhao/Library/CloudStorage/Dropbox/mt_nnmm/nnmm_phenotypes_y1.csv"
y_path="/Users/tianjingzhao/Library/CloudStorage/Dropbox/mt_nnmm/nnmm_phenotypes_y1_w_na.csv"

#if polygenic effects are included in the model, the pedigree relationship should be provided
pedfile    = dataset("pedigree.csv")
pedigree = get_pedigree(pedfile,separator=",",header=true);

#other user-defined relationship matrix
GRM = readdlm(dataset("GRM.csv"),',')
grm_names=GRM[:,1]
GRM=Float32.(Matrix(GRM[:,2:end]))+ I*0.0001

G2=0.5
G3=0.1


layers = [
    Layer(layer_name="geno", data_path=g_path), #g_path -> list if partial connected
    Layer(layer_name="omics", data_path=o_path, missing_value="NA"), # should also include other random&covariate data
    Layer(layer_name="phenotypes",data_path=y_path, missing_value="NA") # should also include other random&covariate data
];

equations = [ 
    Equation(
        from_layer_name="geno",
        to_layer_name="omics",
        equation="omics = intercept + x1 + x2 + x3 + ID + geno",
        omics_name = ["o1","o2","o3","o4","o5"],
        covariate=["x1","x2"],
        random=[
            (name="ID", pedigree=pedigree), 
            (name="x3",)],
        method="BayesC"), #nsnp-by-nomics
    Equation(
        from_layer_name="omics",
        to_layer_name ="phenotypes",
        equation ="phenotypes = intercept + ID + x4 + x5 + omics", 
        phenotype_name = ["y1"],
        covariate=["x4"],
        random=[
            (name="x5", G=G3),
            (name="ID", pedigree=pedigree)],
        method="BayesC",
        activation_function="sigmoid")
];

tt=runNNMM(layers, equations; chain_length=10);


