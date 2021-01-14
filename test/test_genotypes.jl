using CSV,DataFrames,JWAS,JWAS.Datasets,DelimitedFiles,CSV

genofile                  = Datasets.dataset("example","genotypes.txt")
genofile_noheader         = Datasets.dataset("example","genotypes_noheader.txt")
geno                      = readdlm(genofile,',')
geno_array,markerIDs,obsID = map(Float64,geno[2:end,2:end]),geno[1,:],geno[2:end,1]
geno_dataframe            = CSV.read(genofile,DataFrame)[:,2:end]

println("load genotype ...")

println("\nfrom a text file with a header (marker IDs)")
model1 = build_model("y = intercept",1.0);
add_genotypes(model1,genofile,1.0,separator=',',header=true,center=true)

println("\nfrom a text file without a header (no marker IDs)")
model2 = build_model("y = intercept",1.0);
add_genotypes(model2,genofile_noheader,1.0,separator=',',header=false,center=true)

println("\nfrom an Array (matrix) with marker IDs and individual IDs provided")
model3 = build_model("y = intercept",1.0);
add_genotypes(model3,geno_array,1.0,header=markerIDs,rowID=obsID,center=true)

println("\nfrom an Array (matrix) without marker IDs and with individual IDs provided")
model4 = build_model("y = intercept",1.0);
add_genotypes(model4,geno_array,1.0,rowID=obsID,center=true)

println("\nfrom an Array (matrix) without marker IDs and individual IDs")
model5 = build_model("y = intercept",1.0);
add_genotypes(model5,geno_array,1.0,center=true)

println("\nfrom a DataFrame with marker IDs and individual IDs provided")
model6 = build_model("y = intercept",1.0);
add_genotypes(model6,geno_dataframe,1.0,header=markerIDs,rowID=obsID,center=true)
