# using Revise

using DataFrames,CSV,Random,JWAS,DelimitedFiles,Statistics
Random.seed!(123)

nNodes=[2]
L=1
chainLength = 10
burnin=5

@show nNodes,L,chainLength,burnin

############ READ DATA ##########
#read phenotype
data_path =  "/home/tianjing/BNN_data/soy/"
phenofile  = data_path*"soy_pheno.csv"   ##may change
phenotypes = DataFrame!(CSV.File(phenofile))
phenotypes[!,:ID]=string.(phenotypes[!,:ID]);

#read genotype
genofile = data_path*"soy_geno_clean.txt"      ##may change

# read trainID for one specific CV
trainID = vec(Int64.(readdlm("/home/tianjing/HMC_testing_0923/cv1.trainID.txt")))
testID = vec(Int64.(readdlm("/home/tianjing/HMC_testing_0923/cv1.testID.txt")))



############ JWAS HMC ##########
geno = get_genotypes(genofile,method="BayesC",estimatePi=true,quality_control=true);
model_equations_hmc = "HT = intercept + geno";

#HMC
model_hmc = build_model(model_equations_hmc,num_latent_traits = nNodes[1],L=L,nNodes=nNodes,nonlinear_function="Neural Network")
out_hmc  = runMCMC(model_hmc,phenotypes[trainID,:],mega_trait=true,chain_length=chainLength,burnin=burnin);

#accuracy
out_hmc_ebv = out_hmc["EBV_NonLinear"]
out_hmc_ebv[!,:ID] = string.(out_hmc_ebv[!,:ID]);
results_hmc = innerjoin(out_hmc_ebv, phenotypes, on = :ID);
accuruacy_hmc=cor(results_hmc[testID,:EBV],results_hmc[testID,:HT])


############ JWAS MH ##########
geno = get_genotypes(genofile,method="BayesC",estimatePi=true,quality_control=true);
model_equations_mh = "HT = intercept + geno";

#MH
model_mh = build_model(model_equations_mh,num_latent_traits = nNodes[1],nonlinear_function="Neural Network");
out_mh  = runMCMC(model_mh,phenotypes[trainID,:],mega_trait=true,chain_length=chainLength,burnin=burnin);

#accuracy
out_mh_ebv = out_mh["EBV_NonLinear"]
out_mh_ebv[!,:ID] = string.(out_mh_ebv[!,:ID]);
results_mh = innerjoin(out_mh_ebv, phenotypes, on = :ID);
accuruacy_mh=cor(results_mh[testID,:EBV],results_mh[testID,:HT])


######### save accuracy
accuracy_all = [accuruacy_hmc,accuruacy_mh]
open("accuracy.txt", "w") do io
    writedlm(io, accuracy_all)
end
