using Revise

using DataFrames,CSV,Random,JWAS,DelimitedFiles,Statistics,JWAS.Datasets
myseed=123
Random.seed!(myseed)


nNodes=[3,2]
L=2
chainLength = 10000


############ READ DATA ##########
#read data
phenofile  = Datasets.dataset("example","my_y.txt")
genofile   = Datasets.dataset("example","my_x.txt")
phenotypes = CSV.read(phenofile,delim = ',',header=true,missingstrings=["NA"])
phenotypes[!,:ID]=string.(phenotypes[!,:ID]);


trainID = collect(1:1000)
testID = collect(1001:1500)


############ JWAS HMC ##########
global geno = get_genotypes(genofile,method="BayesC",estimatePi=true,quality_control=true);
model_equations_hmc = "y = intercept + geno";

#HMC
model_hmc = build_model(model_equations_hmc,num_latent_traits = nNodes[1],L=L,nNodes=nNodes,nonlinear_function="Neural Network")
out_hmc  = runMCMC(model_hmc,phenotypes[trainID,:],mega_trait=true,chain_length=chainLength);

#accuracy
out_hmc_ebv = out_hmc["EBV_NonLinear"]
out_hmc_ebv[!,:ID] = string.(out_hmc_ebv[!,:ID]);
results_hmc = innerjoin(out_hmc_ebv, phenotypes, on = :ID);
accuruacy_hmc=cor(results_hmc[testID,:EBV],results_hmc[testID,:y])


# ############ JWAS MH ##########
# global geno = get_genotypes(genofile,method="BayesC",estimatePi=true,quality_control=true);
# model_equations_mh = "y = intercept + geno";
#
# #MH
# model_mh = build_model(model_equations_mh,num_latent_traits = nNodes[1],nonlinear_function="Neural Network");
# out_mh  = runMCMC(model_mh,phenotypes[trainID,:],mega_trait=true,chain_length=chainLength);
#
# #accuracy
# out_mh_ebv = out_mh["EBV_NonLinear"]
# out_mh_ebv[!,:ID] = string.(out_mh_ebv[!,:ID]);
# results_mh = innerjoin(out_mh_ebv, phenotypes, on = :ID);
# accuruacy_mh=cor(results_mh[testID,:EBV],results_mh[testID,:y])


############ JWAS Linear ##########
global geno = get_genotypes(genofile,method="BayesC",estimatePi=true,quality_control=true);
model_equations_li = "y = intercept + geno";

#liear
model_li = build_model(model_equations_li)
out_li  = runMCMC(model_li,phenotypes[trainID,:],chain_length=chainLength);

#accuracy
out_li_ebv = out_li["EBV_y"]
out_li_ebv[!,:ID] = string.(out_li_ebv[!,:ID]);
results_li = innerjoin(out_li_ebv, phenotypes, on = :ID);
accuruacy_li=cor(results_li[testID,:EBV],results_li[testID,:y])


######### save accuracy
accuracy_all = [accuruacy_hmc,accuruacy_li]
open("accuracy.txt", "w") do io
    writedlm(io, accuracy_all)
end
