using Revise

using DataFrames,CSV,Random,JWAS,DelimitedFiles,Statistics,Plots


runseed=1
nNodes=[2]
L=1
chainLength = 10_00


############ READ DATA ##########
#read data
phenofile  = "C:/Users/ztjsw/Box/BNN/simulated_data/my_y.n1500.p100.nNodes2.seed1.txt"
genofile   = "C:/Users/ztjsw/Box/BNN/simulated_data/my_x.n1500.p100.nNodes2.seed1.txt"
phenotypes = CSV.read(phenofile,delim = ',',header=true,missingstrings=["NA"])
phenotypes[!,:ID]=string.(phenotypes[!,:ID]);


testID = collect(1:500)
trainID = collect(501:1500)



############ JWAS HMC ##########
Random.seed!(runseed)
geno = get_genotypes(genofile,method="BayesC",estimatePi=true,quality_control=true);
model_equations_hmc = "y = intercept + geno";

#run HMC
model_hmc = build_model(model_equations_hmc,num_latent_traits = nNodes[1],L=L,nNodes=nNodes,nonlinear_function="Neural Network")
out_hmc  = runMCMC(model_hmc,phenotypes[trainID,:],mega_trait=true,chain_length=chainLength,output_folder="HMC_res");

#calculate cummulated accuracy
out_hmc_ebv = out_hmc["EBV_NonLinear"]
out_hmc_ebv[!,:ID] = string.(out_hmc_ebv[!,:ID]);

MCMC_ebv_hmc =CSV.read("HMC_res/MCMC_samples_EBV_NonLinear.txt",delim = ',',header=true)
MCMC_ebv_hmc = Matrix(MCMC_ebv_hmc)'  #(p,1000)

cum_ebv_hmc = cumsum(MCMC_ebv_hmc,dims=2)
for i in 1:1000
    cum_ebv_hmc[:,i]=cum_ebv_hmc[:,i]/i
end

cum_ebv_hmc=DataFrame(cum_ebv_hmc) #(p,1000)
insertcols!(cum_ebv_hmc, 1, :ID => out_hmc_ebv[!,:ID])

results_hmc = innerjoin(cum_ebv_hmc, phenotypes, on = :ID)

accuracy_hmc=zeros(1000)
for i in 1:1000
    accuracy_hmc[i] = cor(results_hmc[testID,Symbol("x$i")],results_hmc[testID,:y])
end


# ############ JWAS MH ##########
Random.seed!(runseed)
global geno = get_genotypes(genofile,method="BayesC",estimatePi=true,quality_control=true);
model_equations_mh = "y = intercept + geno";

#run MH
model_mh = build_model(model_equations_mh,num_latent_traits = nNodes[1],nonlinear_function="Neural Network");
out_mh  = runMCMC(model_mh,phenotypes[trainID,:],mega_trait=true,chain_length=chainLength,output_folder="MH_res");

#calculate cummulated accuracy
out_mh_ebv = out_mh["EBV_NonLinear"]
out_mh_ebv[!,:ID] = string.(out_mh_ebv[!,:ID]);

MCMC_ebv_mh =CSV.read("MH_res/MCMC_samples_EBV_NonLinear.txt",delim = ',',header=true)
MCMC_ebv_mh = Matrix(MCMC_ebv_mh)'  #(p,1000)

cum_ebv_mh = cumsum(MCMC_ebv_mh,dims=2)
for i in 1:1000
    cum_ebv_mh[:,i]=cum_ebv_mh[:,i]/i
end

cum_ebv_mh=DataFrame(cum_ebv_mh) #(p,1000)
insertcols!(cum_ebv_mh, 1, :ID => out_mh_ebv[!,:ID])

results_mh = innerjoin(cum_ebv_mh, phenotypes, on = :ID)

accuracy_mh=zeros(1000)
for i in 1:1000
    accuracy_mh[i] = cor(results_mh[testID,Symbol("x$i")],results_mh[testID,:y])
end



############ JWAS Linear ##########
Random.seed!(runseed)
geno = get_genotypes(genofile,method="BayesC",estimatePi=true,quality_control=true);
model_equations_li = "y = intercept + geno";

#liear
model_li = build_model(model_equations_li)
out_li  = runMCMC(model_li,phenotypes[trainID,:],chain_length=chainLength,output_folder="Linear_res");


#calculate cummulated accuracy
out_li_ebv = out_li["EBV_y"]
out_li_ebv[!,:ID] = string.(out_li_ebv[!,:ID]);

MCMC_ebv_li =CSV.read("Linear_res/MCMC_samples_EBV_y.txt",delim = ',',header=true)
MCMC_ebv_li = Matrix(MCMC_ebv_li)'  #(p,1000)

cum_ebv_li = cumsum(MCMC_ebv_li,dims=2)
for i in 1:1000
    cum_ebv_li[:,i]=cum_ebv_li[:,i]/i
end

cum_ebv_li=DataFrame(cum_ebv_li) #(p,1000)
insertcols!(cum_ebv_li, 1, :ID => out_li_ebv[!,:ID])

results_li = innerjoin(cum_ebv_li, phenotypes, on = :ID)

accuracy_li=zeros(1000)
for i in 1:1000
    accuracy_li[i] = cor(results_li[testID,Symbol("x$i")],results_li[testID,:y])
end



######### save accuracy
nNodes_str = join(string.(nNodes))

#hmc
open("accuracy_hmc.L$L.nNodes$nNodes_str.nIter$chainLength.runseed$runseed.txt", "w") do io
    writedlm(io, accuracy_hmc)
end

#mh
open("accuracy_mh.nNodes$nNodes_str.nIter$chainLength.runseed$runseed.txt", "w") do io
    writedlm(io, accuracy_mh)
end

#linear
open("accuracy_li.nIter$chainLength.runseed$runseed.txt", "w") do io
    writedlm(io, accuracy_li)
end


using Plots
myfig=plot(collect(10:10:10_000),randn(1000),label="hmc",xlabel="iteration",ylabel="accuracy",legend=:bottomright)
plot!(collect(10:10:10_000),randn(1000),label="mh")
plot!(collect(10:10:10_000),randn(1000),label="linear")

savefig(myfig,"accuracy.L$L.nNodes$nNodes_str.runseed$runseed.png")
