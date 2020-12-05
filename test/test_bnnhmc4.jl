# the same file used on server
# using Pkg
# Pkg.add(PackageSpec(url="https://github.com/reworkhow/JWAS.jl.git", rev="TJZ"))

using Revise
using DataFrames,CSV,Random,JWAS,DelimitedFiles,Statistics,Plots

data=2
runseed=1

nNodes=[3,2]
L=2
chainLength = 5_000


############ READ DATA ##########
#read data
phenofile  = "/Users/tianjing/Box/BNN/simulated_data/my_y.n1500.p100.nNodes3.seed$data.txt"
genofile   = "/Users/tianjing/Box/BNN/simulated_data/my_x.n1500.p100.nNodes3.seed$data.txt"
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



######### save accuracy
nNodes_str = join(string.(nNodes))

#hmc
open("/Users/tianjing/Box/BNN/HMC/HMC_test_1205/accuracy_hmc.L$L.nNodes$nNodes_str.nIter$chainLength.runseed$runseed.txt", "w") do io
    writedlm(io, accuracy_hmc)
end

myfig=plot(collect(5:5:5_000),accuracy_hmc,label="hmc",xlabel="iteration",ylabel="accuracy",legend=:bottomright)

savefig(myfig,"/Users/tianjing/Box/BNN/HMC/HMC_test_1205/accuracy.L$L.nNodes$nNodes_str.runseed$runseed.png")
