using JWAS
Random.seed!(1)

#read data
phenofile  = "/Users/tianjing/Box/BNN/simulated_data/y.txt"
genofile   = "/Users/tianjing/Box/BNN/simulated_data/x.txt"
phenotypes = CSV.read(phenofile,delim = ',',header=true,missingstrings=["NA"])
phenotypes[!,:ID]=string.(phenotypes[!,:ID]);


testID = collect(101:120)
trainID = collect(1:100)



geno = get_genotypes(genofile,method="BayesC",estimatePi=true,quality_control=true);
model_equations_hmc = "y = intercept + geno";

nNodes=[2]
L=length(nNodes)
chainLength=1000
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
