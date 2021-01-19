
# using Pkg
# Pkg.add(PackageSpec(url="https://github.com/reworkhow/JWAS.jl.git", rev="TJZ"))


using DataFrames,CSV,Random,JWAS,DelimitedFiles,Statistics,Plots

rep=1
num_latent_traits=2
chainLength = 10_000
fixed_varz=false#[1.0 0; 0 1.0]
@show rep,num_latent_traits,chainLength,fixed_varz

############ READ DATA ##########
#read data
phenofile  = "/Users/tianjing/Box/BNN/HMC/test_all3_0117/my_y.n1500.p100.nNodes2.seed1.txt"
genofile   = "/Users/tianjing/Box/BNN/HMC/test_all3_0117/my_x.n1500.p100.nNodes2.seed1.txt"
phenotypes = CSV.read(phenofile,delim = ',',header=true,missingstrings=["NA"],DataFrame)
phenotypes[!,:ID]=string.(phenotypes[!,:ID]);


testID = collect(1:500)
trainID = collect(501:1500)


############ JWAS HMC ##########
Random.seed!(rep)
geno = get_genotypes(genofile,method="RR-BLUP",estimatePi=false);
model_equations_hmc = "y = intercept + geno";

#run HMC
model_hmc = build_model(model_equations_hmc,num_latent_traits = num_latent_traits,hmc=true,
                        fixed_varz=fixed_varz,nonlinear_function="Neural Network")
out_hmc  = runMCMC(model_hmc,phenotypes[trainID,:],mega_trait=true,chain_length=chainLength,output_folder="HMC_res");

#calculate cummulated accuracy
out_hmc_ebv = out_hmc["EBV_NonLinear"]
out_hmc_ebv[!,:ID] = string.(out_hmc_ebv[!,:ID]);

MCMC_ebv_hmc =CSV.read("HMC_res/MCMC_samples_EBV_NonLinear.txt",delim = ',',header=true,DataFrame)
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
@show accuracy_hmc[end]

myfig=plot(collect(chainLength/1000:chainLength/1000:chainLength),accuracy_hmc,label="hmc",xlabel="iteration",ylabel="accuracy",legend=:bottomright)
