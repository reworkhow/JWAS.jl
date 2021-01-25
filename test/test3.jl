using DataFrames, JWAS, CSV

genofile = "/Users/tianjing/Box/BNN/simulated_data/my_x.n500.p5000.nNodes2.seed3.txt"
phenofile   = "/Users/tianjing/Box/BNN/simulated_data/my_y.n500.p5000.nNodes2.seed3.txt"
phenotypes = CSV.read(phenofile,delim = ',',header=true,missingstrings=["NA"],DataFrame)
phenotypes[!,:ID]=string.(phenotypes[!,:ID]);



geno = get_genotypes(genofile,method="RR-BLUP",estimatePi=false);
model_equations_hmc = "y = intercept + geno";

#run HMC
model_hmc = build_model(model_equations_hmc,num_latent_traits = 20,hmc=true,
                        fixed_varz=false,nonlinear_function="Neural Network")
@time out_hmc  = runMCMC(model_hmc,phenotypes,mega_trait=true,chain_length=1000,output_folder="HMC_res");
