 # Mixed effect neural network: Genotypes -> Unobserved intemediate traits -> Phenotyes
 
 Tips: please center the phenotypes to have zero mean.
 
 ## Example(a): fully-connected neural networks, all intemediate traits are unobserved
- nonlinear function (to define relationship between middle layer and phenotye): tanh (other supported activation functions: "sigmoid", "relu", "leakyrelu", "linear")
- number of nodes in the middle layer: 3
- Bayesian model: multiple independent single-trait BayesC (to sample marker effects on intemediate traits). Note, to use multi-trait Bayesian Alphabet models, please set `mega_trait=false` in `runMCMC()` function.
- sample the unobserved intemediate traits in the middle layer: Hamiltonian Monte Carlo

![](https://github.com/zhaotianjing/figures/blob/main/part2_example.png?raw=true)

```julia
# Step 1: Load packages
using JWAS,DataFrames,CSV,Statistics,JWAS.Datasets,Random
Random.seed!(123)

# Step 2: Read data 
phenofile  = dataset("phenotypes.csv") #get example data path
genofile   = dataset("genotypes.csv")  #get example data path

phenotypes = CSV.read(phenofile,DataFrame,delim = ',',header=true,missingstrings=["NA"]) #read phenotypes (output layer)
genotypes  = get_genotypes(genofile,separator=',',method="BayesC");                      #read genotypes  (input layer)


# Step 3: Build Model Equations 
model_equation  ="y1 = intercept + genotypes"  #name of phenotypes is "y1" in the phenotypes data
                                               #name of genotypes is "genotypes" (user-defined in the previous step)
                                               #the single-trait mixed model used between input and each node in middle layer is: middle node = intercept + genotypes
model = build_model(model_equation,
		    num_hidden_nodes=3,            #number of nodes in middle layer is 3
		    nonlinear_function="tanh");    #tanh function is used to approximate relationship between middle layer and phenotype


# Step 4: Run Analysis
out=runMCMC(model,phenotypes,chain_length=5000); 

# Step 5: Check Accuruacy
results    = innerjoin(out["EBV_NonLinear"], phenotypes, on = :ID) 
accuruacy  = cor(results[!,:EBV],results[!,:bv1])
```


## Example output files

The i-th middle nodes will be named as "trait name"+"i". In our example, the observed trait is named "y1", and there are 3 middle nodes, so the middle nodes are named as "y11", "y12", and "y13", respectively.

Below is a list of files containing estimates and standard deviations for variables of interest. 

| file name      | description |
| -----------    | ----------- |
| EBV_NonLinear.txt | estimated breeding values for observed trait  |
| EBV_y11.txt       | estimated breeding values for middle node 1      |
| EBV_y12.txt       | estimated breeding values for middle node 2      |
| EBV_y13.txt       | estimated breeding values for middle node 3      |
|genetic_variance.txt                   | estimated genetic variance-covariance of all middle nodes |
|heritability.txt                       | estimated heritability of all middle nodes                |
|location_parameters.txt                | estimated bias of all middle nodes                        |
|neural_networks_bias_and_weights.txt.  | estimated bias of phenotypes and weights between middle nodes and phenotypes|
|pi_genotypes.txt                       | estimated pi of all middle nodes                          | 
|marker_effects_genotypes.txt           | estimated marker effects of all middle nodes              |
|residual_variance.txt                  | estimated residual variance-covariance for all middle nodes| 

Below is a list of files containing MCMC samples for variables of interest. 

| file name      | description |
| -----------    | ----------- |
| MCMC_samples_EBV_NonLinear.txt     | MCMC samples from the posterior distribution of breeding values for phenotypes  |
| MCMC_samples_EBV_y11.txt     | MCMC samples from the posterior distribution of breeding values for middle node 1       |
| MCMC_samples_EBV_y12.txt     | MCMC samples from the posterior distribution of breeding values for middle node 2       |
| MCMC_samples_EBV_y13.txt     | MCMC samples from the posterior distribution of breeding values for middle node 3       |
|MCMC_samples_genetic_variance.txt                   | MCMC samples from the posterior distribution of genetic variance-covariance for all middle nodes   |
|MCMC_samples_heritability.txt                       | MCMC samples from the posterior distribution of heritability for all middle nodes                                | 
|MCMC_samples_marker_effects_genotypes_y11            |MCMC samples from the posterior distribution of marker effects for middle node 1           |
|MCMC_samples_marker_effects_genotypes_y12            |MCMC samples from the posterior distribution of marker effects for middle node 2           |
|MCMC_samples_marker_effects_genotypes_y13            |MCMC samples from the posterior distribution of marker effects for middle node 3           |
|MCMC_samples_marker_effects_variances_genotypes.txt | MCMC samples from the posterior distribution of marker effect variance for all middle nodes        |
|MCMC_samples_neural_networks_bias_and_weights.txt.  | MCMC samples from the posterior distribution of bias of observed trait and weights between middle nodes and phenotypes|
|MCMC_samples_pi_genotypes.txt                       | MCMC samples from the posterior distribution of pi for all middle nodes                                          |
|MCMC_samples_residual_variance.txt                  |  MCMC samples from the posterior distribution of residual variance-covariance for all middle nodes |


