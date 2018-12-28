include("outputMCMCsamples.jl")
include("DRY.jl")
include("MCMC_BayesC.jl")
include("MCMC_GBLUP.jl")
include("MT_MCMC_BayesC.jl")
include("MT_PBLUP_constvare.jl")
include("../SSBR/SSBR.jl")
include("output.jl")

"""
    runMCMC(mme,df;
    chain_length=1000,starting_value=false,burnin = 0,output_samples_frequency::Int64 = 0,output_file,
    printout_model_info=true,printout_frequency=100,
    methods="conventional (no markers)",Pi=0.0,estimatePi=false,missing_phenotypes=false,
    single_step_analysis= false,pedigree = false,
    constraint=false,update_priors_frequency::Int64=0,
    outputEBV=true)

Run MCMC for Bayesian Linear Mixed Models with estimation of variance components.

* available **methods** include "conventional (no markers)", "RR-BLUP", "BayesB", "BayesC", "Bayesian Lasso", and "GBLUP".
* **starting_value** can be provided as a vector of numbers for all location parameteres and marker effects, defaulting to `0.0`s.
* the first **burnin** iterations are discarded at the beginning of an MCMC run of length **chain_length**
* save MCMC samples every **output_samples_frequency** iterations to files **output_file** defaulting to `MCMC_samples.txt`.
* **Pi** for single-trait analyses is a number; **Pi** for multi-trait analyses is a dictionary such as `Pi=Dict([1.0; 1.0]=>0.7,[1.0; 0.0]=>0.1,[0.0; 1.0]=>0.1,[0.0; 0.0]=>0.1)`,
  defaulting to `all markers have effects on the trait (0.0)` in single-trait analysis and `all markers have effects on all traits` in multi-trait analysis.
* print out the monte carlo mean in REPL with **printout_frequency**
* **constraint**=true if constrain residual covariances between traits to be zeros.
* Individual EBVs are returned if **outputEBV**=true.
"""
function runMCMC(mme::MME,df;
                Pi                  = 0.0,
                burnin              = 0,
                chain_length        = 100,
                starting_value      = false,
                missing_phenotypes  = false,
                constraint          = false,
                estimatePi          = false,
                estimate_variance   = true,
                methods             = "conventional (no markers)",
                output_file         = "MCMC_samples",
                printout_frequency  = chain_length+1,
                printout_model_info = true,
                output_samples_frequency::Int64 = 0,
                update_priors_frequency::Int64=0,
                #parameters for single-step analysis
                single_step_analysis= false,
                pedigree            = false,
                #output
                #Genomic Prediction
                outputEBV               = true,
                output_genetic_variance = false,
                output_residual_variance= false,
                output_heritability     = false)

    ############################################################################
    # Pre-Check
    ############################################################################
    #check errors in function arguments
    errors_args(mme,methods,Pi)
    #users need to provide high-quality pedigree file
    check_pedigree(mme,df,pedigree)
    #user-defined IDs to return genetic values (EBVs)
    check_outputID(outputEBV,mme)
    #check phenotypes, only use phenotypes for individuals in pedigree or genotyped
    check_phenotypes(mme,df,single_step_analysis)
    ############################################################################
    # Single-Step
    ############################################################################
    #impute genotypes for non-genotyped individuals
    #add ϵ (imputation errors) and J as variables in data for non-genotyped individuals
    if single_step_analysis == true
        SSBRrun(mme,pedigree,df)
    end
    ############################################################################
    # Initiate Mixed Model Equations for Non-marker Parts (after SSBRrun for ϵ & J)
    ############################################################################
    starting_value,df = init_mixed_model_equations(mme,df,starting_value)

    if mme.M!=0
        #align genotypes with 1) phenotypes IDs; 2) output IDs.
        align_genotypes(mme)
        Pi = set_marker_hyperparameters_variances_and_pi(mme,Pi,methods)
    end

    if mme.output_ID!=0
        get_outputX_others(mme,single_step_analysis)
    end

    #printout basic MCMC information
    if printout_model_info == true
      getinfo(mme)
      MCMCinfo(methods,Pi,chain_length,burnin,(starting_value!=zeros(size(mme.mmeLhs,1))),printout_frequency,
              output_samples_frequency,missing_phenotypes,constraint,estimatePi,
              update_priors_frequency,mme)
    end

    if mme.nModels ==1
        if methods in ["conventional (no markers)","BayesC","RR-BLUP","BayesB"]
            res=MCMC_BayesC(chain_length,mme,df,
                            burnin                   = burnin,
                            π                        = Pi,
                            methods                  = methods,
                            estimatePi               = estimatePi,
                            sol                      = starting_value,
                            outFreq                  = printout_frequency,
                            output_samples_frequency = output_samples_frequency,
                            output_file              = output_file,
                            update_priors_frequency  = update_priors_frequency)
        elseif methods =="GBLUP"
            res=MCMC_GBLUP(chain_length,mme,df;
                            burnin                   = burnin,
                            sol                      = starting_value,
                            outFreq                  = printout_frequency,
                            output_samples_frequency = output_samples_frequency,
                            output_file              = output_file,
                            update_priors_frequency  = update_priors_frequency)
        else
            error("No options!!!")
        end
    elseif mme.nModels > 1
        if methods == "conventional (no markers)" && estimate_variance == false
          res=MT_MCMC_PBLUP_constvare(chain_length,mme,df,
                            sol    = starting_value,
                            outFreq= printout_frequency,
                            missing_phenotypes=missing_phenotypes,
                            estimate_variance = estimate_variance,
                            output_samples_frequency=output_samples_frequency,
                            output_file=output_file,
                            update_priors_frequency=update_priors_frequency)
        elseif methods in ["BayesL","BayesC","BayesCC","BayesB","RR-BLUP","conventional (no markers)"]
          res=MT_MCMC_BayesC(chain_length,mme,df,
                          Pi     = Pi,
                          sol    = starting_value,
                          outFreq= printout_frequency,
                          missing_phenotypes=missing_phenotypes,
                          constraint = constraint,
                          estimatePi = estimatePi,
                          estimate_variance = estimate_variance,
                          methods    = methods,
                          output_samples_frequency=output_samples_frequency,
                          output_file=output_file,
                          update_priors_frequency=update_priors_frequency)
        else
            error("No methods options!!!")
        end
    else
        error("No options!")
    end
  mme.output = res
  return res
end
