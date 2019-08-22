include("outputMCMCsamples.jl")
include("DRY.jl")
include("MCMC_BayesC.jl")
include("MCMC_GBLUP.jl")
include("MT_MCMC_BayesC.jl")
include("MT_PBLUP_constvare.jl")
include("../SSBR/SSBR.jl")
include("output.jl")

"""
    runMCMC(model::MME,df::DataFrame;
            ### MCMC
            chain_length             = 1000,
            starting_value           = false,
            burnin                   = 0,
            output_samples_frequency = 0,
            output_samples_file      = "MCMC_samples",
            update_priors_frequency  = 0,
            ### Methods
            methods                  = "conventional (no markers)",
            estimate_variance        = true,
            Pi                       = 0.0,
            estimatePi               = false,
            estimateScale            = false,
            single_step_analysis     = false,
            pedigree                 = false,
            categorical_trait        = false,
            missing_phenotypes       = false,
            constraint               = false,
            causal_structure         = false,
            ### Genomic Prediction
            outputEBV                = true,
            output_PEV               = false,
            output_heritability      = false,
            ### MISC
            seed                     = false,
            printout_model_info      = true,
            printout_frequency       = 0)

**Run MCMC for Bayesian Linear Mixed Models with or without estimation of variance components.**

* Markov chain Monte Carlo
    * The first **burnin** iterations are discarded at the beginning of a MCMC chain of length **chain_length**.
    * The **starting_value** can be provided as a vector for all location parameteres and marker effects, defaulting to `0.0`s.
    * Save MCMC samples every **output_samples_frequency** iterations, defaulting to `chain_length/1000`, to files **output_samples_file**,
      defaulting to `MCMC_samples.txt`. MCMC samples for hyperparametes (variance componets) and marker effects are saved by default.
      MCMC samples for location parametes can be saved using `output_MCMC_samples()`. Note that saving MCMC samples too frequently slows
      down the computation.
    * Miscellaneous Options
        * Priors are updated every **update_priors_frequency** iterations, defaulting to `0`.
* Methods
    * Available **methods** include "conventional (no markers)", "RR-BLUP", "BayesB", "BayesC", "Bayesian Lasso", and "GBLUP".
    * Single step analysis is allowed if **single_step_analysis** = `true` and **pedigree** is provided.
    * In Bayesian variable selection methods, **Pi** for single-trait analyses is a number; **Pi** for multi-trait analyses is
      a dictionary such as `Pi=Dict([1.0; 1.0]=>0.7,[1.0; 0.0]=>0.1,[0.0; 1.0]=>0.1,[0.0; 0.0]=>0.1)`, defaulting to `all markers
      have effects (Pi = 0.0)` in single-trait analysis and `all markers have effects on all traits
      (Pi=Dict([1.0; 1.0]=>1.0,[0.0; 0.0]=>0.0))` in multi-trait analysis. **Pi** is estimated if **estimatePi** = true
    * Variance components are estimated if **estimate_variance**=true, defaulting to `true`.
    * Scale parameter for prior of marker effect variance is estimated if **estimateScale** = true
    * Miscellaneous Options
        * Missing phenotypes are allowed in multi-trait analysis with **missing_phenotypes**=true, defaulting to `true`.
        * Catogorical Traits are allowed if **categorical_trait**=true, defaulting to false.
        * If **constraint**=true, defaulting to `false`, constrain residual covariances between traits to be zeros.
        * If **causal_structure** is provided, e.g., causal_structure = [0.0,0.0,0.0;1.0,0.0,0.0;1.0,0.0,0.0] for
          trait 2 -> trait 1 and trait 3 -> trait 1, phenotypic causal networks will be incorporated using structure equation models.
* Genomic Prediction
    * Individual estimted breeding values (EBVs) are returned if **outputEBV**=`true`, defaulting to `true`. Heritability and genetic
    variances are returned if **output_heritability**=`true`, defaulting to `false`. Note that estimation of heritability is computaionally intensive.
* Miscellaneous Options
  * Print out the model information in REPL if `printout_model_info=true`; print out the monte carlo mean in REPL with **printout_frequency**,
    defaulting to `false`.
  * If **seed**, defaulting to `false`, is provided, a reproducible sequence of numbers will be generated for random number generation.
"""
function runMCMC(mme::MME,df;
                #MCMC
                chain_length::Int64             = 100,
                starting_value                  = false,
                burnin                          = 0,
                output_samples_frequency::Int64 = chain_length>1000 ? div(chain_length,1000) : 1,
                output_samples_file             = "MCMC_samples",
                update_priors_frequency::Int64  = 0,
                #Methods
                methods                         = "conventional (no markers)",
                estimate_variance               = true,
                Pi                              = 0.0,
                estimatePi                      = false,
                estimateScale                   = false,
                single_step_analysis            = false, #parameters for single-step analysis
                pedigree                        = false, #parameters for single-step analysis
                categorical_trait               = false,
                missing_phenotypes              = true,
                constraint                      = false,
                causal_structure                = false,
                #Genomic Prediction
                outputEBV                       = true,
                output_heritability             = false, #complete or incomplete genomic data
                output_PEV                      = false,
                #MISC
                seed                            = false,
                printout_model_info             = true,
                printout_frequency              = chain_length+1)

    ############################################################################
    # Pre-Check
    ############################################################################
    if seed != false
        Random.seed!(seed)
    end
    if causal_structure != false
        missing_phenotypes = false
        constraint         = true
        #warning 1)be lower triangular
    end
    mme.MCMCinfo = MCMCinfo(chain_length,
                            starting_value,
                            burnin,
                            output_samples_file,
                            output_samples_frequency,
                            printout_model_info,
                            printout_frequency,
                            methods,
                            Pi,
                            estimatePi,
                            estimateScale,
                            single_step_analysis, #pedigree,
                            missing_phenotypes,
                            constraint,
                            estimate_variance,
                            update_priors_frequency,
                            outputEBV,
                            output_heritability,
                            output_PEV,
                            categorical_trait)
    #check errors in function arguments
    errors_args(mme,methods)
    #users need to provide high-quality pedigree file
    check_pedigree(mme,df,pedigree)
    #user-defined IDs to return genetic values (EBVs)
    check_outputID(mme)
    #check phenotypes, only use phenotypes for individuals in pedigree or genotyped
    df = check_phenotypes(mme,df)
    ############################################################################
    # Single-Step
    ############################################################################
    #impute genotypes for non-genotyped individuals
    #add ϵ (imputation errors) and J as variables in data for non-genotyped individuals
    if single_step_analysis == true
        SSBRrun(mme,pedigree,df)
    end
    ############################################################################
    # Initiate Mixed Model Equations for Non-marker Parts (run after SSBRrun for ϵ & J)
    ############################################################################
    starting_value,df = init_mixed_model_equations(mme,df,starting_value)

    if mme.M!=0
        #align genotypes with 1) phenotypes IDs; 2) output IDs.
        align_genotypes(mme,output_heritability,single_step_analysis)
        Pi = set_marker_hyperparameters_variances_and_pi(mme,Pi,methods)
    end

    if mme.output_ID!=0
        get_outputX_others(mme,single_step_analysis)
    end

    #printout basic MCMC information
    if printout_model_info == true
      getinfo(mme)
      getMCMCinfo(methods,Pi,chain_length,burnin,(starting_value!=zeros(size(mme.mmeLhs,1))),printout_frequency,
                  output_samples_frequency,missing_phenotypes,constraint,estimatePi, estimateScale,
                  update_priors_frequency,mme)
    end

    if mme.nModels ==1
        if methods in ["conventional (no markers)","BayesC","RR-BLUP","BayesB"]
            res=MCMC_BayesC(chain_length,mme,df,
                            burnin                   = burnin,
                            π                        = Pi,
                            methods                  = methods,
                            estimatePi               = estimatePi,
                            estimateScale            = estimateScale,
                            sol                      = starting_value,
                            outFreq                  = printout_frequency,
                            output_samples_frequency = output_samples_frequency,
                            output_file              = output_samples_file,
                            update_priors_frequency  = update_priors_frequency,
                            categorical_trait        = categorical_trait)
        elseif methods =="GBLUP"
            res=MCMC_GBLUP(chain_length,mme,df;
                            burnin                   = burnin,
                            sol                      = starting_value,
                            outFreq                  = printout_frequency,
                            output_samples_frequency = output_samples_frequency,
                            output_file              = output_samples_file)
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
                            output_file=output_samples_file,
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
                          output_file=output_samples_file,
                          update_priors_frequency=update_priors_frequency,
                          causal_structure = causal_structure)
        else
            error("No methods options!!!")
        end
    else
        error("No options!")
    end
  mme.output = res
  return res
end
