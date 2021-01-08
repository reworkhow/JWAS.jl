using Distributions,Printf,Random
using DelimitedFiles
using InteractiveUtils #for versioninfo
using DataFrames,CSV
using SparseArrays
using LinearAlgebra
using ProgressMeter
using .PedModule

import StatsBase: describe #a new describe is exported

#Models
include("types.jl")
include("build_MME.jl")
include("random_effects.jl")
include("residual.jl")
include("variance_components.jl")

#Iterative Solver
include("iterative_solver/solver.jl")

#Markov chain Monte Carlo
include("MCMC/MCMC_BayesianAlphabet.jl")

#Genomic Markers
include("markers/tools4genotypes.jl")
include("markers/readgenotypes.jl")
include("markers/BayesianAlphabet/BayesABC.jl")
include("markers/BayesianAlphabet/BayesC0L.jl")
include("markers/BayesianAlphabet/GBLUP.jl")
include("markers/BayesianAlphabet/MTBayesABC.jl")
include("markers/BayesianAlphabet/MTBayesC0L.jl")
include("markers/Pi.jl")

#Incomplete Genomic Data (Single-step Methods)
include("single_step/SSBR.jl")
include("single_step/SSGBLUP.jl")

#Categorical trait
include("categorical_trait/categorical_trait.jl")
include("censored_trait/censored_trait.jl")

#Structure Equation Models
include("structure_equation_model/SEM.jl")

#Latent Traits
include("Nonlinear/nonlinear.jl")

#input
include("input_data_validation.jl")

#output
include("output.jl")

export build_model,set_covariate,set_random,add_genotypes,get_genotypes
export outputMCMCsamples,outputEBV,getEBV
export solve,runMCMC
export showMME,describe
#Pedmodule
export get_pedigree,get_info
#misc
export GWAS,QC
export get_additive_genetic_variances,get_breeding_values
export get_correlations,get_heritability
#others
#export adjust_phenotypes,LOOCV

"""
    runMCMC(model::MME,df::DataFrame;
            weight                   = false,
            ### MCMC
            chain_length             = 1000,
            starting_value           = false,
            burnin                   = 0,
            output_samples_frequency = chain_length/1000,
            output_folder            = results,
            update_priors_frequency  = 0,
            ### Methods
            estimate_variance        = true,
            estimateScale            = false,
            single_step_analysis     = false,
            pedigree                 = false,
            categorical_trait        = false,
            missing_phenotypes       = true,
            constraint               = false,
            causal_structure         = false,
            ### Genomic Prediction
            outputEBV                = true,
            output_heritability      = true,
            ### MISC
            seed                     = false,
            printout_model_info      = true,
            printout_frequency       = 0)

**Run MCMC for Bayesian Linear Mixed Models with or without estimation of variance components.**

* Markov chain Monte Carlo
    * The first `burnin` iterations are discarded at the beginning of a MCMC chain of length `chain_length`.
    * Save MCMC samples every `output_samples_frequency` iterations, defaulting to `chain_length/1000`, to a folder `output_folder`,
      defaulting to `results`. MCMC samples for hyperparametes (variance componets) and marker effects are saved by default.
      MCMC samples for location parametes can be saved using `output_MCMC_samples()`. Note that saving MCMC samples too frequently slows
      down the computation.
    * The `starting_value` can be provided as a vector for all location parameteres and marker effects, defaulting to `0.0`s.
      The order of starting values for location parameters and marker effects should be the order of location parameters in
      the Mixed Model Equation for all traits (This can be obtained by getNames(model)) and then markers for all traits (all
      markers for trait 1 then all markers for trait 2...).
    * Miscellaneous Options
        * Priors are updated every `update_priors_frequency` iterations, defaulting to `0`.
* Methods
    * Single step analysis is allowed if `single_step_analysis` = `true` and `pedigree` is provided.
    * Variance components are estimated if `estimate_variance`=true, defaulting to `true`.
    * Miscellaneous Options
        * Missing phenotypes are allowed in multi-trait analysis with `missing_phenotypes`=true, defaulting to `true`.
        * Catogorical Traits are allowed if `categorical_trait`=true, defaulting to `false`. Phenotypes should be coded as 1,2,3...
        * Censored traits are allowed if the upper bounds are provided in `censored_trait` as an array, and lower bounds are provided as phenotypes.
        * If `constraint`=true, defaulting to `false`, constrain residual covariances between traits to be zeros.
        * If `causal_structure` is provided, e.g., causal_structure = [0.0,0.0,0.0;1.0,0.0,0.0;1.0,0.0,0.0] for
          trait 2 -> trait 1 and trait 3 -> trait 1, phenotypic causal networks will be incorporated using structure equation models.
* Genomic Prediction
    * Predicted values for individuals of interest can be obtained based on an user-defined prediction equation `prediction_equation`, e.g., "y1:animal + y1:geno + y1:age".
    For now, genomic data is always included. Genetic values including effects defined with genotype and pedigre information are returned if `prediction_equation`= false, defaulting to `false`.
    * Individual estimted genetic values and prediction error variances (PEVs) are returned if `outputEBV`=true, defaulting to `true`. Heritability and genetic
    variances are returned if `output_heritability`=`true`, defaulting to `true`. Note that estimation of heritability is computaionally intensive.
* Miscellaneous Options
  * Print out the model information in REPL if `printout_model_info`=true; print out the monte carlo mean in REPL with `printout_frequency`,
    defaulting to `false`.
  * If `seed`, defaulting to `false`, is provided, a reproducible sequence of numbers will be generated for random number generation.
"""
function runMCMC(mme::MME,df;
                #Data
                heterogeneous_residuals           = false,
                #MCMC
                chain_length::Integer             = 100,
                starting_value                    = false,
                burnin::Integer                   = 0,
                output_samples_frequency::Integer = chain_length>1000 ? div(chain_length,1000) : 1,
                update_priors_frequency::Integer  = 0,
                #Methods
                estimate_variance               = true,
                single_step_analysis            = false, #parameters for single-step analysis
                pedigree                        = false, #parameters for single-step analysis
                fitting_J_vector                = true,
                categorical_trait               = false,
                censored_trait                  = false,
                missing_phenotypes              = true,
                constraint                      = false,
                causal_structure                = false,
                mega_trait                      = false,
                #Genomic Prediction
                outputEBV                       = true,
                output_heritability             = true,  #complete or incomplete genomic data
                prediction_equation             = false,
                #MISC
                seed                            = false,
                printout_model_info             = true,
                printout_frequency              = chain_length+1,
                big_memory                      = false,
                double_precision                = false,
                #MCMC samples (defaut to marker effects and hyperparametes (variance componets))
                output_folder                   = "results",
                output_samples_for_all_parameters = false,
                #for deprecated JWAS
                methods                         = "conventional (no markers)",
                Pi                              = 0.0,
                estimatePi                      = false,
                estimateScale                   = false) #calculare lsmeans
    #for deprecated JWAS fucntions
    if mme.M != 0
        for Mi in mme.M
            if Mi.name == false
                Mi.name              = "geno"
                Mi.π                 = Pi
                Mi.estimatePi        = estimatePi
                Mi.estimateScale     = estimateScale
                Mi.method            = methods
            end
        end
    end
    ############################################################################
    # Set a seed in the random number generator
    ############################################################################
    if seed != false
        Random.seed!(seed)
    end
    ############################################################################
    # Create an folder to save outputs
    ############################################################################
    myfolder,folderi = output_folder, 1
    while ispath(output_folder)
        printstyled("The folder $output_folder already exists.\n",bold=false,color=:red)
        output_folder = myfolder*string(folderi)
        folderi += 1
    end
    mkdir(output_folder)
    printstyled("The folder $output_folder is created to save results.\n",bold=false,color=:green)
    ############################################################################
    # Save MCMC argumenets in MCMCinfo
    ############################################################################
    mme.MCMCinfo = MCMCinfo(chain_length,burnin,output_samples_frequency,
                   printout_model_info,printout_frequency, single_step_analysis,
                   fitting_J_vector,missing_phenotypes,constraint,mega_trait,estimate_variance,
                   update_priors_frequency,outputEBV,output_heritability,prediction_equation,
                   categorical_trait,censored_trait,
                   seed,double_precision,output_folder)
    ############################################################################
    # Check 1)Arguments; 2)Input Pedigree,Genotype,Phenotypes,
    #       3)output individual IDs; 4)Priors  5)Prediction equation
    #       in the Whole Dataset.
    ############################################################################
    errors_args(mme)       #check errors in function arguments
    df=check_pedigree_genotypes_phenotypes(mme,df,pedigree)
    prediction_setup(mme)  #set prediction equation, defaulting to genetic values
    check_outputID(mme)    #check individual of interest for prediction
    df,df_whole = make_dataframes(df,mme)
    set_default_priors_for_variance_components(mme,df)  #check priors (set default priors)

    if mme.M!=0
        if single_step_analysis == true
            SSBRrun(mme,df,df_whole,big_memory)
        end
        set_marker_hyperparameters_variances_and_pi(mme)
    end
    ############################################################################
    # Adhoc functions
    ############################################################################
    #save MCMC samples for all parameters (?seperate function user call)
    if output_samples_for_all_parameters == true
        allparameters=[term[2] for term in split.(keys(mme.modelTermDict),":")]
        allparameters=unique(allparameters)
        for parameter in allparameters
            outputMCMCsamples(mme,parameter)
        end
    end
    #structure equation model
    mme.causal_structure = causal_structure
    if causal_structure != false
        #no missing phenotypes and residual covariance for identifiability
        missing_phenotypes, constraint = false, true
        if !istril(causal_structure)
            error("The causal structue needs to be a lower triangular matrix.")
        end
    end
    #Nonlinear
    if mme.latent_traits == true
        yobs = df[!,Symbol(string(Symbol(mme.lhsVec[1]))[1:(end-1)])]
        for i in mme.lhsVec
            df[!,i]= yobs
        end
    end
    # Double Precision
    if double_precision == true
        if mme.M != 0
            for Mi in mme.M
                Mi.genotypes = map(Float64,Mi.genotypes)
                Mi.G         = map(Float64,Mi.G)
                Mi.α         = map(Float64,Mi.α)
            end
        end
        for random_term in mme.rndTrmVec
            random_term.Vinv  = map(Float64,random_term.Vinv)
            random_term.GiOld = map(Float64,random_term.GiOld)
            random_term.GiNew = map(Float64,random_term.GiNew)
            random_term.Gi    = map(Float64,random_term.Gi)
        end
        mme.Gi = map(Float64,mme.Gi)
    end
    #mega_trait
    if mme.MCMCinfo.mega_trait == true || mme.MCMCinfo.constraint == true
        if mme.nModels == 1
            error("more than 1 trait is required for MegaLMM analysis.")
        end
        mme.MCMCinfo.constraint = true
        ##sample from scale-inv-⁠χ2, not InverseWishart
        mme.df.residual  = mme.df.residual - mme.nModels
        mme.scaleR       = diag(mme.scaleR/(mme.df.residual - 1))*(mme.df.residual-2)/mme.df.residual #diag(R_prior_mean)*(ν-2)/ν
        if mme.M != 0
            for Mi in mme.M
                Mi.df        = Mi.df - mme.nModels
                Mi.scale    = diag(Mi.scale/(Mi.df - 1))*(Mi.df-2)/Mi.df
            end
        end
    end
    ############################################################################
    #Make incidence matrices and genotype covariates for training observations
    #and individuals of interest
    ############################################################################
    #make incidence matrices (non-genomic effects) (after SSBRrun for ϵ & J)
    df=make_incidence_matrices(mme,df,df_whole,heterogeneous_residuals)
    #align genotypes with 1) phenotypes IDs; 2) output IDs.
    if mme.M != false
        align_genotypes(mme,output_heritability,single_step_analysis)
    end
    # initiate Mixed Model Equations and check starting values
    init_mixed_model_equations(mme,df,starting_value)

    #printout basic MCMC information
    if printout_model_info == true
    end
    ############################################################################
    # MCMC
    ############################################################################
    describe(mme)
    mme.output=MCMC_BayesianAlphabet(mme,df)

    ############################################################################
    # Save output to text files
    ############################################################################
    for (key,value) in mme.output
      CSV.write(output_folder*"/"*replace(key," "=>"_")*".txt",value)
    end
    if mme.M != 0
        for Mi in mme.M
            if Mi.name == "GBLUP"
                mv("marker_effects_variance_"*Mi.name*".txt","genetic_variance(REML)_"*Mi.name*".txt")
            end
        end
    end
    printstyled("\n\nThe version of Julia and Platform in use:\n\n",bold=true)
    versioninfo()
    printstyled("\n\nThe analysis has finished. Results are saved in the returned ",bold=true)
    printstyled("variable and text files. MCMC samples are saved in text files.\n\n\n",bold=true)
    return mme.output
end
################################################################################
# Print out Model or MCMC information
################################################################################
"""
    describe(model::MME)

* Print out model information.
"""
function describe(model::MME;data=false)
  if model.MCMCinfo == false || model.MCMCinfo.printout_model_info == false
      return "Model information is not printed."
  end
  printstyled("\nA Linear Mixed Model was build using model equations:\n\n",bold=true)
  for i in model.modelVec
    println(i)
  end
  println()
  printstyled("Model Information:\n\n",bold=true)
  @printf("%-15s %-12s %-10s %11s\n","Term","C/F","F/R","nLevels")

  random_effects=Array{AbstractString,1}()
  if model.pedTrmVec != 0
    for i in model.pedTrmVec
        push!(random_effects,split(i,':')[end])
    end
  end
  for i in model.rndTrmVec
      for j in i.term_array
          push!(random_effects,split(j,':')[end])
      end
  end

  terms=[]
  for i in model.modelTerms
    term    = split(i.trmStr,':')[end]
    if term in terms
        continue
    else
        push!(terms,term)
    end

    nLevels = i.nLevels
    fixed   = (term in random_effects) ? "random" : "fixed"
    factor  = (nLevels==1) ? "covariate" : "factor"

    if size(model.mmeRhs)==()&&data==false
      nLevels="NA"
    elseif size(model.mmeRhs)==()
      getMME(model,data)
    end

    if term =="intercept"
        factor="factor"
    elseif length(split(term,'*'))!=1
        factor="interaction"
    end

    @printf("%-15s %-12s %-10s %11s\n",term,factor,fixed,nLevels)
  end
  println()
  getMCMCinfo(model)
end

"""
    getMCMCinfo(model::MME)

* (internal function) Print out MCMC information.
"""
function getMCMCinfo(mme)
    if mme.MCMCinfo == 0
        printstyled("MCMC information is not available\n\n",bold=true)
        return
    end
    MCMCinfo = mme.MCMCinfo
    printstyled("MCMC Information:\n\n",bold=true)
    @printf("%-30s %20s\n","chain_length",MCMCinfo.chain_length)
    @printf("%-30s %20s\n","burnin",MCMCinfo.burnin)
    @printf("%-30s %20s\n","starting_value",mme.sol != false ? "true" : "false")
    @printf("%-30s %20d\n","printout_frequency",MCMCinfo.printout_frequency)
    @printf("%-30s %20d\n","output_samples_frequency",MCMCinfo.output_samples_frequency)
    @printf("%-30s %20s\n","constraint",MCMCinfo.constraint ? "true" : "false")
    @printf("%-30s %20s\n","missing_phenotypes",MCMCinfo.missing_phenotypes ? "true" : "false")
    @printf("%-30s %20d\n","update_priors_frequency",MCMCinfo.update_priors_frequency)
    @printf("%-30s %20s\n","seed",MCMCinfo.seed)

    printstyled("\nHyper-parameters Information:\n\n",bold=true)
    for i in mme.rndTrmVec
        thisterm= join(i.term_array, ",")
        if mme.nModels == 1
            @printf("%-30s %20s\n","random effect variances ("*thisterm*"):",Float64.(round.(inv(i.GiNew),digits=3)))
        else
            @printf("%-30s\n","random effect variances ("*thisterm*"):")
            Base.print_matrix(stdout,round.(inv(i.GiNew),digits=3))
            println()
        end
    end
    if mme.pedTrmVec!=0
        @printf("%-30s\n","genetic variances (polygenic):")
        Base.print_matrix(stdout,round.(inv(mme.Gi),digits=3))
        println()
    end
    if mme.nModels == 1
        @printf("%-30s %20.3f\n","residual variances:", (MCMCinfo.categorical_trait ? 1.0 : mme.R))
    else
        @printf("%-30s\n","residual variances:")
        Base.print_matrix(stdout,round.(mme.R,digits=3))
        println()
    end

    printstyled("\nGenomic Information:\n\n",bold=true)
    if mme.M != false
        print(MCMCinfo.single_step_analysis ? "incomplete genomic data" : "complete genomic data")
        println(MCMCinfo.single_step_analysis ? " (i.e., single-step analysis)" : " (i.e., non-single-step analysis)")
        for Mi in mme.M
            println()
            @printf("%-30s %20s\n","Genomic Category", Mi.name)
            @printf("%-30s %20s\n","Method",Mi.method)
            for Mi in mme.M
                if Mi.genetic_variance != false
                    if mme.nModels == 1
                        @printf("%-30s %20.3f\n","genetic variances (genomic):",Mi.genetic_variance)
                    else
                        @printf("%-30s\n","genetic variances (genomic):")
                        Base.print_matrix(stdout,round.(Mi.genetic_variance,digits=3))
                        println()
                    end
                end
                if !(Mi.method in ["GBLUP"])
                    if mme.nModels == 1
                        @printf("%-30s %20.3f\n","marker effect variances:",Mi.G)
                    else
                        @printf("%-30s\n","marker effect variances:")
                        Base.print_matrix(stdout,round.(Mi.G,digits=3))
                        println()
                    end
                end
                if !(Mi.method in ["RR-BLUP","BayesL","GBLUP"])
                    if mme.nModels == 1
                        @printf("%-30s %20s\n","π",Mi.π)
                    else
                        println("\nΠ: (Y(yes):included; N(no):excluded)\n")
                        print(string.(mme.lhsVec))
                        @printf("%20s\n","probability")
                        for (i,j) in Mi.π
                            i = replace(string.(i),"1.0"=>"Y","0.0"=>"N")
                            print(i)
                            @printf("%20s\n",j)
                        end
                        println()
                    end
                    @printf("%-30s %20s\n","estimatePi",Mi.estimatePi ? "true" : "false")
                end
                @printf("%-30s %20s\n","estimateScale",Mi.estimateScale ? "true" : "false")
            end
        end
    end
    printstyled("\nDegree of freedom for hyper-parameters:\n\n",bold=true)
    @printf("%-30s %20.3f\n","residual variances:",mme.df.residual)
    for randomeffect in mme.rndTrmVec
        if randomeffect.randomType != "A"
            @printf("%-30s %20.3f\n","random effect variances:",randomeffect.df)
        end
    end
    if mme.pedTrmVec!=0
        @printf("%-30s %20.3f\n","polygenic effect variances:",mme.df.polygenic)
    end
    if mme.M!=0
        @printf("%-30s %20.3f\n","marker effect variances:",mme.df.marker)
    end
    @printf("\n\n\n")
end
