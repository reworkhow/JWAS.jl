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

#Categorical and Censored Traits
include("categorical_and_censored_trait/categorical_and_censored_trait.jl")

#Structure Equation Models
include("structure_equation_model/SEM.jl")

#Random Regression Model
include("RRM/MCMC_BayesianAlphabet_RRM.jl")
include("RRM/RRM.jl")

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
export GWAS
#dataset
export dataset
#export adjust_phenotypes,LOOCV

"""
    runMCMC(model::MME,df::DataFrame;
            #Data
            heterogeneous_residuals           = false,
            #MCMC
            chain_length::Integer             = 100,
            starting_value                    = false,
            burnin::Integer                   = 0,
            output_samples_frequency::Integer = chain_length>1000 ? div(chain_length,1000) : 1,
            update_priors_frequency::Integer  = 0,
            #Methods
            single_step_analysis            = false, #parameters for single-step analysis
            pedigree                        = false, #parameters for single-step analysis
            fitting_J_vector                = true,  #parameters for single-step analysis
            categorical_trait               = false,
            censored_trait                  = false,
            causal_structure                = false,
            mega_trait                      = false,
            missing_phenotypes              = true,
            constraint                      = false,
            #Genomic Prediction
            outputEBV                       = true,
            output_heritability             = true,
            prediction_equation             = false,
            #MISC
            seed                            = false,
            printout_model_info             = true,
            printout_frequency              = chain_length+1,
            big_memory                      = false,
            double_precision                = false,
            ##MCMC samples (defaut to marker effects and hyperparametes (variance components))
            output_folder                     = "results",
            output_samples_for_all_parameters = false,
            ##for deprecated JWAS
            methods                         = "conventional (no markers)",
            Pi                              = 0.0,
            estimatePi                      = false)

**Run MCMC for Bayesian Linear Mixed Models with or without estimation of variance components.**

* Markov chain Monte Carlo
    * The first `burnin` iterations are discarded at the beginning of a MCMC chain of length `chain_length`.
    * Save MCMC samples every `output_samples_frequency` iterations, defaulting to `chain_length/1000`, to a folder `output_folder`,
      defaulting to `results`. MCMC samples for hyperparametes (variance componets) and marker effects are saved by default.
      MCMC samples for location parametes can be saved using function `output_MCMC_samples()`. Note that saving MCMC samples too
      frequently slows down the computation.
    * The `starting_value` can be provided as a vector for all location parameteres and marker effects, defaulting to `0.0`s.
      The order of starting values for location parameters and marker effects should be the order of location parameters in
      the Mixed Model Equation for all traits (This can be obtained by getNames(model)) and then markers for all traits (all
      markers for trait 1 then all markers for trait 2...).
    * Miscellaneous Options
        * Priors are updated every `update_priors_frequency` iterations, defaulting to `0`.
* Methods
    * Single step analysis is allowed if `single_step_analysis` = `true` and `pedigree` is provided.
    * Miscellaneous Options
        * Missing phenotypes are allowed in multi-trait analysis with `missing_phenotypes`=true, defaulting to `true`.
        * Catogorical Traits are allowed if `categorical_trait`=true, defaulting to `false`. Phenotypes should be coded as 1,2,3...
        * Censored traits are allowed if the upper bounds are provided in `censored_trait` as an array, and lower bounds are provided as phenotypes.
        * If `constraint`=true, defaulting to `false`, constrain residual covariances between traits to be zeros.
        * If `causal_structure` is provided, e.g., causal_structure = [0.0 0.0 0.0;1.0 0.0 0.0;1.0 0.0 0.0] for
          trait 1 -> trait 2 and trait 1 -> trait 3 (column index affacts row index, and a lower triangular matrix is required), phenotypic causal networks will be incorporated using structure equation models.
* Genomic Prediction
    * Predicted values for individuals of interest can be obtained based on a user-defined prediction equation `prediction_equation`, e.g., "y1:animal + y1:age".
    For now, genomic data is always included. Genetic values including effects defined with genotype and pedigree information are returned if `prediction_equation`= false, defaulting to `false`.
    * Individual estimted genetic values and prediction error variances (PEVs) are returned if `outputEBV`=true, defaulting to `true`. Heritability and genetic
    variances are returned if `output_heritability`=`true`, defaulting to `true`. Note that estimation of heritability is computaionally intensive.
* Miscellaneous Options
  * Print out the model information in REPL if `printout_model_info`=true; print out the monte carlo mean in REPL with `printout_frequency`,
    defaulting to `false`.
  * If `seed`, defaulting to `false`, is provided, a reproducible sequence of numbers will be generated for random number generation.
  * If `big_memory`=true, defaulting to `false`, a machine with  lots of memory is assumed which may speed up the analysis.
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
                single_step_analysis            = false, #parameters for single-step analysis
                pedigree                        = false, #parameters for single-step analysis
                fitting_J_vector                = true,  #parameters for single-step analysis
                causal_structure                = false,
                missing_phenotypes              = true,
                # constraint                      = false,
                RRM                             = false, #  RRM -> false or a matrix (Phi: orthogonalized time covariates)
                #Genomic Prediction
                outputEBV                       = true,
                output_heritability             = true,
                prediction_equation             = false,
                #MISC
                seed                            = false,
                printout_model_info             = true,
                printout_frequency              = chain_length+1,
                big_memory                      = false,
                double_precision                = false,
                fast_blocks                     = false,
                #MCMC samples (defaut to marker effects and hyperparametes (variance componets))
                output_folder                     = "results",
                output_samples_for_all_parameters = false,
                #for deprecated JWAS
                methods                         = "conventional (no markers)",
                Pi                              = 0.0,
                estimatePi                      = false,
                estimate_scale                   = false, #this has been moved to get_genotypes()
                estimate_variance               = true,   #this has been moved to build_MME() and set_random()
                categorical_trait               = false,  #this has been moved to build_model()
                censored_trait                  = false)  #this has been moved to build_model()

    #for deprecated JWAS fucntions
    ############################################################################
    if mme.M != 0
        for Mi in mme.M
            if Mi.name == false
                Mi.name              = "geno"
                Mi.π                 = Pi
                Mi.estimatePi        = estimatePi
                Mi.G.estimate_scale  = estimate_scale
                Mi.method            = methods
            end
        end
    end
    if categorical_trait != false || censored_trait != false
        print_single_categorical_censored_trait_example()
        error("The arguments 'categorical_trait' and  'censored_trait' has been moved to build_model(). Please check our latest example.")
    end
    if estimate_scale != false #user set estimate_variance=true in runMCMC()
        error("The argument 'estimate_scale' for marker effect variance has been moved to get_genotypes().")
    end
    if estimate_variance != true #user set estimate_variance=false in runMCMC()
        error("The argument 'estimate_variance' for non-marker variance components has been moved to build_MME() for residual variance,
               and set_random() for random terms.")
    end
    ############################################################################
    # censored traits
    ############################################################################
    # Goal: add the column named "traitname" using trait's lower bound and upper
    #       bound, to avoid the error that mme.lhsVec is not in df
    if "censored" ∈ mme.traits_type
        add_censored_trait_column!(mme,df) #changed: df
    end
    ############################################################################
    # Set a seed in the random number generator
    ############################################################################
    if seed != false
        Random.seed!(seed)
    end
    #when using multi-thread, make sure the results are reproducible for users
    nThread = Threads.nthreads()
    if nThread>1
        Threads.@threads for i = 1:nThread
              Random.seed!(seed+i)
        end
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
    mme.MCMCinfo = MCMCinfo(heterogeneous_residuals,
                   chain_length,burnin,output_samples_frequency,
                   printout_model_info,printout_frequency, single_step_analysis,
                   fitting_J_vector,missing_phenotypes,
                   update_priors_frequency,outputEBV,output_heritability,prediction_equation,
                   seed,double_precision,output_folder,RRM,fast_blocks)
    ############################################################################
    # Check 1)Arguments; 2)Input Pedigree,Genotype,Phenotypes,
    #       3)output individual IDs; 4)Priors  5)Prediction equation
    #       in the Whole Dataset.
    ############################################################################
    errors_args(mme)       #check errors in function arguments
    df=check_pedigree_genotypes_phenotypes(mme,df,pedigree)
    prediction_setup(mme)  #set prediction equation, defaulting to genetic values
    check_outputID(mme)    #check individual of interest for prediction
    df_whole,train_index = make_dataframes(df,mme)
    set_default_priors_for_variance_components(mme,df_whole)  #check priors (set default priors)

    if mme.M!=0
        if single_step_analysis == true
            SSBRrun(mme,df_whole,train_index,big_memory)
        end
        set_marker_hyperparameters_variances_and_pi(mme)
    end
    ############################################################################
    #fast blocks #now only work for one geno
    ############################################################################
    if fast_blocks != false
        if fast_blocks == true
            block_size = Int(floor(sqrt(mme.M[1].nObs)))
        elseif typeof(fast_blocks) <: Number
            block_size = Int(floor(fast_blocks))
        end
        mme.MCMCinfo.fast_blocks  = collect(range(1, step=block_size, stop=mme.M[1].nMarkers))
        mme.MCMCinfo.chain_length = Int(floor(chain_length/(mme.MCMCinfo.fast_blocks[2]-mme.MCMCinfo.fast_blocks[1]))) #number of outer loop
        println("BLOCK SIZE: $block_size")
        flush(stdout)
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
        mme.MCMCinfo.missing_phenotypes, mme.MCMCinfo.constraint = false, true
        if !istril(causal_structure)
            error("The causal structue needs to be a lower triangular matrix.")
        end
    end
    # Double Precision
    # Only convert when needed to avoid unnecessary full-size allocations.
    tofloat64_if_needed = x -> begin
        if x isa AbstractArray
            return (eltype(x) == Float64) ? x : map(Float64, x)
        elseif x isa Float64
            return x
        else
            return map(Float64, x)
        end
    end
    if double_precision == true
        if mme.M != 0
            for Mi in mme.M
                Mi.genotypes = tofloat64_if_needed(Mi.genotypes)
                Mi.G.val     = tofloat64_if_needed(Mi.G.val)
                Mi.α         = tofloat64_if_needed(Mi.α)
            end
        end
        for random_term in mme.rndTrmVec
            random_term.Vinv      = tofloat64_if_needed(random_term.Vinv)
            random_term.GiOld.val = tofloat64_if_needed(random_term.GiOld.val)
            random_term.GiNew.val = tofloat64_if_needed(random_term.GiNew.val)
            random_term.Gi.val    = tofloat64_if_needed(random_term.Gi.val)
        end
        mme.Gi = tofloat64_if_needed(mme.Gi)
    end

    #constraint on covariance matrix
    if mme.R.constraint==true
        R_constraint!(mme) #modify mme.R.df and mme.R.scale; scale is a diagonal matrix
    end
    if mme.M != 0 && mme.M[1].G.constraint==true
        G_constraint!(mme) #modify Mi.G.df and Mi.G.scale; scale is a diagonal matrix
    end

    ############################################################################
    #Make incidence matrices and genotype covariates for training observations
    #and individuals of interest
    ############################################################################
    #make incidence matrices (non-genomic effects) (after SSBRrun for ϵ & J)
    df=make_incidence_matrices(mme,df_whole,train_index)
    #align genotypes with 1) phenotypes IDs; 2) output IDs.
    if mme.M != false
        align_genotypes(mme,output_heritability,single_step_analysis)
    end
    # initiate Mixed Model Equations and check starting values
    init_mixed_model_equations(mme,df,starting_value)
    ############################################################################
    # MCMC
    ############################################################################
    describe(mme)
    if RRM == false
        mme.output=MCMC_BayesianAlphabet(mme,df)
    else
        mme.output=MCMC_BayesianAlphabet_RRM(mme,df,
                                  Φ                        = RRM,
                                  burnin                   = burnin,
                                  outFreq                  = printout_frequency,
                                  output_samples_frequency = output_samples_frequency,
                                  update_priors_frequency  = update_priors_frequency,
                                  output_folder            = mme.MCMCinfo.output_folder,
                                  nIter                    = chain_length)
    end

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

    # make MCMC samples for indirect marker effect
    if causal_structure != false

        #generate the MCMC sample file for indirect and direct effect file.
        generate_indirect_marker_effect_sample(mme.lhsVec,output_folder,causal_structure,"structure_coefficient_MCMC_samples.txt")
        generate_overall_marker_effect_sample(mme.lhsVec,output_folder,causal_structure)

        # generate marker effet file for direct, indirect, and overall effect
        generate_marker_effect(mme.lhsVec, output_folder,causal_structure,"direct")
        generate_marker_effect(mme.lhsVec, output_folder,causal_structure,"indirect")
        generate_marker_effect(mme.lhsVec, output_folder,causal_structure,"overall")
    end

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
    if size(model.mmeRhs)==() && data != false
      solve(model,data)
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

    if term =="intercept"
        factor="factor"
    elseif length(split(term,'*'))!=1
        factor="interaction"
    end

    @printf("%-15s %-12s %-10s %11s\n",term,factor,fixed,nLevels)
    end
    println()
    if model.MCMCinfo != false && model.MCMCinfo.printout_model_info == true
        getMCMCinfo(model)
    end
end

"""
    getMCMCinfo(model::MME)

* (internal function) Print out MCMC information.
"""
function getMCMCinfo(mme)
    if mme.MCMCinfo == false
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
    # @printf("%-30s %20s\n","constraint",MCMCinfo.constraint ? "true" : "false")
    @printf("%-30s %19s\n","constraint on residual variance",mme.R.constraint ? "true" : "false")
    if mme.M != 0
        for Mi in mme.M
            geno_name = Mi.name
            @printf("%-30s %5s\n","constraint on marker effect variance for $geno_name",Mi.G.constraint ? "true" : "false")
        end
    end
    @printf("%-30s %20s\n","missing_phenotypes",MCMCinfo.missing_phenotypes ? "true" : "false")
    @printf("%-30s %20d\n","update_priors_frequency",MCMCinfo.update_priors_frequency)
    @printf("%-30s %20s\n","seed",MCMCinfo.seed)

    printstyled("\nHyper-parameters Information:\n\n",bold=true)
    for i in mme.rndTrmVec
        thisterm= join(i.term_array, ",")
        if mme.nModels == 1
            @printf("%-30s %20s\n","random effect variances ("*thisterm*"):",Float64.(round.(inv(i.GiNew.val),digits=3)))
        else
            @printf("%-30s\n","random effect variances ("*thisterm*"):")
            Base.print_matrix(stdout,round.(inv(i.GiNew.val),digits=3))
            println()
        end
    end
    if mme.pedTrmVec!=0
        polygenic_pos = findfirst(i -> i.randomType=="A", mme.rndTrmVec)
    end
    if mme.pedTrmVec!=0
        @printf("%-30s\n","genetic variances (polygenic):")
        Base.print_matrix(stdout,round.(inv(mme.rndTrmVec[polygenic_pos].Gi.val),digits=3))
        println()
    end
    if mme.nModels == 1
        @printf("%-30s %20.3f\n","residual variances:", mme.R.val)
    else
        @printf("%-30s\n","residual variances:")
        Base.print_matrix(stdout,round.(mme.R.val,digits=3))
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
                if Mi.genetic_variance.val != false
                    if mme.nModels == 1 && mme.MCMCinfo.RRM == false
                        @printf("%-30s %20.3f\n","genetic variances (genomic):",Mi.genetic_variance.val)
                    else
                        @printf("%-30s\n","genetic variances (genomic):")
                        Base.print_matrix(stdout,round.(Mi.genetic_variance.val,digits=3))
                        println()
                    end
                end
                if !(Mi.method in ["GBLUP"])
                    if mme.nModels == 1 && mme.MCMCinfo.RRM == false
                        @printf("%-30s %20.3f\n","marker effect variances:",Mi.G.val)
                    else
                        @printf("%-30s\n","marker effect variances:")
                        Base.print_matrix(stdout,round.(Mi.G.val,digits=3))
                        println()
                    end
                end
                if !(Mi.method in ["RR-BLUP","BayesL","GBLUP"])
                    if mme.nModels == 1 && mme.MCMCinfo.RRM == false
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
                @printf("%-30s %20s\n","estimate_scale",Mi.G.estimate_scale ? "true" : "false")
            end
        end
    end
    printstyled("\nDegree of freedom for hyper-parameters:\n\n",bold=true)
    @printf("%-30s %20.3f\n","residual variances:",mme.R.df)
    for randomeffect in mme.rndTrmVec
        if randomeffect.randomType != "A"
            @printf("%-30s %20.3f\n","random effect variances:",randomeffect.Gi.df)
        end
    end
    if mme.pedTrmVec!=0
        # polygenic_pos = findfirst(i -> i.randomType=="A", mme.rndTrmVec)
        @printf("%-30s %20.3f\n","polygenic effect variances:",mme.rndTrmVec[polygenic_pos].Gi.df)
    end
    if mme.M!=0
        for Mi in mme.M
            @printf("%-30s %20.3f\n","marker effect variances:",Mi.G.df)
        end
    end
    @printf("\n\n\n")
end
