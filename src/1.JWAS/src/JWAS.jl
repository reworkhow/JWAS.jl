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
include("buildMME/types.jl")
include("buildMME/build_MME.jl")
include("buildMME/random_effects.jl")
include("buildMME/residual.jl")
include("buildMME/sample_variances.jl")
include("buildMME/solver.jl")

#Markov chain Monte Carlo
include("MCMC/MCMC_BayesianAlphabet.jl")
include("MCMC/MT_MCMC_BayesianAlphabet.jl")

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
include("SSBR/SSBR.jl")

#Categorical trait
include("CategoricalTrait/categorical_trait.jl")

#Structure Equation Models
include("StructureEquationModel/SEM.jl")

#output
include("output.jl")

export build_model,set_covariate,set_random,add_genotypes,get_genotypes
export outputMCMCsamples,outputEBV,getEBV
export solve,runMCMC
export showMME,describe
#Pedmodule
export get_pedigree,getinfo
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
            output_samples_file      = "MCMC_samples",
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
    * Save MCMC samples every `output_samples_frequency` iterations, defaulting to `chain_length/1000`, to files `output_samples_file`,
      defaulting to `MCMC_samples.txt`. MCMC samples for hyperparametes (variance componets) and marker effects are saved by default.
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
        * If `constraint`=true, defaulting to `false`, constrain residual covariances between traits to be zeros.
        * If `causal_structure` is provided, e.g., causal_structure = [0.0,0.0,0.0;1.0,0.0,0.0;1.0,0.0,0.0] for
          trait 2 -> trait 1 and trait 3 -> trait 1, phenotypic causal networks will be incorporated using structure equation models.
* Genomic Prediction
    * Individual estimted breeding values (EBVs) and prediction error variances (PEVs) are returned if `outputEBV`=true, defaulting to `true`. Heritability and genetic
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
                categorical_trait               = false,
                missing_phenotypes              = true,
                constraint                      = false,
                causal_structure                = false,
                #Genomic Prediction
                outputEBV                       = true,
                output_heritability             = true,  #complete or incomplete genomic data
                #MISC
                seed                            = false,
                printout_model_info             = true,
                printout_frequency              = chain_length+1,
                big_memory                      = false,
                double_precision                = false,
                #MCMC samples (defaut to marker effects and hyperparametes (variance componets))
                output_samples_file             = "MCMC_samples",
                output_samples_for_all_parameters = false,
                #for deprecated JWAS
                methods                         = "conventional (no markers)",
                Pi                              = 0.0,
                estimatePi                      = false,
                estimateScale                   = false) #calculare lsmeans

    ############################################################################
    # Pre-Check
    ############################################################################
    if seed != false
        Random.seed!(seed)
    end
    if output_samples_for_all_parameters == true
        allparameters=[term[2] for term in split.(keys(mme.modelTermDict),":")]
        allparameters=unique(allparameters)
        for parameter in allparameters
            outputMCMCsamples(mme,parameter)
        end
    end
    if causal_structure != false
        #no missing phenotypes and residual covariance for identifiability
        missing_phenotypes, constraint = false, true
        if !istril(causal_structure)
            error("The causal structue needs to be a lower triangular matrix.")
        end
    end
    mme.MCMCinfo = MCMCinfo(chain_length,
                            starting_value,
                            burnin,
                            output_samples_frequency,
                            printout_model_info,
                            printout_frequency,
                            single_step_analysis, #pedigree,
                            missing_phenotypes,
                            constraint,
                            estimate_variance,
                            update_priors_frequency,
                            outputEBV,
                            output_heritability,
                            categorical_trait,
                            seed,
                            double_precision,
                            output_samples_file)
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
    # Check Arguments, Pedigree, Phenotypes, and output individual IDs (before align_genotypes)
    ############################################################################
    #check errors in function arguments
    errors_args(mme)
    #users need to provide high-quality pedigree file
    #println("calling check_pedigree(mme,df,pedigree)")
    check_pedigree(mme,df,pedigree)
    check_genotypes(mme)
    #user-defined IDs to return genetic values (EBVs), defaulting to all genotyped
    #individuals in genomic analysis
    #println("calling check_outputID(mme)")
    check_outputID(mme)
    #check phenotypes, only use phenotypes for individuals in pedigree
    #(incomplete genomic data,PBLUP) or with genotypes (complete genomic data)
    #check heterogeneous residuals
    df = check_phenotypes(mme,df,heterogeneous_residuals)
    #check priors
    #make default covariance matrices if not provided
    set_default_priors_for_variance_components(mme,df)
    #double precision
    if double_precision == true
        if mme.M != 0
            for Mi in mme.M
                Mi.genotypes = map(Float64,Mi.genotypes)
                Mi.G         = map(Float64,Mi.G)
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
    ############################################################################
    #align genotypes with 1) phenotypes IDs; 2) output IDs.
    ############################################################################
    if mme.M!=0
        align_genotypes(mme,output_heritability,single_step_analysis)
    end
    ############################################################################
    # Incomplete Genomic Data (Single-Step)
    ############################################################################
    #1)reorder in A (pedigree) as ids for genotyped then non-genotyped inds
    #2)impute genotypes for non-genotyped individuals
    #3)add ϵ (imputation errors) and J as variables in data for non-genotyped inds
    if single_step_analysis == true
        SSBRrun(mme,df,big_memory)
    end
    ############################################################################
    # Initiate Mixed Model Equations for Non-marker Parts (run after SSBRrun for ϵ & J)
    ############################################################################
    # initiate Mixed Model Equations and check starting values
    mme.MCMCinfo.starting_value,df = init_mixed_model_equations(mme,df,starting_value)

    if mme.M!=0
        set_marker_hyperparameters_variances_and_pi(mme)
    end

    if mme.output_ID!=0
        get_outputX_others(mme,single_step_analysis)
    end

    #printout basic MCMC information
    if printout_model_info == true
      describe(mme)
    end
    ############################################################################
    # Double Precision or not
    ############################################################################
    if mme.nModels ==1 #single-trait analysis
        mme.output=MCMC_BayesianAlphabet(chain_length,mme,df,
                        Rinv                     = mme.invweights,
                        burnin                   = burnin,
                        estimate_variance        = estimate_variance,
                        starting_value           = mme.MCMCinfo.starting_value,
                        outFreq                  = printout_frequency,
                        output_samples_frequency = output_samples_frequency,
                        output_file              = output_samples_file,
                        update_priors_frequency  = update_priors_frequency,
                        categorical_trait        = categorical_trait)
    else #multi-trait analysis
        mme.output=MT_MCMC_BayesianAlphabet(chain_length,mme,df,
                      Rinv   = mme.invweights,
                      sol    = mme.MCMCinfo.starting_value,
                      outFreq= printout_frequency,
                      missing_phenotypes=missing_phenotypes,
                      constraint = constraint,
                      estimate_variance = estimate_variance,
                      output_samples_frequency=output_samples_frequency,
                      output_file=output_samples_file,
                      update_priors_frequency=update_priors_frequency,
                      causal_structure = causal_structure)
    end
    for (key,value) in mme.output
      CSV.write(replace(key," "=>"_")*".txt",value)
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
#
# Pre-check before running MCMC for
# 1) arguments
# 2) phenotypes
# 3) genotypes
# 4) pedigree
# 5) starting values
# 6) output IDs
#
################################################################################
function errors_args(mme)
    if mme.mmePos != 1
      error("Please build your model again using the function build_model().")
    end

    if mme.MCMCinfo.output_samples_frequency <= 0
        error("output_samples_frequency should be an integer > 0.")
    end

    if mme.M != 0
        for Mi in mme.M
            if !(Mi.method in ["BayesL","BayesC","BayesB","BayesA","RR-BLUP","GBLUP"])
                error(Mi.method," is not available in JWAS. Please read the documentation.")
            end

            if Mi.method in ["RR-BLUP","BayesL","GBLUP","BayesA"]
                if Mi.π != false
                    error(Mi.method," runs with π = false.")
                elseif Mi.estimatePi == true
                    error(Mi.method," runs with estimatePi = false.")
                end
            end
            if Mi.method == "BayesA"
                Mi.method = "BayesB"
                println("BayesA is equivalent to BayesB with known π=0. BayesB with known π=0 runs.")
            end
            if Mi.method == "GBLUP"
                if Mi.genetic_variance == false
                    error("Please provide values for the genetic variance for GBLUP analysis")
                end
                if mme.MCMCinfo.single_step_analysis == true
                    error("SSGBLUP is not available")
                end
            end
            if mme.nModels > 1 && Mi.π != 0.0
                if round(sum(values(Mi.π)),digits=2)!=1.0
                  error("Summation of probabilities of Pi is not equal to one.")
                end
                if typeof(Mi.π) <: Number
                    error("Pi cannot be a number in multi-trait analysis.")
                end
            end
        end
    end
    if mme.MCMCinfo.single_step_analysis == true && mme.M == 0
        error("Genomic information is required for single-step analysis.")
    end
end

function check_pedigree(mme,df,pedigree)
    if mme.ped == 0 && pedigree != false #check whether mme.ped exists
        mme.ped = deepcopy(pedigree)
    end
    if mme.ped != 0
        pedID=map(string,collect(keys(mme.ped.idMap)))
        if mme.M!=0
            for Mi in mme.M
                if !issubset(Mi.obsID,pedID)
                    error("Not all genotyped individuals are found in pedigree!")
                end
            end
        end
    end
end

function check_genotypes(mme)
    if mme.M != 0
        obsID = mme.M[1].obsID
        for Mi in mme.M
            if obsID != Mi.obsID
                error("genotypic information is not provided for same individuals")
            end
        end
    end
    if mme.M != 0 && mme.MCMCinfo.single_step_analysis == true
        if length(mme.M) != 1
            error("Now only one genomic category is allowed in single-step analysis.")
        end
    end
end

function check_outputID(mme)
    #Genotyped individuals are usaully not many, and are used in GWAS (complete
    #and incomplete), thus are used as default output_ID if not provided
    if mme.M == 0 && mme.pedTrmVec == 0
        mme.MCMCinfo.outputEBV = false #no EBV in non-genetic analysis
    end

    if mme.MCMCinfo.outputEBV == false
        mme.output_ID = 0
    elseif mme.output_ID == 0 && mme.M != 0 #all genotyped inds if no output ID
        mme.output_ID = copy(mme.M[1].obsID)
    elseif mme.output_ID == 0 && mme.M == 0 && mme.pedTrmVec != 0 #all inds in PBLUP
        mme.output_ID = copy(mme.ped.IDs)
    end

    single_step_analysis = mme.MCMCinfo.single_step_analysis
    if single_step_analysis == false && mme.M != 0 #complete genomic data
        if mme.output_ID!=0 && !issubset(mme.output_ID,mme.M[1].obsID)
            printstyled("Testing individuals are not a subset of ",
            "genotyped individuals (complete genomic data,non-single-step). ",
            "Only output EBV for tesing individuals with genotypes.\n",bold=false,color=:red)
            mme.output_ID = intersect(mme.output_ID,mme.M[1].obsID)
        end
    elseif mme.ped != false #1)incomplete genomic data 2)PBLUP
        pedID = map(string,collect(keys(mme.ped.idMap)))
        if mme.output_ID!=0 && !issubset(mme.output_ID,pedID)
            printstyled("Testing individuals are not a subset of ",
            "individuals in pedigree (incomplete genomic data (single-step) or PBLUP). ",
            "Only output EBV for tesing individuals in the pedigree.",bold=false,color=:red)
            mme.output_ID = intersect(mme.output_ID,pedID)
        end
    end
    #Set ouput IDs to all genotyped inds in complete genomic data analysis for h² estimation
    if mme.MCMCinfo.output_heritability == true && mme.MCMCinfo.single_step_analysis == false && mme.M != 0
        mme.output_ID = copy(mme.M[1].obsID)
    end

end

function check_phenotypes(mme,df,heterogeneous_residuals)
    printstyled("Checking phenotypes...\n" ,bold=false,color=:green)
    df[!,1] = strip.(map(string,df[!,1])) #make IDs stripped string
    printstyled("Individual IDs (strings) are provided in the first column of the phenotypic data.\n" ,bold=false,color=:green)
    writedlm("IDs_for_individuals_with_phenotypes.txt",df[!,1])

    #remove records whose phenotypes are missing for all traits fitted in the model
    missingdf  = ismissing.(convert(Matrix,df[!,mme.lhsVec]))
    allmissing = fill(true,mme.nModels)
    nonmissingindex = Array{Int64,1}()
    for i in 1:size(missingdf,1)
        if missingdf[i,:] != allmissing
            push!(nonmissingindex,i)
        else
            printstyled("Phenotypes for all traits included in the model for individual ",df[!,1][i], " in the row ",i," are missing. This record is deleted.\n" ,bold=false,color=:red)
        end
    end
    if length(nonmissingindex) != 0
        df = df[nonmissingindex,:]
    end

    phenoID = df[!,1]
    single_step_analysis = (mme.MCMCinfo != false ? mme.MCMCinfo.single_step_analysis : false)
    if single_step_analysis == false && mme.M != 0 #complete genomic data
        if !issubset(phenoID,mme.M[1].obsID)
            printstyled("Phenotyped individuals are not a subset of ",
            "genotyped individuals (complete genomic data,non-single-step). ",
            "Only use phenotype information for genotyped individuals.",bold=false,color=:red)
            index = [phenoID[i] in mme.M[1].obsID for i=1:length(phenoID)]
            df    = df[index,:]
        end
        printstyled("The number of observations with both genotypes and phenotypes used ",
        "in the analysis is ",size(df,1),".\n",bold=false,color=:green)
    end
    if mme.ped != false #1)incomplete genomic data 2)PBLUP or 3)complete genomic data with polygenic effect
        pedID = map(string,collect(keys(mme.ped.idMap)))
        if !issubset(phenoID,pedID)
            printstyled("Phenotyped individuals are not a subset of ",
            "individuals in pedigree (incomplete genomic data (single-step) or PBLUP). ",
            "Only use phenotype information for individuals in the pedigree.",bold=false,color=:red)
            index = [phenoID[i] in pedID for i=1:length(phenoID)]
            df    = df[index,:]
        end
        printstyled("The number of observations with both phenotype and pedigree information ",
        "used in the analysis is ",size(df,1),".\n",bold=false,color=:green)
    end
    #***************************************************************************
    # set IDs for phenotypes
    #***************************************************************************
    mme.obsID  = map(string,df[!,1])
    #***************************************************************************
    # set Riv for heterogeneous residuals
    #***************************************************************************
    if heterogeneous_residuals == true
        invweights = 1 ./ convert(Array,df[!,Symbol("weights")])
    else
        invweights = ones(length(mme.obsID))
    end
    mme.invweights = (mme.MCMCinfo == false || mme.MCMCinfo.double_precision ? Float64.(invweights) : Float32.(invweights))
    return df
end

#set default covariance matrices for variance components if not provided
function set_default_priors_for_variance_components(mme,df)
  myvar     = [var(skipmissing((df[!,mme.lhsVec[i]]))) for i=1:size(mme.lhsVec,1)]
  phenovar  = diagm(0=>myvar)
  h2        = 0.5

  genetic_random_count    = (mme.M!=0 ? length(mme.M) : 0)
  nongenetic_random_count = 1
  for random_term in mme.rndTrmVec
    if random_term.randomType == "A"
        genetic_random_count += 1
    elseif split(random_term.term_array[1],":")[2] != "ϵ"
        nongenetic_random_count +=1
    end
  end
  if genetic_random_count != 0
      varg = phenovar*h2/genetic_random_count
  end
  vare = phenovar*h2/nongenetic_random_count

  #genetic variance or marker effect variance
  if mme.M!=0
    for Mi in mme.M
      if Mi.G == false && Mi.genetic_variance == false
          printstyled("Prior information for genomic variance is not provided and is generated from the data.\n",bold=false,color=:green)
          if mme.nModels==1
              Mi.genetic_variance = varg[1,1]
          elseif mme.nModels>1
              Mi.genetic_variance = varg
          end #mme.M.G and its scale parameter will be reset in function set_marker_hyperparameters_variances_and_pi
      end
    end
  end
  #residual effects
  if mme.nModels==1 && isposdef(mme.RNew) == false #single-trait
    printstyled("Prior information for residual variance is not provided and is generated from the data.\n",bold=false,color=:green)
    mme.RNew = mme.ROld = vare[1,1]
    mme.scaleRes = mme.RNew*(mme.df.residual-2)/mme.df.residual
  elseif mme.nModels>1 && isposdef(mme.R) == false #multi-trait
    printstyled("Prior information for residual variance is not provided and is generated from the data.\n",bold=false,color=:green)
    mme.R = vare
    mme.scaleRes = mme.R*(mme.df.residual - mme.nModels - 1)
  end
  #random effects
  if length(mme.rndTrmVec) != 0
    for randomEffect in mme.rndTrmVec
      if isposdef(randomEffect.Gi) == false
        printstyled("Prior information for random effect variance is not provided and is generated from the data.\n",bold=false,color=:green)
        myvarout  = [split(i,":")[1] for i in randomEffect.term_array]
        myvarin   = string.(mme.lhsVec)
        Zdesign   = mkmat_incidence_factor(myvarout,myvarin)
        if randomEffect.randomType == "A" || split(randomEffect.term_array[1],":")[2] == "ϵ"
            G = diagm(Zdesign*diag(varg))
        else
            G = diagm(Zdesign*diag(vare))
        end
        randomEffect.Gi = randomEffect.GiOld = randomEffect.GiNew = Symmetric(inv(G))
        randomEffect.scale = G*(randomEffect.df-length(randomEffect.term_array)-1)
        if randomEffect.randomType == "A"
          mme.Gi = randomEffect.Gi
          mme.scalePed = randomEffect.scale
        end
      end
    end
  end
end


function init_mixed_model_equations(mme,df,sol)
    getMME(mme,df)
    #starting value for sol can be provided
    if mme.M == 0                          #PBLUP
        nsol = size(mme.mmeLhs,1)
    else
        nsol = size(mme.mmeLhs,1)+sum((Mi.method != "GBLUP" ? Mi.nMarkers : Mi.nObs) for Mi in mme.M)*mme.nModels
    end
    if sol == false #no starting values
        sol = zeros((mme.MCMCinfo.double_precision ? Float64 : Float32),nsol)
    else            #besure type is Float64
        printstyled("Starting values are provided. The order of starting values for location parameters and\n",
        "marker effects should be the order of location parameters in the Mixed Model Equation for all traits (This can be\n",
        "obtained by getNames(model)) and then markers for all traits (all markers for trait 1 then all markers for trait 2...)\n",bold=false,color=:green)
        if length(sol) != nsol || typeof(sol) <: AbstractString
            error("length or type of starting values is wrong.")
        end
        sol = map((mme.MCMCinfo.double_precision ? Float64 : Float32),sol)
    end
    return sol,df
end

################################################################################
#
#
# Print out Model or MCMC information
#
#
################################################################################
"""
    describe(model::MME)

* Print out model information.
"""
function describe(model::MME;data=false)
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
    @printf("%-30s %20s\n","starting_value",MCMCinfo.starting_value != zero(MCMCinfo.starting_value) ? "true" : "false")
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
        @printf("%-30s %20.3f\n","residual variances:", (MCMCinfo.categorical_trait ? 1.0 : mme.RNew))
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
