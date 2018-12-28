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
    # Single-Step
    ############################################################################
    #impute genotypes for non-genotyped individuals
    #add imputation errors and J as variables in data for non-genotyped individuals
    if single_step_analysis == true
        SSBRrun(mme,pedigree,df)
    end
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
    #make mixed model equations for non-marker parts
    starting_value,df = pre_check(mme,df,starting_value)

    if mme.M!=0
        #align genotypes with 1) phenotypes IDs; 2) output IDs.
        align_genotypes(mme)
    end
    if mme.M!=0
        #set up Pi
        if mme.nModels !=1 && Pi==0.0
            printstyled("Pi (Π) is not provided.\n",bold=false,color=:green)
            printstyled("Pi (Π) is generated assuming all markers have effects on all traits.\n",bold=false,color=:green)
            mykey=Array{Float64}(undef,0)
            ntraits=mme.nModels
            Pi=Dict{Array{Float64,1},Float64}()
            #for i in [ bin(n,ntraits) for n in 0:2^ntraits-1 ] `bin(n, pad)` is deprecated, use `string(n, base=2, pad=pad)
            for i in [ string(n,base=2,pad=ntraits) for n in 0:2^ntraits-1 ]
              Pi[parse.(Float64, split(i,""))]=0.0
            end
            Pi[ones(ntraits)]=1.0
        end
        #set up marker effect variances
        if mme.M.G_is_marker_variance==false && methods!="GBLUP"
            genetic2marker(mme.M,Pi)
            println()
            if mme.nModels != 1
              if !isposdef(mme.M.G) #also work for scalar
                error("Marker effects covariance matrix is not postive definite! Please modify the argument: Pi.")
              end
              print("The prior for marker effects covariance matrix is calculated from ")
              print("genetic covariance matrix and Π. The mean of the prior for the marker effects ")
              print("covariance matrix is: \n")
              Base.print_matrix(stdout,round.(mme.M.G,digits=6))
            else
              if !isposdef(mme.M.G) #positive scalar (>0)
                error("Marker effects variance is negative!")
              end
              print("The prior for marker effects variance is calculated from ")
              print("the genetic variance and π. The mean of the prior for the marker effects variance ")
              print("is: ",round.(mme.M.G,digits=6))
            end
        end
        println("\n\n")
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

################################################################################
#Print out MCMC information #Make all info a member of MME?
################################################################################
function MCMCinfo(methods,Pi,chain_length,burnin,starting_value,printout_frequency,
                  output_samples_frequency,missing_phenotypes,constraint,
                  estimatePi,update_priors_frequency,mme)

    println("MCMC Information:\n")
    @printf("%-30s %20s\n","methods",methods)
#    @printf("%-20s %10s\n","seed",seed)
    @printf("%-30s %20s\n","chain_length",chain_length)
    @printf("%-30s %20s\n","burnin",burnin)
    if !(methods in ["conventional (no markers)", "GBLUP"])
      @printf("%-30s %20s\n","estimatePi",estimatePi ? "true" : "false")
    end
    @printf("%-30s %20s\n","starting_value",starting_value ? "true" : "false")
    @printf("%-30s %20d\n","printout_frequency",printout_frequency)
    @printf("%-30s %20d\n","output_samples_frequency",output_samples_frequency)

    @printf("%-30s %20s\n","constraint",constraint ? "true" : "false")
    @printf("%-30s %20s\n","missing_phenotypes",missing_phenotypes ? "true" : "false")
    @printf("%-30s %20d\n","update_priors_frequency",update_priors_frequency)


    @printf("\n%-30s\n\n","Hyper-parameters Information:")

    if mme.nModels==1
        for i in mme.rndTrmVec
            thisterm=split(i.term_array[1],':')[end]
            @printf("%-30s %20s\n","random effect variances ("*thisterm*"):",round.(inv(i.GiNew),digits=3))
        end
        @printf("%-30s %20.3f\n","residual variances:",mme.RNew)
        if mme.pedTrmVec!=0
            @printf("%-30s\n %50s\n","genetic variances (polygenic):",round.(inv(mme.GiNew),digits=3))
        end
        if !(methods in ["conventional (no markers)", "GBLUP"])
            if mme.M == 0
                error("Please add genotypes using add_genotypes().")
            end
            #@printf("%-30s %20.3f\n","genetic variances (genomic):",mme.M.G)
            @printf("%-30s %20.3f\n","marker effect variances:",mme.M.G)
            @printf("%-30s %20s\n","π",Pi)
        end
    else
        for i in mme.rndTrmVec
            thisterm=split(i.term_array[1],':')[end]
            @printf("%-30s\n","random effect variances ("*thisterm*"):")
            Base.print_matrix(stdout,round.(inv(i.GiNew),digits=3))
            println()
        end
        @printf("%-30s\n","residual variances:")
        Base.print_matrix(stdout,round.(mme.R,digits=3))
        println()
        if mme.pedTrmVec!=0
            @printf("%-30s\n","genetic variances (polygenic):")
            Base.print_matrix(stdout,round.(inv(mme.Gi),digits=3))
            println()
        end
        if !(methods in ["conventional (no markers)", "GBLUP"])
            @printf("%-30s\n","genetic variances (genomic):")
            Base.print_matrix(stdout,round.(mme.M.G,digits=3))
            println()
            @printf("%-30s\n","marker effect variances:")
            Base.print_matrix(stdout,round.(mme.M.G,digits=3))
            println()
            println("\nΠ: (Y(yes):included; N(no):excluded)")
            @printf("%-20s %12s\n",string.(mme.lhsVec),"probability")
            for (i,j) in Pi
                i = replace(string.(i),"1.0"=>"Y","0.0"=>"N")
                @printf("%-20s %12s\n",i,j)
            end
        end
    end

    @printf("\n%-30s\n\n","Degree of freedom for hyper-parameters:")
    @printf("%-30s %20.3f\n","residual variances:",mme.df.residual)
    @printf("%-30s %20.3f\n","iid random effect variances:",mme.df.random)
    if mme.pedTrmVec!=0
        @printf("%-30s %20.3f\n","polygenic effect variances:",mme.df.polygenic)
    end
    if mme.M!=0
        @printf("%-30s %20.3f\n","marker effect variances:",mme.df.marker)
    end
    @printf("\n\n\n")
end
