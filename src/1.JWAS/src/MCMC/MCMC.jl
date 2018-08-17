include("outputMCMCsamples.jl")
include("DRY.jl")
include("MCMC_BayesC.jl")
include("MCMC_GBLUP.jl")
include("MT_MCMC_BayesC.jl")
include("../SSBR/SSBR.jl")
include("output.jl")

"""
    runMCMC(mme,df;Pi=0.0,estimatePi=false,chain_length=1000,burnin = 0,starting_value=false,printout_frequency=100,missing_phenotypes=false,constraint=false,methods="conventional (no markers)",output_samples_frequency::Int64 = 0,printout_model_info=true)

Run MCMC (marker information included or not) with sampling of variance components.

* available **methods** include "conventional (no markers)", "RR-BLUP", "BayesB", "BayesC".
* save MCMC samples every **output_samples_frequency** iterations to files **output_file** defaulting to `MCMC_samples`.
* the first **burnin** iterations are discarded at the beginning of an MCMC run
* **Pi** for single-trait analyses is a number; **Pi** for multi-trait analyses is a dictionary such as `Pi=Dict([1.0; 1.0]=>0.7,[1.0; 0.0]=>0.1,[0.0; 1.0]=>0.1,[0.0; 0.0]=>0.1)`,
    * if Pi (Π) is not provided in multi-trait analysis, it will be generated assuming all markers have effects on all traits.
* **starting_value** can be provided as a vector for all location parameteres except marker effects.
* print out the monte carlo mean in REPL with **printout_frequency**
* **constraint**=true if constrain residual covariances between traits to be zeros.
"""
function runMCMC(mme,df;
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
                #Genomic Prediction
                outputEBV           = false,
                #parameters for single-step analysis
                single_step_analysis= false,
                pedigree            = false,
                #output
                output_genetic_variance = false,
                output_residual_variance= false,
                output_heritability     = false)

    ############################################################################
    # Pre-Check
    ############################################################################
    if methods =="conventional (no markers)" && mme.M!=0
        error("Conventional analysis runs without genotypes!")
    end

    #users need to provide high-quality pedigree file
    if mme.ped!=0 || pedigree!=false
        check_pedigree(mme,df,pedigree)
    end


    #Genotyped individuals are usaully not many, and are used in GWAS (complete
    #and incomplete), thus those IDs are default output_ID if not provided
    if outputEBV == false
        mme.output_ID = 0
    elseif mme.output_ID == 0 && mme.M != 0
        mme.output_ID = mme.M.obsID
    end

    #single-step
    #impute genotypes for non-genotyped individuals
    #add imputation error and J as varaibles in data
    if single_step_analysis==true
        SSBRrun(mme,pedigree,df)
    end
    #make mixed model equations for non-marker parts
    #assign IDs for observations
    starting_value,df = pre_check(mme,df,starting_value)


    if mme.M!=0
        #align genotypes with phenotypes IDs and align genotypes with output IDs
        align_genotypes(mme)
    end
    if mme.M!=0
        #set up Pi
        if mme.nModels !=1 && Pi==0.0
            info("Pi (Π) is not provided.","\n")
            info("Pi (Π) is generated assuming all markers have effects on all traits.","\n")
            mykey=Array{Float64}(0)
            ntraits=mme.nModels
            Pi=Dict{Array{Float64,1},Float64}()
            for i in [ bin(n,ntraits) for n in 0:2^ntraits-1 ]
              Pi[float(split(i,""))]=0.0
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
              println("The prior for marker effects covariance matrix is calculated from ")
              println("genetic covariance matrix and Π. The prior for the marker effects ")
              println("covariance matrix is: \n")
              Base.print_matrix(STDOUT,round.(mme.M.G,6))
            else
              if !isposdef(mme.M.G) #positive scalar (>0)
                error("Marker effects variance is negative!")
              end
              println("The prior for marker effects variance is calculated from ")
              println("the genetic variance and π. The prior for the marker effects variance ")
              println("is: ",round.(mme.M.G,6))
            end
        elseif methods=="GBLUP" && mme.M.G_is_marker_variance==true
            error("Please provide genetic variance for GBLUP analysis")
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
                            output_file              = output_file)
        elseif methods =="GBLUP"
            res=MCMC_GBLUP(chain_length,mme,df;
                            burnin                   = burnin,
                            sol                      = starting_value,
                            outFreq                  = printout_frequency,
                            output_samples_frequency = output_samples_frequency,
                            output_file              = output_file)
        else
            error("No options!!!")
        end
    elseif mme.nModels > 1
        if Pi != 0.0 && round(sum(values(Pi)),2)!=1.0
          error("Summation of probabilities of Pi is not equal to one.")
        end
        if methods in ["BayesC","BayesCC","BayesB","RR-BLUP","conventional (no markers)"]
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
  res
end

################################################################################
#Print out MCMC information
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
      @printf("%-30s %20s\n","estimatePi",estimatePi?"true":"false")
    end
    @printf("%-30s %20s\n","starting_value",starting_value?"true":"false")
    @printf("%-30s %20d\n","printout_frequency",printout_frequency)
    @printf("%-30s %20d\n","output_samples_frequency",output_samples_frequency)

    @printf("%-30s %20s\n","constraint",constraint?"true":"false")
    @printf("%-30s %20s\n","missing_phenotypes",missing_phenotypes?"true":"false")
    @printf("%-30s %20d\n","update_priors_frequency",update_priors_frequency)


    @printf("\n%-30s\n\n","Hyper-parameters Information:")

    if mme.nModels==1
        for i in mme.rndTrmVec
            thisterm=split(i.term_array[1],':')[end]
            @printf("%-30s %20s\n","random effect variances ("*thisterm*"):",round.(inv(i.GiNew),3))
        end
        @printf("%-30s %20.3f\n","residual variances:",mme.RNew)
        if mme.pedTrmVec!=0
            @printf("%-30s\n %50s\n","genetic variances (polygenic):",round.(inv(mme.GiNew),3))
        end
        if !(methods in ["conventional (no markers)", "GBLUP"])
            if mme.M == 0
                error("Please add genotypes using add_genotypes().")
            end
            @printf("%-30s %20.3f\n","genetic variances (genomic):",mme.M.G)
            @printf("%-30s %20.3f\n","marker effect variances:",mme.M.G)
            @printf("%-30s %20s\n","π",Pi)
        end
    else
        for i in mme.rndTrmVec
            thisterm=split(i.term_array[1],':')[end]
            @printf("%-30s\n","random effect variances ("*thisterm*"):")
            Base.print_matrix(STDOUT,round.(inv(i.GiNew),3))
            println()
        end
        @printf("%-30s\n","residual variances:")
        Base.print_matrix(STDOUT,round.(mme.R,3))
        println()
        if mme.pedTrmVec!=0
            @printf("%-30s\n","genetic variances (polygenic):")
            Base.print_matrix(STDOUT,round.(inv(mme.Gi),3))
            println()
        end
        if !(methods in ["conventional (no markers)", "GBLUP"])
            @printf("%-30s\n","genetic variances (genomic):")
            Base.print_matrix(STDOUT,round.(mme.M.G,3))
            println()
            @printf("%-30s\n","marker effect variances:")
            Base.print_matrix(STDOUT,round.(mme.M.G,3))
            println()
            println("Π:")
            @printf("%-20s %12s\n","combinations","probability")
            for (i,j) in Pi
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
