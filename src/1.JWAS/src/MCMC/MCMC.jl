include("outputMCMCsamples.jl")
include("DRY.jl")
#include("MCMC_BayesB.jl")
include("MCMC_BayesC.jl")
include("MCMC_GBLUP.jl")
#include("MT_MCMC_BayesB.jl")
include("MT_MCMC_BayesC.jl")

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
                methods             = "conventional (no markers)",
                output_file         = "MCMC_samples",
                printout_frequency  = chain_length+1,
                printout_model_info = true,
                output_samples_frequency::Int64 = 0,
                update_priors_frequency::Int64=0)

    ############################################################################
    # Pre-Check
    ############################################################################
    starting_value = pre_check(mme,df,starting_value)

    if methods!="conventional (no markers)"
        #set up Pi
        if mme.M != 0 && mme.nModels !=1 && Pi==0.0
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
        if mme.M != 0 && mme.M.G_is_marker_variance==false && methods!="GBLUP"
            genetic2marker(mme.M,Pi)
            println()
            if mme.nModels != 1
              if !isposdef(mme.M.G) #also work for scalar
                error("Marker effects covariance matrix is not postive definite! Please modify the argument: Pi.")
              end
              println("The prior for marker effects covariance matrix is calculated ")
              println("from genetic covariance matrix and Π. The marker effects ")
              println("covariance matrix is: \n")
              Base.print_matrix(STDOUT,round.(mme.M.G,6))
            else
              if !isposdef(mme.M.G) #positive scalar (>0)
                error("Marker effects variance is negative!")
              end
              println("The prior for marker effects variance is calculated from ")
              println("the genetic variance and π. The marker effects variance ")
              println("is: ",round.(mme.M.G,6))
            end
        elseif methods=="GBLUP" && mme.M.G_is_marker_variance==true
            error("Please provide genetic variance for GBLUP analysis")
        end
        println("\n\n")
    end

    #set up starting values for location parameters (not markers)
    have_starting_value=false
    if starting_value != false
      starting_value=vec(starting_value)
      have_starting_value=true
    end

    #printout basic MCMC information
    if printout_MCMCinfo == true
      getinfo(mme)
      MCMCinfo(methods,Pi,chain_length,burnin,have_starting_value,printout_frequency,
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
        # elseif methods =="BayesB"
        #     res=MCMC_BayesB(chain_length,mme,df,
        #                     burnin                   = burnin,
        #                     π                        = Pi,
        #                     estimatePi               = estimatePi,
        #                     sol                      = starting_value,
        #                     outFreq                  = printout_frequency,
        #                     output_samples_frequency = output_samples_frequency,
        #                     output_file              = output_file)
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

    if !(methods in ["conventional (no markers)", "GBLUP"])
      @printf("\n%-30s\n\n","Hyper-parameters Information:")
      if issubtype(typeof(Pi),Number)
          @printf("%-30s %20s\n","π",Pi)
      else
          println("Π")
          @printf("%-20s %12s\n","combinations","probability")

          for (i,j) in Pi
              @printf("%-20s %12s\n",i,j)
          end
          println()
      end
      # @printf("%-30s %20.3f\n","genetic variances (genomic):",mme.df.residual)
      # @printf("%-30s %20.3f\n","genetic variances (polygenic):",mme.df.residual)
      # @printf("%-30s %20.3f\n","marker effect variances:",mme.df.residual)
      # @printf("%-30s %20.3f\n","residual variances:",mme.df.residual)
      # @printf("%-30s %20.3f\n","random effect variances (x1):",mme.df.random)
      # @printf("%-30s %20.3f\n","random effect variances (x2):",mme.df.random)
    end

    @printf("\n%-30s\n\n","Degree of freedom for hyper-parameters:")
    @printf("%-30s %20.3f\n","residual variances:",mme.df.residual)
    @printf("%-30s %20.3f\n","iid random effect variances:",mme.df.random)
    if mme.ped!=0
        @printf("%-30s %20.3f\n","polygenic effect variances:",mme.df.polygenic)
    end
    if mme.M!=0
        @printf("%-30s %20.3f\n","marker effect variances:",mme.df.marker)
    end
    @printf("\n\n\n")
end
