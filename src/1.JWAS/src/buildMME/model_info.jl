################################################################################
#Print out model information
################################################################################
"""
    getinfo(model::MME)

* Print out model information.
"""
function getinfo(model::MME;data=false)
  printstyled("A Linear Mixed Model was build using model equations:\n\n",bold=true)
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
  #incidence matrix , #elements non-zero elements
  #"incomplete or complete",genomic data
end

################################################################################
#Print out MCMC information #Make all info a member of MME?
################################################################################
function getMCMCinfo(mme)
    MCMCinfo = mme.MCMCinfo
    printstyled("MCMC Information:\n\n",bold=true)

    @printf("%-30s %20s\n","methods",MCMCinfo.methods)
    if !(MCMCinfo.methods in ["conventional (no markers)"])
      @printf("%-30s %20s\n","genomic data",MCMCinfo.single_step_analysis ? "incomplete genomic data (single-step analysis)" : "complete genomic data (non-single-step analysis)")
    end
    @printf("%-30s %20s\n","chain_length",MCMCinfo.chain_length)
    @printf("%-30s %20s\n","burnin",MCMCinfo.burnin)
    if !(MCMCinfo.methods in ["conventional (no markers)", "GBLUP"])
      @printf("%-30s %20s\n","estimatePi",MCMCinfo.estimatePi ? "true" : "false")
    end
    @printf("%-30s %20s\n","estimateScale",MCMCinfo.estimateScale ? "true" : "false")
    @printf("%-30s %20s\n","starting_value",MCMCinfo.starting_value ? "true" : "false")
    @printf("%-30s %20d\n","printout_frequency",MCMCinfo.printout_frequency)
    @printf("%-30s %20d\n","output_samples_frequency",MCMCinfo.output_samples_frequency)
    @printf("%-30s %20s\n","constraint",MCMCinfo.constraint ? "true" : "false")
    @printf("%-30s %20s\n","missing_phenotypes",MCMCinfo.missing_phenotypes ? "true" : "false")
    @printf("%-30s %20d\n","update_priors_frequency",MCMCinfo.update_priors_frequency)
    @printf("%-30s %20s\n","seed",MCMCinfo.seed)

    printstyled("\nHyper-parameters Information:\n\n",bold=true)

    if mme.nModels==1
        for i in mme.rndTrmVec
            thisterm=split(i.term_array[1],':')[end]
            @printf("%-30s %20s\n","random effect variances ("*thisterm*"):",round.(inv(i.GiNew),digits=3))
        end
        @printf("%-30s %20.3f\n","residual variances:",mme.RNew)
        if mme.pedTrmVec!=0
            @printf("%-30s\n %50s\n","genetic variances (polygenic):",round.(inv(mme.GiNew),digits=3))
        end
        if !(MCMCinfo.methods in ["conventional (no markers)", "GBLUP"])
            if mme.M == 0
                error("Please add genotypes using add_genotypes().")
            end
            @printf("%-30s %20.3f\n","genetic variances (genomic):",mme.M.genetic_variance)
            @printf("%-30s %20.3f\n","marker effect variances:",mme.M.G)
            @printf("%-30s %20s\n","π",MCMCinfo.Pi)
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
        if !(MCMCinfo.methods in ["conventional (no markers)", "GBLUP"])
            @printf("%-30s\n","genetic variances (genomic):")
            if mme.M.genetic_variance != false
                Base.print_matrix(stdout,round.(mme.M.genetic_variance,digits=3))
            end
            println()
            @printf("%-30s\n","marker effect variances:")
            Base.print_matrix(stdout,round.(mme.M.G,digits=3))
            println()
            println("\nΠ: (Y(yes):included; N(no):excluded)\n")
            print(string.(mme.lhsVec))
            @printf("%20s\n","probability")
            for (i,j) in MCMCinfo.Pi
                i = replace(string.(i),"1.0"=>"Y","0.0"=>"N")
                print(i)
                @printf("%20s\n",j)
            end
        end
    end

    printstyled("\nDegree of freedom for hyper-parameters:\n\n",bold=true)
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
