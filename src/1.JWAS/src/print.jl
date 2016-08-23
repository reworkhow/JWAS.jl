function MCMCinfo(methods,chain_length,starting_value,printout_frequency,
    output_samples_frequency,missing_phenotypes,constraint,estimatePi,update_priors_frequency,mme)

    println("MCMC Information:")
    @printf("%-30s %20s\n","methods",methods)
#    @printf("%-20s %10s\n","seed",seed)
    @printf("%-30s %20s\n","chain_length",chain_length)
    @printf("%-30s %20s\n","estimatePi",estimatePi?"true":"false")
    @printf("%-30s %20s\n","constraint",constraint?"true":"false")
    @printf("%-30s %20s\n","missing_phenotypes",missing_phenotypes?"true":"false")
    @printf("%-30s %20s\n","starting_value",starting_value)
    @printf("%-30s %20d\n","output_samples_frequency",output_samples_frequency)
    @printf("%-30s %20d\n","printout_frequency",printout_frequency)
    @printf("%-30s %20d\n","update_priors_frequency",update_priors_frequency)

    @printf("\n%-30s\n","Degree of freedom for hyper-parameters:")
    @printf("%-30s %20.3f\n","residual variances:",mme.df.residual)
    @printf("%-30s %20.3f\n","iid random effect variances:",mme.df.random)
    @printf("%-30s %20.3f\n","polygenic effect variances:",mme.df.polygenic)
    @printf("%-30s %20.3f\n","marker effect variances:",mme.df.marker)
    @printf("\n\n\n")
end

"""
    showMME(mme::MME,df::DataFrame)

* Show left-hand side and right-hand side of mixed model equations (no markers).
"""
function showMME(mme::MME,df::DataFrame)
   if size(mme.mmeRhs)==()
     getMME(mme,df)
   end
   return [getNames(mme) mme.mmeLhs],[getNames(model) mme.mmeRhs]
end

function getNames(mme::MME)
    names = Array(AbstractString,0)
    for trm in mme.modelTerms
        for name in trm.names
            push!(names,trm.trmStr*" : "*name)
        end
    end
    return names
end

function getNames(trm::ModelTerm)
    names = Array(AbstractString,0)
    for name in trm.names
        push!(names,trm.trmStr*" : "*name)
    end
    return names
end
