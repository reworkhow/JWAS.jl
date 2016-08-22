function MCMCinfo(methods,chain_length,starting_value,printout_frequency,
    output_samples_frequency,missing_phenotypes,constraint,estimatePi,mme)

    println("MCMC Information:")
    @printf("%-30s %20s\n","methods",methods)
#    @printf("%-20s %10s\n","seed",seed)
    @printf("%-30s %20s\n","chain length",chain_length)
    @printf("%-30s %20s\n","estimatePi",estimatePi?"true":"false")
    @printf("%-30s %20s\n","constraint",constraint?"true":"false")
    @printf("%-30s %20s\n","missing_phenotypes",missing_phenotypes?"true":"false")
    @printf("%-30s %20s\n","starting value",starting_value)
    @printf("%-30s %20d\n","output sample frequency",output_samples_frequency)
    @printf("%-30s %20d\n","printout frequency",printout_frequency)

    @printf("\n%-30s\n","Degree of freedom for hyper-parametes:")
    @printf("%-30s %20.3f\n","residual variances:",mme.df.residual)
    @printf("%-30s %20.3f\n","iid random effect variances:",mme.df.random)
    @printf("%-30s %20.3f\n","polygenic effect variances:",mme.df.polygenic)
    @printf("%-30s %20.3f\n","marker effect variances:",mme.df.marker)
    @printf("\n\n\n")
end
