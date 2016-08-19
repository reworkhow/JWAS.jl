function ssBayesC0(matrices::HybridMatrices,
                 geno::misc.Genotypes,fixed::misc.FixedMatrix,
                 ped::PedModule.Pedigree,
                 input::misc.InputParameters;outFreq=5000)

    y   = matrices.y.full
    X   = matrices.X.full
    W   = matrices.W.full
    Zn  = matrices.Z.n
    Ai_nn = matrices.Ai.nn

    current   = misc.Current(input,geno,fixed,y)
    current.fixed_effects = zeros(length(current.fixed_effects)+1) #add one for μ_g
    current.imputation_residual = zeros(matrices.num.pedn)

    output                 = misc.Output(input,geno,fixed)
    output.mean_fixed_effects = zeros(length(output.mean_fixed_effects)+1) #add one for μ_g
    output.mean_imputation_residual        = zeros(matrices.num.pedn)

    wGibbs = misc.GibbsMats(W)
    xGibbs = misc.GibbsMats(X)

    meanVare  = 0.0
    meanVara  = 0.0
    meanVarg  = 0.0

    println("running ",input.method," with a MCMC of length ",input.chainLength)

    for iter = 1:input.chainLength

      current.iter += 1
      # sample fixed effects
      misc.sample_fixed!(xGibbs,current,output)
      # sample marker effects
      misc.sample_random_ycorr!(wGibbs,current,output)
      # sample epsilon
      SSBR.sampleEpsi!(matrices,current,output)
      # sample residual vairance
      current.varResidual=misc.sample_variance(current.yCorr, matrices.num.y, input.nuRes, current.scaleRes)
      # sample marker vairance
      current.varEffect = misc.sample_variance(current.α, geno.nMarkers, input.dfEffectVar, current.scaleVar)
      # sample marker vairance
      current.varGenotypic = misc.sample_epsilon_variance(current.imputation_residual,
                                                     Ai_nn,
                                                     matrices.num.pedn,
                                                     input.nuGen,
                                                     current.scaleGen)
      #print out values to check convergence
      meanVare += (current.varResidual - meanVare)/current.iter
      meanVara += (current.varEffect - meanVara)/current.iter
      meanVarg += (current.varGenotypic - meanVarg)/current.iter

      #save monte carlo mean for all samples of variance components at each iteration
      output.resVar[iter] = meanVare
      output.genVar[iter] = meanVarg

      if (iter%outFreq ==0)
        @printf("Iteration %d with mean residual/marker effect/genetic(imputation) variance %6.3f/%6.3f/%6.3f.\n",
               iter, meanVare, meanVara, meanVarg)
      end
    end

    estimatedMarkerEffects = output.meanMarkerEffects

    mu_g = output.mean_fixed_effects[end]
    EBV = matrices.J.full*mu_g+matrices.M.full*estimatedMarkerEffects
    EBV[1:matrices.num.pedn,:] += output.mean_imputation_residual

    IDs=PedModule.getIDs(ped);
    EBV=DataFrame(ID=IDs,EBV=vec(EBV))

    return Results(matrices,output,EBV)
end
