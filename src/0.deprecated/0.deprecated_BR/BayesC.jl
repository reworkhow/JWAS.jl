function BayesC(y::Array{Float64,1},
                geno::Genotypes,
                fixed::FixedMatrix,
                input::QTL.InputParameters;outFreq=5000)

    y   = y
    X   = fixed.C
    W   = geno.genotypes

    current   = Current(input,geno,fixed,y)
    output    = QTL.Output(input,geno,fixed)

    wGibbs = GibbsMats(W)
    xGibbs = GibbsMats(X)

    meanVare  = 0.0
    meanVara  = 0.0
    meanPi    = 0.0

    println("running ",input.method," with a MCMC of length ",input.chainLength)

    for iter = 1:input.chainLength

      current.iter += 1
      # sample fixed effects
      sample_fixed!(xGibbs,current,output)
      # sample marker effects
      sampleEffectsBayesC!(wGibbs,current,output)
      # sample residual vairance
      current.varResidual=sample_variance(current.yCorr, geno.nObs, input.nuRes, current.scaleRes)
      # sample marker vairance
      current.varEffect = sample_variance(current.α, current.nLoci, input.dfEffectVar, current.scaleVar)

      if (input.estimatePi == true)
            current.π = samplePi(current.nLoci, geno.nMarkers)
      end

      #print out values to check convergence
      meanVare += (current.varResidual - meanVare)/current.iter
      meanVara += (current.varEffect - meanVara)/current.iter
      meanPi   += (current.π - meanPi)/current.iter

       if (iter%outFreq ==0)
        @printf("Iteration %d with %d loci included in the model, mean residual/marker effect %6.3f/%6.3f with mean π= %6.3f.\n",
               iter,current.nLoci, meanVare, meanVara, meanPi)
       end
    end

    estimatedMarkerEffects = output.meanMarkerEffects

    return estimatedMarkerEffects
end
