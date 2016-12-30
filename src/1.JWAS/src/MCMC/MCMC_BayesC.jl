function MCMC_BayesC(nIter,mme,df;
                     π         =0.0,
                     estimatePi=false,
                     sol       =false,
                     outFreq   =1000,
                     methods   ="conventional analyses",
                     output_samples_frequency=0)

    ############################################################################
    # Pre-Check
    ############################################################################
    #starting values for location parameters(no marker) are sol
    sol,solMean = pre_check(mme,sol)

    if π==0.0
      methods="BayesC0"
    end
    if methods=="BayesC0" && π != 0.0
      error("BayesC0 runs with π=0.")
    end
    if methods=="BayesC0" && estimatePi == true
      error("BayesC0 runs with estimatePi = false.")
    end
    ############################################################################
    # PRIORS
    ############################################################################
    #prior for residual variance
    vRes        = mme.RNew
    nuRes       = mme.df.residual
    scaleRes    = vRes*(nuRes-2)/nuRes
    meanVare    = 0.0

    #priors for marker effect variance
    mGibbs      = GibbsMats(mme.M.genotypes)
    nObs,nMarkers,mArray,mpm,M = mGibbs.nrows,mGibbs.ncols,mGibbs.xArray,mGibbs.xpx,mGibbs.X
    dfEffectVar = mme.df.marker
    vEff        = mme.M.G
    scaleVar    = vEff*(dfEffectVar-2)/dfEffectVar #scale factor for locus effects
    meanVara    = 0.0 #variable to save variance for marker effect
    #vectors to save solutions for marker effects
    α           = zeros(nMarkers)#starting values for marker effeccts are zeros
    δ           = zeros(nMarkers)#inclusion indicator for marker effects
    meanAlpha   = zeros(nMarkers)#vectors to save solutions for marker effects
    mean_pi     = 0.0

    #priors for genetic variance (polygenic effects;A) e.g Animal+ Maternal
    if mme.ped != 0
       ν         = mme.df.polygenic
       pedTrmVec = mme.pedTrmVec
       k         = size(pedTrmVec,1)  #2
       νG0       = ν + k
       G0        = inv(mme.GiNew)
       P         = G0*(νG0 - k - 1)
       S         = zeros(Float64,k,k)
       G0Mean    = zeros(Float64,k,k)
    end
    ############################################################################
    #  WORKING VECTORS (ycor)
    ############################################################################
    #adjust y for starting values
    ycorr       = vec(full(mme.ySparse)-mme.X*sol)

    ############################################################################
    #  SET UP OUTPUT MCMC samples
    ############################################################################
    if output_samples_frequency != 0
      out_i,outfile,pi=output_MCMC_samples_setup(mme,nIter,output_samples_frequency)
    end

    ############################################################################
    # MCMC
    ############################################################################
    @showprogress "running MCMC for "*methods*"..." for iter=1:nIter

        ########################################################################
        # 1.1. Non-Marker Location Parameters
        ########################################################################
        ycorr = ycorr + mme.X*sol
        rhs = mme.X'ycorr

        Gibbs(mme.mmeLhs,sol,rhs,vRes)

        ycorr = ycorr - mme.X*sol
        solMean += (sol - solMean)/iter

        ########################################################################
        # 1.2 Marker Effects
        ########################################################################
        if methods=="BayesC"
          nLoci = sampleEffectsBayesC!(mArray, mpm, ycorr, α, δ,vRes, vEff, π)
        elseif methods=="BayesC0"
          sampleEffectsBayesC0!(mArray,mpm,ycorr,α,vRes,vEff)
          nLoci = nMarkers
        end
        meanAlpha += (α - meanAlpha)/iter

        #sample Pi
        if estimatePi == true
          π = samplePi(nLoci, nMarkers)
          mean_pi += (π-mean_pi)/iter
        end

        ########################################################################
        # 2.1 Genetic Covariance Matrix (Polygenic Effects) (variance.jl)
        ########################################################################
        if mme.ped != 0
          G0=sample_variance_pedigree(mme,pedTrmVec,sol,P,S,νG0)
          G0Mean  += (G0  - G0Mean )/iter
        end
        ########################################################################
        # 2.2 varainces for (iid) random effects;not required(empty)=>jump out
        ########################################################################
        sampleVCs(mme,sol)
        addLambdas(mme)
        ########################################################################
        # 2.3 Residual Variance
        ########################################################################
        mme.ROld = mme.RNew
        vRes     = sample_variance(ycorr, length(ycorr), nuRes, scaleRes)
        mme.RNew = vRes
        meanVare += (vRes - meanVare)/iter
        ########################################################################
        # 2.4 Marker Effects Variance
        ########################################################################
        if mme.M!=0
          vEff  = sample_variance(α, nLoci, dfEffectVar, scaleVar)
          meanVara += (vEff - meanVara)/iter
        end

        ########################################################################
        # 3.1 Save MCMC samples
        ########################################################################
        if output_samples_freqcy != 0 && iter%output_samples_frequency==0
            out_i=output_MCMC_samples(mme,out_i,pi,sol,α,vRes,G0,outfile,estimatePi)
        end
        ########################################################################
        # 3.2 Printout
        ########################################################################
        if iter%outFreq==0
            println("\nPosterior means at iteration: ",iter)
            println("Residual variance: ",round(meanVare,6))
            if mme.ped !=0
              println("Polygenic effects covariance matrix \n",round(G0Mean,3))
            end
            println("Marker effects variance: ",round(meanVara,6))
            if estimatePi == true
              println("π: ", round(mean_pi,3))
            end
        end
    end

    #######################################################
    # After MCMC
    #######################################################
    if output_samples_frequency != 0
        close(outfile)
    end

    output=output_result(mme,solMean,output_samples_frequency,estimatePi)
    return output
end
