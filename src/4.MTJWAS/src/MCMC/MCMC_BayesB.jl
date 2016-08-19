function MCMC_BayesB(nIter,mme,df,π;
                     estimatePi=false,
                     sol       =false,
                     outFreq   =100,
                     output_marker_effects_frequency =0)
    if size(mme.mmeRhs)==()
        getMME(mme,df)
    end

    if sol == false
        sol=zeros(size(mme.mmeLhs,1))
    end #starting value for sol can be provided


    #######################################################
    # PRIORS
    #######################################################
    #prior for residual variance
    vRes        = mme.RNew
    nuRes       = 4
    scaleRes    = vRes*(nuRes-2)/nuRes
    #priors for marker effect variance
    mGibbs      = GibbsMats(mme.M.genotypes)
    nObs        = mGibbs.nrows
    nMarkers    = mGibbs.ncols
    mArray      = mGibbs.xArray
    mpm         = mGibbs.xpx
    M           = mGibbs.X
    dfEffectVar = 4.0
    vEff        = mme.M.G
    locusEffectVar = fill(vEff,nMarkers)
    scaleVar       = vEff*(dfEffectVar-2)/dfEffectVar  #scale factor for locus effects
    #priors for genetic variance (polygenic effects;A) e.g Animal+ Maternal
    if mme.ped != 0
       ν = 4
       pedTrmVec = mme.pedTrmVec
       k         = size(pedTrmVec,1)  #2
       νG0       = ν + k
       G0        = inv(mme.GiNew)
       P         = G0*(νG0 - k - 1)
       S         = zeros(Float64,k,k)
       G0Mean    = zeros(Float64,k,k)
    end

    #####################################################
    # WORKING VECTORS (ycor, saving values)
    #####################################################
    #vectors to save solutions for conventional MME part
    p           = size(mme.mmeRhs,1)
    solMean     = zeros(p)
    #vectors to save solutions for marker effects
    α           = zeros(nMarkers)       #starting values for marker effeccts are zeros
    δ           = zeros(nMarkers)       # inclusion indicator for marker effects
    u           = zeros(nMarkers)       # inclusion indicator for marker effects
    meanu       = zeros(nMarkers)
    if output_marker_effects_frequency != 0  #write samples for marker effects to a txt file
      outfile   = open("MCMC samples for marker effects"*"_$(now()).txt","w")
      if mme.M.markerID[1]!="NA"
        writedlm(outfile,mme.M.markerID')
      end
    end

    #initiate vectors to save samples of MCMC
    initSampleArrays(mme,nIter)
    #variables to save variance for marker effects or residual
    meanVare    = 0.0
    meanVara    = 0.0
    #vector to save π
    pi          = zeros(nIter)
    #adjust y for strating values
    ycorr       = vec(full(mme.ySparse)-mme.X*sol)   #starting values for location parameters(no marker) are sol

    #######################################################
    # MCMC
    #######################################################
    @showprogress "running MCMC " for iter=1:nIter

        #sample non-marker part
        ycorr = ycorr + mme.X*sol
        rhs = mme.X'ycorr

        Gibbs(mme.mmeLhs,sol,rhs,vRes)

        ycorr = ycorr - mme.X*sol
        solMean += (sol - solMean)/iter

        #sample marker effects
        nLoci = sampleEffectsBayesB!(mArray,mpm,ycorr,u,α,δ,vRes,locusEffectVar,π)
        meanu += (u - meanu)/iter

        #sample variances for polygenic effects
        if mme.ped != 0
            for (i,trmi) = enumerate(pedTrmVec)
                pedTrmi   = mme.modelTermDict[trmi]
                startPosi = pedTrmi.startPos
                endPosi   = startPosi + pedTrmi.nLevels - 1
                for (j,trmj) = enumerate(pedTrmVec)
                    pedTrmj   = mme.modelTermDict[trmj]
                    startPosj = pedTrmj.startPos
                    endPosj   = startPosj + pedTrmj.nLevels - 1
                    S[i,j]    = (sol[startPosi:endPosi]'*mme.Ai*sol[startPosj:endPosj])[1,1]
                end
            end

            pedTrm1 = mme.modelTermDict[pedTrmVec[1]]
            q  = pedTrm1.nLevels
            G0 = rand(InverseWishart(νG0 + q, P + S)) #ν+q?

            mme.GiOld = copy(mme.GiNew)
            mme.GiNew = inv(G0)

            G0Mean  += (G0  - G0Mean )/iter
            mme.genVarSampleArray[iter,:] = vec(G0)

            addA(mme) #add Ainverse*lambda
        end

        #sample varainces for (iid) random effects;not required(empty)=>jump out
        sampleVCs(mme,sol,iter)
        addLambdas(mme)

        #sample residual variance
        mme.ROld = mme.RNew
        vRes     = sample_variance(ycorr, nObs, nuRes, scaleRes)
        mme.RNew = vRes
        meanVare += (vRes - meanVare)/iter
        mme.resVarSampleArray[iter] = vRes

        # sample locus effect variance
        for j=1:nMarkers
            locusEffectVar[j] = sample_variance(α[j],1,dfEffectVar, scaleVar)
        end

        ##sample Pi
        #if estimatePi == true
        #  π = samplePi(nLoci, nMarkers)
        #  pi[iter] = π
        #end

        #output samples for different effects
        outputSamples(mme,sol,iter)
        if output_marker_effects_frequency != 0  #write samples for marker effects to a txt file
          if iter%output_marker_effects_frequency==0
            writedlm(outfile,u')
          end
        end

        if iter%outFreq==0
            println("\nPosterior means at iteration: ",iter)
            println("Residual variance: ",round(meanVare,3))
            if estimatePi == true
              println("π: ", round(mean_pi,3))
            end
        end

    end

    #######################################################
    # OUTPUT
    #######################################################
    if output_marker_effects_frequency != 0  #write samples for marker effects to a txt file
      close(outfile)
    end

    output = Dict()
    output["Posterior Mean of Location Parameters"] = [getNames(mme) solMean]
    output["MCMC samples for residual variance"]    = mme.resVarSampleArray
    if mme.ped != 0
      output["MCMC samples for polygenic effects var-cov parameters"] = mme.genVarSampleArray
    end
    if mme.M != 0
        output["Posterior Mean of Marker Effects"] = meanu
    end
    if estimatePi == true
        output["MCMC samples for: π"] = pi
    end
    for i in  mme.outputSamplesVec
        trmi   = i.term
        trmStr = trmi.trmStr
        output["MCMC samples for: "*trmStr] = i.sampleArray
    end
    for i in  mme.rndTrmVec
        trmi   = i.term
        trmStr = trmi.trmStr
        output["MCMC samples for: variance of "*trmStr] = i.sampleArray
    end

    return output
end
