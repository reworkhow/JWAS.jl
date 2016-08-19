function MCMC_BayesC(nIter,mme,df;
                     π         =0.0,
                     estimatePi=false,
                     sol       =false,
                     outFreq   =1000,
                     methods   ="conventional analyses",
                     output_marker_effects_frequency =0)

   #######################################################
   # Pre-Check
   #######################################################
    if mme.M!=0
      if π==0.0
        methods="BayesC0"
      end
      if methods=="BayesC0" && π != 0.0
        error("BayesC0 runs only with π=0
        Please remove argument: Pi")
      end #π for BayesC0 = 0.0
      if methods=="BayesC0" && estimatePi == true
        error("BayesC0 runs with estimatePi = false.
        Please remove argument: estimatePi")
      end #π for BayesC0 = 0.0
    end

    if size(mme.mmeRhs)==()
        getMME(mme,df)
    end

    if sol == false #no starting values
       sol = zeros(size(mme.mmeLhs,1))
    else #besure type is Float64
       sol = map(Float64,sol)
    end
    solMean = fill(0.0,size(mme.mmeLhs,1))

    #######################################################
    # PRIORS
    #######################################################
    #prior for residual variance
    vRes        = mme.RNew
    nuRes       = 4
    scaleRes    = vRes*(nuRes-2)/nuRes
    meanVare    = 0.0
    #priors for genetic variance (polygenic effects;A) e.g Animal+ Maternal
    if mme.ped != 0
       ν         =  4
       pedTrmVec = mme.pedTrmVec
       k         = size(pedTrmVec,1)  #2
       νG0       = ν + k
       G0        = inv(mme.GiNew)
       P         = G0*(νG0 - k - 1)
       S         = zeros(Float64,k,k)
       G0Mean    = zeros(Float64,k,k)
    end
    #priors for marker effect variance
    if mme.M != 0
      mGibbs      = GibbsMats(mme.M.genotypes)
      nObs        = mGibbs.nrows
      nMarkers    = mGibbs.ncols
      mArray      = mGibbs.xArray
      mpm         = mGibbs.xpx
      M           = mGibbs.X
      dfEffectVar = 4.0
      vEff        = mme.M.G
      scaleVar    = vEff*(dfEffectVar-2)/dfEffectVar   #scale factor for locus effects
      meanVara    = 0.0


      α           = zeros(nMarkers)#starting values for marker effeccts are zeros
      δ           = zeros(nMarkers)#inclusion indicator for marker effects (for BayesC)
      meanAlpha   = zeros(nMarkers)#vectors to save solutions for marker effects

      pi          = zeros(nIter)#vector to save π (for BayesC)
      mean_pi     = 0.0
    end

    #####################################################
    #  WORKING VECTORS (ycor)
    #####################################################
    ##starting values for location parameters(no marker) are sol
    #and adjust y for starting values
    ycorr       = vec(full(mme.ySparse)-mme.X*sol)

    #######################################################
    #  SET UP OUTPUT
    #######################################################
    #initiate vectors to save samples of MCMC
    initSampleArrays(mme,nIter)

    if output_marker_effects_frequency != 0  #write samples for marker effects to a txt file
      if mme.M != 0
        outfile   = open("MCMC samples for marker effects"*"_$(now()).txt","w")
        if mme.M.markerID[1]!="NA"
          writedlm(outfile,mme.M.markerID')
        end
      end
    end

    #######################################################
    # MCMC
    #######################################################
    @showprogress "running MCMC for "*methods*"..." for iter=1:nIter

        #####################################
        # 1.1. Non-Marker Location Parameters
        #####################################
        ycorr = ycorr + mme.X*sol
        rhs = mme.X'ycorr

        Gibbs(mme.mmeLhs,sol,rhs,vRes)

        ycorr = ycorr - mme.X*sol
        solMean += (sol - solMean)/iter

        #####################################
        # 1.2 Marker Effects
        #####################################
        if mme.M != 0
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
            pi[iter] = π
            mean_pi += (π-mean_pi)/iter
          end
        end

        ##########################################################################
        # 2.1 Genetic Covariance Matrix (Polygenic Effects)
        ##########################################################################
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
        ##########################################################################
        # 2.2 varainces for (iid) random effects;not required(empty)=>jump out
        ##########################################################################
        sampleVCs(mme,sol,iter)
        addLambdas(mme)

        ###############################################
        # 2.3 Residual Variamce
        ###############################################
        mme.ROld = mme.RNew
        vRes     = sample_variance(ycorr, length(ycorr), nuRes, scaleRes)
        mme.RNew = vRes
        meanVare += (vRes - meanVare)/iter
        mme.resVarSampleArray[iter] = vRes

        ###############################################
        # 2.4 Marker Effects Variamce
        ###############################################
        if mme.M!=0
          vEff  = sample_variance(α, nLoci, dfEffectVar, scaleVar)
          meanVara += (vEff - meanVara)/iter
        end

        #output samples for different effects
        outputSamples(mme,sol,iter)
        if output_marker_effects_frequency != 0  #write samples for marker effects to a txt file
          if mme.M != 0
            if iter%output_marker_effects_frequency==0
              writedlm(outfile,α')
            end
          end
        end

        if iter%outFreq==0
            println("\nPosterior means at iteration: ",iter)
            println("Residual variance: ",round(meanVare,3))
            if mme.ped !=0
              println("Polygenic effects covariance matrix \n",round(G0Mean,3))
            end
            if mme.M != 0
              println("Marker effects variance: ",round(meanVara,3))
              if estimatePi == true
                println("π: ", round(mean_pi,3))
              end
            end
        end
    end

    #######################################################
    # OUTPUT
    #######################################################
    if output_marker_effects_frequency != 0  #write samples for marker effects to a txt file
      if mme.M !=0
        close(outfile)
      end
    end

    output = Dict()
    output["Posterior Mean of Location Parameters"] = [getNames(mme) solMean]
    output["MCMC samples for residual variance"]    = mme.resVarSampleArray
    if mme.ped != 0
        output["MCMC samples for polygenic effects var-cov parameters"] = mme.genVarSampleArray
    end
    if mme.M != 0
        output["Posterior Mean of Marker Effects"] = meanAlpha
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
