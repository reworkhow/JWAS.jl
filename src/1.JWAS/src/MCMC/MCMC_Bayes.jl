function MCMC_BayesC(nIter,mme,df;
                     π         =0.0,
                     estimatePi=false,
                     sol       =false,
                     outFreq   =1000,
                     methods   ="conventional analyses",
                     output_samples_frequency=0)

    #######################################################
    # Pre-Check
    #######################################################
    pre_check(mme,π,sol)

    if mme.M!=0
      if π==0.0
        methods="BayesC0"
      end
      if methods=="BayesC0" && π != 0.0
        error("BayesC0 runs with π=0.")
      end
      if methods=="BayesC0" && estimatePi == true
        error("BayesC0 runs with estimatePi = false.")
      end
    end

    #######################################################
    # PRIORS
    #######################################################
    #prior for residual variance
    vRes        = mme.RNew
    nuRes       = mme.df.residual
    scaleRes    = vRes*(nuRes-2)/nuRes
    meanVare    = 0.0
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
    #priors for marker effect variance
    if mme.M != 0
      mGibbs      = GibbsMats(mme.M.genotypes)
      nObs        = mGibbs.nrows
      nMarkers    = mGibbs.ncols
      mArray      = mGibbs.xArray
      mpm         = mGibbs.xpx
      M           = mGibbs.X
      dfEffectVar = mme.df.marker
      vEff        = mme.M.G
      scaleVar    = vEff*(dfEffectVar-2)/dfEffectVar   #scale factor for locus effects
      meanVara    = 0.0

      α           = zeros(nMarkers)#starting values for marker effeccts are zeros
      δ           = zeros(nMarkers)#inclusion indicator for marker effects (for BayesC)
      meanAlpha   = zeros(nMarkers)#vectors to save solutions for marker effects
      mean_pi     = 0.0
    end

    #####################################################
    #  WORKING VECTORS (ycor)
    #####################################################
    ##starting values for location parameters(no marker) are sol
    #and adjust y for starting values
    ycorr       = vec(full(mme.ySparse)-mme.X*sol)
    solMean     = fill(0.0,size(mme.mmeLhs,1))

    #######################################################
    #  SET UP OUTPUT
    #######################################################
    if output_samples_frequency != 0
      #initialize arrays to save MCMC samples
      num_samples = Int(floor(nIter/output_samples_frequency))
      init_sample_arrays(mme,num_samples)

      if mme.M != 0 #write samples for marker effects to a txt file
        file_count = 1
        file_name="MCMC_samples_for_marker_effects.txt"
        while isfile(file_name)
          file_name="MCMC_samples_for_marker_effects"*"_$(file_count)"*".txt"
          file_count += 1
        end
        outfile=open(file_name,"w")

        if mme.M.markerID[1]!="NA"
            writedlm(outfile,mme.M.markerID')
        end
        pi = zeros(num_samples)#vector to save π (for BayesC)
      end
      out_i = 1
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
            G0 = rand(InverseWishart(νG0 + q, round(P + S,7))) #ν+q?

            mme.GiOld = copy(mme.GiNew)
            mme.GiNew = inv(G0)

            G0Mean  += (G0  - G0Mean )/iter

            addA(mme) #add Ainverse*lambda
        end
        ##########################################################################
        # 2.2 varainces for (iid) random effects;not required(empty)=>jump out
        ##########################################################################
        sampleVCs(mme,sol)
        addLambdas(mme)

        ###############################################
        # 2.3 Residual Variamce
        ###############################################
        mme.ROld = mme.RNew
        vRes     = sample_variance(ycorr, length(ycorr), nuRes, scaleRes)
        mme.RNew = vRes
        meanVare += (vRes - meanVare)/iter

        ###############################################
        # 2.4 Marker Effects Variamce
        ###############################################
        if mme.M!=0
          vEff  = sample_variance(α, nLoci, dfEffectVar, scaleVar)
          meanVara += (vEff - meanVara)/iter
        end

        ###############################################
        # OUTPUT
        ###############################################
        if output_samples_frequency != 0
          if iter%output_samples_frequency==0
            outputSamples(mme,sol,out_i)
            mme.samples4R[:,out_i]=vRes
            if mme.ped != 0
              mme.samples4G[:,out_i]=vec(G0)
            end
            if mme.M != 0
                writedlm(outfile,α')
                pi[out_i] = π
            end
            out_i +=1
          end
        end

        if iter%outFreq==0
            println("\nPosterior means at iteration: ",iter)
            println("Residual variance: ",round(meanVare,6))
            if mme.ped !=0
              println("Polygenic effects covariance matrix \n",round(G0Mean,3))
            end
            if mme.M != 0
              println("Marker effects variance: ",round(meanVara,6))
              if estimatePi == true
                println("π: ", round(mean_pi,3))
              end
            end
        end
    end

    #######################################################
    # After MCMC
    #######################################################
    if output_samples_frequency != 0 && mme.M !=0
        close(outfile)
    end

    output = Dict()
    output["Posterior mean of location parameters"] = [getNames(mme) solMean]
    if output_samples_frequency != 0
        output["MCMC samples for residual variance"]    = mme.samples4R
        if mme.ped != 0
            output["MCMC samples for polygenic effects var-cov parameters"] = mme.samples4G
        end
        for i in  mme.outputSamplesVec
            trmi   = i.term
            trmStr = trmi.trmStr
            output["MCMC samples for: "*trmStr] = [getNames(trmi) i.sampleArray]
        end
        for i in  mme.rndTrmVec
            trmi   = i.term
            trmStr = trmi.trmStr
            output["MCMC samples for: variance of "*trmStr] = i.sampleArray
        end
    end

    if mme.M != 0
        if mme.M.markerID[1]!="NA"
            markerout=[mme.M.markerID meanAlpha]
        else
            markerout= meanAlpha
        end

        output["Posterior mean of marker effects"] = markerout
        if estimatePi == true
            output["Posterior mean of Pi"] = mean_pi
            if  output_samples_frequency != 0
                output["MCMC samples for: π"] = pi
            end
        end
    end
    return output
end
