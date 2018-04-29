function MCMC_BayesB(nIter,mme,df,π;
                     burnin    =0,
                     sol       =false,
                     outFreq   =100,
                     output_samples_frequency =0,
                     MCMC_marker_effects_file="MCMC_samples_for_marker_effects.txt")

    ############################################################################
    # Pre-Check
    ############################################################################
    #starting values for location parameters(no marker) are sol
    sol,solMean = pre_check(mme,df,sol)

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
    locusEffectVar = fill(vEff,nMarkers)
    scaleVar       = vEff*(dfEffectVar-2)/dfEffectVar  #scale factor for locus effects
    meanVara       = 0.0 #variable to save variance for marker effect
    #vectors to save solutions for marker effects
    α           = zeros(nMarkers) #starting values for partial marker effeccts are zeros
    δ           = zeros(nMarkers) #inclusion indicator for marker effects
    u           = zeros(nMarkers) #marker effects
    meanu       = zeros(nMarkers) #vectors to save solutions for marker effects

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
    # WORKING VECTORS (ycor, saving values)
    ############################################################################
    #adjust y for starting values
    ycorr       = vec(full(mme.ySparse)-mme.X*sol)

    ############################################################################
    #  SET UP OUTPUT MCMC samples
    ############################################################################
    if output_samples_frequency != 0
      out_i,outfile,sample4π=output_MCMC_samples_setup(mme,nIter-burnin,output_samples_frequency,MCMC_marker_effects_file=MCMC_marker_effects_file)
    end #sample4π is not used in MME type since π is BayesC-specific

    #######################################################
    # MCMC
    #######################################################
    @showprogress "running MCMC for BayesB ..." for iter=1:nIter

        ########################################################################
        # 1.1. Non-Marker Location Parameters
        ########################################################################
        ycorr = ycorr + mme.X*sol
        rhs = mme.X'ycorr

        Gibbs(mme.mmeLhs,sol,rhs,vRes)

        ycorr = ycorr - mme.X*sol
        if iter > burnin
            solMean += (sol - solMean)/(iter-burnin)
        end
        ########################################################################
        # 1.2 Marker Effects
        ########################################################################
        nLoci = sampleEffectsBayesB!(mArray,mpm,ycorr,u,α,δ,vRes,locusEffectVar,π)
        if iter > burnin
            meanu += (u - meanu)/(iter-burnin)
        end
        ########################################################################
        # 2.1 Genetic Covariance Matrix (Polygenic Effects) (variance.jl)
        ########################################################################
        if mme.ped != 0
          G0=sample_variance_pedigree(mme,pedTrmVec,sol,P,S,νG0)
          if iter > burnin
              G0Mean  += (G0  - G0Mean )/(iter-burnin)
          end
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
        vRes     = sample_variance(ycorr, nObs, nuRes, scaleRes)
        mme.RNew = vRes
        if iter > burnin
            meanVare += (vRes - meanVare)/(iter-burnin)
        end
        ###############################################
        # 2.4 Marker Effects Variance
        ###############################################
        for j=1:nMarkers
            locusEffectVar[j] = sample_variance(α[j],1,dfEffectVar, scaleVar)
        end

        ########################################################################
        # 3.1 Save MCMC samples
        ########################################################################
        if output_samples_frequency != 0 && iter%output_samples_frequency==0 && iter>burnin
          out_i=output_MCMC_samples(mme,out_i,sol,vRes,(mme.ped!=0?G0:false),0.0,u,locusEffectVar,false,outfile,false)
        end
        ########################################################################
        # 3.2 Printout
        ########################################################################
        if iter%outFreq==0 && iter>burnin
            println("\nPosterior means at iteration: ",iter)
            println("Residual variance: ",round(meanVare,6))
            if mme.ped !=0
              println("Polygenic effects covariance matrix \n",round.(G0Mean,3))
            end
        end
    end

    ############################################################################
    # After MCMC
    ############################################################################
    if output_samples_frequency != 0
      for outfilei in outfile
        close(outfilei)
      end
    end

    output=output_result(mme,solMean,output_samples_frequency,meanu)
    return output
end
