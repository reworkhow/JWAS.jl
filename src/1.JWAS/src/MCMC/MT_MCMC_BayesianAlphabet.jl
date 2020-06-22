function MT_MCMC_BayesianAlphabet(nIter,mme,df;
                        Rinv                       = false,
                        burnin                     = 0,
                        sol                        = false,
                        outFreq                    = 1000,
                        missing_phenotypes         = false,
                        constraint                 = false,
                        estimate_variance          = true,
                        output_samples_frequency   = 0,
                        update_priors_frequency    = 0,
                        output_file                = "MCMC_samples",
                        causal_structure           = false)

    ############################################################################
    # Pre-Check
    ############################################################################
    #starting values for location parameters(no marker) are sol
    sol,α              = sol[1:size(mme.mmeLhs,1)],sol[(size(mme.mmeLhs,1)+1):end]
    solMean, solMean2  = zero(sol),zero(sol)
    #if methods == "BayesCC"  labels,BigPi,BigPiMean=setPi(Pi)  end
    ############################################################################
    # PRIORS
    ############################################################################
    nTraits        = mme.nModels
    nObs           = length(mme.obsID)
    #save posterior mean for residual variance
    meanVare  = zero(mme.R)
    meanVare2 = zero(mme.R)
    #save poseterior mean for variances explained by polygenic effects (A) e.g Animal+ Maternal
    if mme.pedTrmVec != 0
        G0Mean,G0Mean2  = zero(mme.Gi),zero(mme.Gi)
    end
    #priors for marker effect variance
    if mme.M != 0
        for Mi in mme.M
            ########################################################################
            #Priors for marker covaraince matrix
            ########################################################################
            mGibbs                              = GibbsMats(Mi.genotypes,Rinv)
            Mi.mArray,Mi.mRinvArray,Mi.mpRinvm  = mGibbs.xArray,mGibbs.xRinvArray,mGibbs.xpRinvx
            if Mi.method=="BayesB"
                Mi.G          = fill(Mi.G,Mi.nMarkers)#locus specific covaraince matrix
                Mi.nLoci      = zeros(Mi.ntraits)
            end
            if Mi.method=="BayesL"# in BayesL mme.M.G is the scale Matrix, Sigma, in MTBayesLasso paper
                Mi.G         /= 4*(Mi.ntraits+1)
                Mi.scale     /= 4*(Mi.ntraits+1)
                gammaDist     = Gamma((Mi.ntraits+1)/2, 8) #8 (1/8): the scale (rate) parameter
                Mi.gammaArray = rand(gammaDist,Mi.nMarkers)
            end
            if Mi.method=="GBLUP"
                GBLUP_setup(Mi)
            end
            ########################################################################
            ##WORKING VECTORS
            wArray         = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,nTraits)#wArray is list reference of ycor
            ########################################################################
            #Arrays to save solutions for marker effects
            #In RR-BLUP and BayesL, only alphaArray is used.
            #In BayesA, B and C, alphaArray,deltaArray, betaArray are used.
            ########################################################################
            #starting values for marker effects(zeros) and location parameters (sol)
            Mi.β          = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,nTraits) #BayesC
            Mi.δ          = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,nTraits) #BayesC
            Mi.meanDelta  = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,nTraits) #BayesC
            Mi.α          = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,nTraits) #BayesC,BayesC0
            Mi.meanAlpha  = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,nTraits) #BayesC
            Mi.meanAlpha2 = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,nTraits) #BayesC
            for traiti = 1:Mi.ntraits
                Mi.α[traiti]          = copy(α[(traiti-1)*Mi.nMarkers+1:traiti*Mi.nMarkers])
                Mi.β[traiti]          = copy(α[(traiti-1)*Mi.nMarkers+1:traiti*Mi.nMarkers])
                Mi.δ[traiti]          = ones(typeof(α[1]),Mi.nMarkers)
                Mi.meanDelta[traiti]  = zero(Mi.δ[traiti])
                Mi.meanAlpha[traiti]  = zero(Mi.α[traiti])
                Mi.meanAlpha2[traiti] = zero(Mi.α[traiti])
            end
            Mi.π,Mi.mean_pi,Mi.mean_pi2 = copy(Mi.π),copy(Mi.π),copy(Mi.π)
            if Mi.estimatePi == true
              for key in keys(Mi.mean_pi)
                Mi.mean_pi[key]=0.0
                Mi.mean_pi2[key]=0.0
              end
            end
            Mi.meanVara       = zero(mme.R)  #variable to save variance for marker effect
            Mi.meanVara2      = zero(mme.R)  #variable to save variance for marker effect
            Mi.meanScaleVara  =  zero(mme.R) #variable to save Scale parameter for prior of marker effect variance
            Mi.meanScaleVara2 =  zero(mme.R) #variable to save Scale parameter for prior of marker effect variance
        end
    end
    ##Phenotypes CORRECTED for all effects
    ycorr       = vec(Matrix(mme.ySparse)-mme.X*sol)
    if mme.M != 0 && α!=zero(α)
        nObs     = mme.M.nObs
        nMarkers = mme.M.nMarkers
        for traiti in 1:mme.nModels
            ycorr[(traiti-1)*nObs+1 : traiti*nObs] = ycorr[(traiti-1)*nObs+1 : traiti*nObs] - mme.M.genotypes*α[(traiti-1)*nMarkers+1 : traiti*nMarkers]
        end
    end
    wArray = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,nTraits)
    for traiti = 1:mme.nModels
        startPosi             = (traiti-1)*nObs  + 1
        ptr                   = pointer(ycorr,startPosi)
        wArray[traiti]        = unsafe_wrap(Array,ptr,nObs)
    end

    #if starting values for marker effects are provided,
    #re-calculate mmeLhs and mmeRhs (non-genomic mixed model equation)
    if mme.M != 0 && α!=zero(α)
        ycorr[:]   = ycorr[:] + X*sol
        Ri         = mkRi(mme,df,mme.invweights)
        mme.mmeLhs = mme.X'Ri*mme.X
        mme.mmeRhs = mme.X'Ri*ycorr
        for random_term in mme.rndTrmVec
          random_term.GiOld = zero(random_term.GiOld)
        end
        addVinv(mme)
        for random_term in mme.rndTrmVec
          random_term.GiOld = copy(random_term.GiNew)
        end
    end
    ############################################################################
    # Starting values for SEM
    ############################################################################
    if causal_structure != false
        Y,Λy,causal_structure_outfile = SEM_setup(wArray,causal_structure,mme)
    end

    ############################################################################
    # SET UP OUTPUT MCMC samples
    ############################################################################
    if output_samples_frequency != 0
          outfile=output_MCMC_samples_setup(mme,nIter-burnin,output_samples_frequency,output_file)
    end

    ############################################################################
    #MCMC
    ############################################################################
    @showprogress "running MCMC..." for iter=1:nIter
        ########################################################################
        # 1.1. Non-Marker Location Parameters
        ########################################################################
        Gibbs(mme.mmeLhs,sol,mme.mmeRhs)
        ycorr[:] = ycorr[:] - mme.X*sol
        ########################################################################
        # 1.2 Marker Effects
        ########################################################################
        if mme.M != 0
            for Mi in mme.M
                if Mi.method in ["BayesC","BayesB","BayesA"]
                    locus_effect_variances = (Mi.method=="BayesC" ? fill(Mi.G,Mi.nMarkers) : Mi.G)
                    MTBayesABC!(Mi,wArray,mme.R,locus_effect_variances)
                elseif Mi.method == "RR-BLUP"
                    MTBayesC0!(Mi,wArray,mme.R)
                elseif Mi.method == "BayesL"
                    MTBayesL!(Mi,wArray,mme.R)
                elseif Mi.method == "GBLUP"
                    MTGBLUP!(Mi,wArray,ycorr,mme.R,Rinv)
                 end
                if Mi.estimatePi == true
                    samplePi(Mi.δ,Mi.π) #samplePi(deltaArray,Mi.π,labels)
                end
            end
        end
        ########################################################################
        # 2.1 Residual Covariance Matrix
        ########################################################################
        resVec = ycorr #here resVec is alias for ycor ***

        if missing_phenotypes==true
          resVec[:]=sampleMissingResiduals(mme,resVec)
        end

        if estimate_variance == true
            sample_variance(mme,resVec,Rinv,constraint=constraint)
        end
        Ri = kron(inv(mme.R),spdiagm(0=>Rinv))
        ########################################################################
        # -- LHS and RHS for conventional MME (No Markers)
        # -- Position: between new Ri and new Ai
        ########################################################################
        X          = mme.X
        mme.mmeLhs = X'Ri*X #LHS for normal equation (no random effects)
        ycorr[:]   = ycorr[:] + X*sol #same to ycorr[:]=resVec+X*sol
        mme.mmeRhs = X'Ri*ycorr
        ########################################################################
        # 2.1 Genetic Covariance Matrix (Polygenic Effects) (variance.jl)
        # 2.2 varainces for (iid) random effects;not required(empty)=>jump out
        ########################################################################
        if estimate_variance == true
            sampleVCs(mme,sol)
            addVinv(mme)
        end
        ########################################################################
        # 2.3 Marker Covariance Matrix
        ########################################################################

        if mme.M != 0 && mme.MCMCinfo.estimate_variance == true
            for Mi in mme.M
                SM    = zero(Mi.scale)
                if Mi.method == "BayesC"
                    for traiti = 1:Mi.ntraits
                        for traitj = traiti:Mi.ntraits
                            SM[traiti,traitj]   = (Mi.β[traiti]'Mi.β[traitj])
                            SM[traitj,traiti]   = SM[traiti,traitj]
                        end
                    end
                    Mi.G = rand(InverseWishart(Mi.df + Mi.nMarkers, convert(Array,Symmetric(Mi.scale + SM))))
                elseif Mi.method == "RR-BLUP"
                    for traiti = 1:Mi.ntraits
                        for traitj = traiti:Mi.ntraits
                            SM[traiti,traitj]   = (Mi.α[traiti]'Mi.α[traitj])
                            SM[traitj,traiti]   = SM[traiti,traitj]
                        end
                    end
                    Mi.G = rand(InverseWishart(Mi.df + Mi.nMarkers, convert(Array,Symmetric(Mi.scale + SM))))
                elseif Mi.method == "BayesL"
                    for traiti = 1:Mi.ntraits
                        alphai = Mi.α[traiti]./Mi.gammaArray
                        for traitj = traiti:Mi.ntraits
                            SM[traiti,traitj]   = (alphai'Mi.α[traitj])
                            SM[traitj,traiti]   = SM[traiti,traitj]
                        end
                    end
                    Mi.G = rand(InverseWishart(Mi.df + Mi.nMarkers, convert(Array,Symmetric(Mi.scale + SM))))
                    sampleGammaArray!(Mi.gammaArray,Mi.α,Mi.G)# MH sampler of gammaArray (Appendix C in paper)
                elseif Mi.method == "BayesB"
                    marker_effects_matrix = Mi.β[1]
                    for traiti = 2:Mi.ntraits
                        marker_effects_matrix = [marker_effects_matrix Mi.β[traiti]]
                    end
                    marker_effects_matrix = marker_effects_matrix'
                    beta2 = [marker_effects_matrix[:,i]*marker_effects_matrix[:,i]' for i=1:size(marker_effects_matrix,2)]
                    for markeri = 1:Mi.nMarkers
                        Mi.G[markeri] = rand(InverseWishart(Mi.df + 1, convert(Array,Symmetric(Mi.scale + beta2[markeri]))))
                    end
                elseif Mi.method == "GBLUP"
                    for traiti = 1:Mi.ntraits
                        for traitj = traiti:Mi.ntraits
                            alphaArrayi         = Mi.α[traiti]
                            alphaArrayj         = Mi.α[traitj]
                            SM[traiti,traitj]   = alphaArrayi'*Diagonal(1 ./Mi.D)*alphaArrayj
                            SM[traitj,traiti]   = SM[traiti,traitj]
                        end
                    end
                    Mi.G = rand(InverseWishart(Mi.df + Mi.nMarkers, convert(Array,Symmetric(Mi.scale + SM))))
                else
                    error("Sampling of marker effect variances is not available")
                end
            end
        end

        ########################################################################
        # 2.4 Causal Relationships among phenotypes (Structure Equation Model)
        ########################################################################
        if causal_structure != false
            sample4λ = get_Λ(Y,mme.R,ycorr,Λy,mme.ySparse,causal_structure) #no missing phenotypes
        end
        ########################################################################
        # 2.5 Update priors using posteriors (empirical)
        ########################################################################
        if update_priors_frequency !=0 && iter%update_priors_frequency==0
            if mme.M!=0
                mme.M.scale = meanVara*(mme.df.marker-nTraits-1)
            end
            if mme.pedTrmVec != 0
                mme.scalePed  = G0Mean*(mme.df.polygenic - size(mme.pedTrmVec,1) - 1)
            end
            mme.scaleRes  =  meanVare*(mme.df.residual-2)/mme.df.residual
            println("\n Update priors from posteriors.")
        end
        ########################################################################
        # 3.1 Save MCMC samples
        ########################################################################
        if iter>burnin && (iter-burnin)%output_samples_frequency == 0
            output_MCMC_samples(mme,sol,mme.R,(mme.pedTrmVec!=0 ? inv(mme.Gi) : false),outfile)
            if causal_structure != false
                writedlm(causal_structure_outfile,sample4λ',',')
            end

            nsamples = (iter-burnin)/output_samples_frequency
            solMean   += (sol - solMean)/nsamples
            solMean2  += (sol .^2 - solMean2)/nsamples
            meanVare  += (mme.R - meanVare)/nsamples
            meanVare2 += (mme.R .^2 - meanVare2)/nsamples

            if mme.pedTrmVec != 0
                G0Mean  += (inv(mme.Gi)  - G0Mean)/nsamples
                G0Mean2 += (inv(mme.Gi) .^2  - G0Mean2) /nsamples
            end
            if mme.M != 0
                for Mi in mme.M
                    for trait in 1:nTraits
                        Mi.meanAlpha[trait] += (Mi.α[trait] - Mi.meanAlpha[trait])/nsamples
                        Mi.meanAlpha2[trait]+= (Mi.α[trait].^2 - Mi.meanAlpha2[trait])/nsamples
                        Mi.meanDelta[trait] += (Mi.δ[trait] - Mi.meanDelta[trait])/nsamples
                    end
                    if Mi.estimatePi == true
                        for i in keys(Mi.π)
                          Mi.mean_pi[i] += (Mi.π[i]-Mi.mean_pi[i])/nsamples
                          Mi.mean_pi2[i] += (Mi.π[i].^2-Mi.mean_pi2[i])/nsamples
                        end
                    end
                    if Mi.method != "BayesB"
                        Mi.meanVara += (Mi.G - Mi.meanVara)/nsamples
                        Mi.meanVara2 += (Mi.G .^2 - Mi.meanVara2)/nsamples
                    end
                end
            end
        end
        ########################################################################
        # 3.2 Printout
        ########################################################################
        if iter%outFreq==0 && iter>burnin
            println("\nPosterior means at iteration: ",iter)
            println("Residual covariance matrix: \n",round.(meanVare,digits=6))
        end
    end

    ############################################################################
    # After MCMC
    ############################################################################
    if output_samples_frequency != 0
      for (key,value) in outfile
        close(value)
      end
      if causal_structure != false
        close(causal_structure_outfile)
      end
    end
    output=output_result(mme,output_file,
                         solMean,meanVare,
                         mme.pedTrmVec!=0 ? G0Mean : false,
                         solMean2,meanVare2,
                         mme.pedTrmVec!=0 ? G0Mean2 : false)
    return output
end
