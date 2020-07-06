function MT_MCMC_BayesianAlphabet(mme,df;causal_structure=false)
    ############################################################################
    chain_length             = mme.MCMCinfo.chain_length
    burnin                   = mme.MCMCinfo.burnin
    output_samples_frequency = mme.MCMCinfo.output_samples_frequency
    output_samples_file      = mme.MCMCinfo.output_samples_file
    estimate_variance        = mme.MCMCinfo.estimate_variance
    Rinv                     = mme.invweights
    update_priors_frequency  = mme.MCMCinfo.update_priors_frequency
    missing_phenotypes       = mme.MCMCinfo.missing_phenotypes
    constraint               = mme.MCMCinfo.constraint
    causal_structure         = causal_structure
    ############################################################################
    # Working Variables
    # 1) samples at current iteration (starting values default to zeros)
    # 2) posterior mean and variance at current iteration (zeros at the beginning)
    # 3) ycorr: phenotypes corrected for all effects
    ############################################################################
    #location parameters
    #mme.sol (starting values were set in runMCMC)
    mme.solMean, mme.solMean2  = zero(mme.sol),zero(mme.sol)
    #if methods == "BayesCC"  labels,BigPi,BigPiMean=setPi(Pi)  end
    ############################################################################
    # PRIORS
    ############################################################################
    ntraits        = mme.nModels
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
            wArray         = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,ntraits)#wArray is list reference of ycor
            ########################################################################
            #Arrays to save solutions for marker effects
            #In RR-BLUP and BayesL, only alphaArray is used.
            #In BayesA, B and C, alphaArray,deltaArray, betaArray are used.
            ########################################################################
            #starting values for marker effects(zeros) and location parameters (sol)
            #Mi.α  (starting values were set in get_genotypes)
            Mi.β          = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,ntraits) #BayesC
            Mi.δ          = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,ntraits) #BayesC
            Mi.meanDelta  = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,ntraits) #BayesC
            Mi.meanAlpha  = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,ntraits) #BayesC
            Mi.meanAlpha2 = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,ntraits) #BayesC

            for traiti = 1:Mi.ntraits
                Mi.β[traiti]          = copy(Mi.α[traiti])
                Mi.δ[traiti]          = ones(typeof(Mi.α[traiti][1]),Mi.nMarkers)
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
    ycorr       = vec(Matrix(mme.ySparse)-mme.X*mme.sol)
    if mme.M != 0
        for Mi in mme.M
            for traiti in 1:mme.nModels
                if Mi.α[traiti] != zero(Mi.α[traiti])
                    ycorr[(traiti-1)*Mi.nObs+1 : traiti*Mi.nObs] = ycorr[(traiti-1)*Mi.nObs+1 : traiti*Mi.nObs]
                                                                 - Mi.genotypes*Mi.α[traiti]
                end
            end
        end
    end
    wArray = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,ntraits)
    for traiti = 1:mme.nModels
        startPosi             = (traiti-1)*nObs  + 1
        ptr                   = pointer(ycorr,startPosi)
        wArray[traiti]        = unsafe_wrap(Array,ptr,nObs)
    end

    #Starting value for Ri is made based on missing value pattern
    #(imputed phenotypes will not used to compute first mmeRhs)
    Ri         = mkRi(mme,df,mme.invweights)
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
          outfile=output_MCMC_samples_setup(mme,chain_length-burnin,output_samples_frequency,output_samples_file)
    end

    ############################################################################
    #MCMC
    ############################################################################
    @showprogress "running MCMC..." for iter=1:chain_length
        ########################################################################
        # 1. Non-Marker Location Parameters
        ########################################################################
        if mme.MCMCinfo.missing_phenotypes==true
          ycorr[:]=sampleMissingResiduals(mme,ycorr)
        end
        # 1.1 Update Left-hand-side of MME
        mme.mmeLhs =  mme.X'Ri* mme.X #LHS for normal equation (no random effects)
        addVinv(mme)
        # 1.2 Update Right-hand-side of MME
        ycorr[:]   = ycorr[:] + mme.X*mme.sol
        mme.mmeRhs =  mme.X'Ri*ycorr

        Gibbs(mme.mmeLhs,mme.sol,mme.mmeRhs)
        ycorr[:] = ycorr[:] - mme.X*mme.sol
        ########################################################################
        # 2. Marker Effects
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
                 ########################################################################
                 # Marker Inclusion Probability
                 ########################################################################
                if Mi.estimatePi == true
                    samplePi(Mi.δ,Mi.π) #samplePi(deltaArray,Mi.π,labels)
                end
                ########################################################################
                # Marker Effects Variance
                ########################################################################
                if Mi.estimateVariance == true
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
                    end
                end
            end
        end
        ########################################################################
        # 3. Non-marker Variance Components
        ########################################################################
        if mme.MCMCinfo.estimate_variance == true
            ########################################################################
            # 3.1 Variance of Non-marker Random Effects
            # e.g, iid; polygenic effects (pedigree)
            ########################################################################
            sampleVCs(mme,mme.sol)
            ########################################################################
            # 3.2 Residual Variance
            ########################################################################
            sample_variance(mme,ycorr,Rinv,constraint=constraint)
            Ri = kron(inv(mme.R),spdiagm(0=>Rinv))
        end
        ########################################################################
        # 4. Causal Relationships among Phenotypes (Structure Equation Model)
        ########################################################################
        if causal_structure != false
            sample4λ = get_Λ(Y,mme.R,ycorr,Λy,mme.ySparse,causal_structure) #no missing phenotypes
        end
        ########################################################################
        # 2.5 Update priors using posteriors (empirical)
        ########################################################################
        if update_priors_frequency !=0 && iter%update_priors_frequency==0
            if mme.M!=0
                mme.M.scale = meanVara*(mme.df.marker-ntraits-1)
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
            output_MCMC_samples(mme,mme.sol,mme.R,(mme.pedTrmVec!=0 ? inv(mme.Gi) : false),outfile)
            if causal_structure != false
                writedlm(causal_structure_outfile,sample4λ',',')
            end

            nsamples = (iter-burnin)/output_samples_frequency
            mme.solMean   += (mme.sol - mme.solMean)/nsamples
            mme.solMean2  += (mme.sol .^2 - mme.solMean2)/nsamples
            meanVare  += (mme.R - meanVare)/nsamples
            meanVare2 += (mme.R .^2 - meanVare2)/nsamples

            if mme.pedTrmVec != 0
                G0Mean  += (inv(mme.Gi)  - G0Mean)/nsamples
                G0Mean2 += (inv(mme.Gi) .^2  - G0Mean2) /nsamples
            end
            if mme.M != 0
                for Mi in mme.M
                    for trait in 1:ntraits
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
        if iter%output_samples_frequency==0 && iter>burnin
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
    output=output_result(mme,output_samples_file,
                         mme.solMean,meanVare,
                         mme.pedTrmVec!=0 ? G0Mean : false,
                         mme.solMean2,meanVare2,
                         mme.pedTrmVec!=0 ? G0Mean2 : false)
    return output
end
