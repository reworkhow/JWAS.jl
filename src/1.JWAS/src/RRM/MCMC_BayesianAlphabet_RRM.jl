################################################################################
#MCMC for RR-BLUP, BayesC, BayesCπ and "conventional (no markers).
#
################################################################################
function MCMC_BayesianAlphabet_RRM(mme,df;
                     Φ                          = false,
                     burnin                     = 0,
                     outFreq                    = 1000,
                     output_samples_frequency   = 0,
                     update_priors_frequency    = 0,
                     nIter                      = 0,
                     output_folder              = false)
    ############################################################################
    # Working Variables
    # 1) samples at current iteration (starting values at the beginning, defaulting to zeros)
    # 2) posterior mean and variance at current iteration (zeros at the beginning)
    # 3) ycorr, phenotypes corrected for all effects
    ############################################################################
    mme.lhsVec = Symbol.(string.(1:size(mme.MCMCinfo.RRM,2)))
    #transition matrix from yfull to yobs
    T,whichzeros = matrix_yfull_to_yobs(df[!,1],mme.ySparse,df[!,:time])
    Z  = mkmat_incidence_factor(unique(df[!,1]) ,mme.M[1].obsID)
    if mme.M != 0
        for Mi in mme.M
            Mi.genotypes      = Z*Mi.genotypes
            Mi.obsID          = unique(df[!,1])
            Mi.nObs           = length(Mi.obsID)
            Mi.meanVara       = zero(Mi.G) #variable to save variance for marker effect
            Mi.meanVara2      = zero(Mi.G)
            Mi.meanScaleVara  = zero(Mi.G) #variable to save Scale parameter for prior of marker effect variance
            Mi.meanScaleVara2 = zero(Mi.G)
        end
    end

    #size of Φ
    ntime, ncoeff      = size(Φ)
    mme.solMean, mme.solMean2  = zero(mme.sol),zero(mme.sol)
    #residual variance
    mme.meanVare = mme.meanVare2 = 0.0


    #polygenic effects (A), e.g, Animal+ Maternal
    if mme.pedTrmVec != 0
       mme.G0Mean,mme.G0Mean2  = zero(mme.Gi),zero(mme.Gi)
    end

    if mme.M != 0
        for Mi in mme.M
            #marker effects
            mGibbs                        = GibbsMats(Mi.genotypes, ones(size(Mi.genotypes,1)))
            Mi.nObs,Mi.nMarkers, M              = mGibbs.nrows,mGibbs.ncols,mGibbs.X
            Mi.mArray,Mi.mRinvArray,Mi.mpRinvm = mGibbs.xArray,mGibbs.xRinvArray,mGibbs.xpRinvx
            Mi.β          = [copy(Mi.α[1]) for coeffi = 1:ncoeff] #partial marker effeccts used in BayesB
            Mi.δ          = [ones(typeof(Mi.α[1][1]),Mi.nMarkers) for coeffi = 1:ncoeff] #inclusion indicator for marker effects
            Mi.meanDelta  = [zero(Mi.β[coeffi]) for coeffi = 1:ncoeff]
            Mi.α          = [copy(Mi.α[1]) for coeffi = 1:ncoeff]
            Mi.meanAlpha  = [zero(Mi.α[coeffi]) for coeffi = 1:ncoeff]
            Mi.meanAlpha2 = [zero(Mi.α[coeffi]) for coeffi = 1:ncoeff] #marker effects
        end
    end

    #phenotypes corrected for all effects #assumming starting values zeros
    ycorr    = vec(Matrix(mme.ySparse))
    yfull    = T'ycorr
    wArray   = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,ntime)#wArray is list reference of ycor
    if mme.M != 0
        for Mi in mme.M
            for timei = 1:ntime
                startPosi              = (timei-1)*Mi.nObs  + 1
                ptr                    = pointer(yfull,startPosi)
                wArray[timei]         = unsafe_wrap(Array,ptr,Mi.nObs)
            end
        end
    end

    #Pi value
    if mme.M != 0
        for Mi in mme.M
            Mi.π,Mi.mean_pi,Mi.mean_pi2 = copy(Mi.π),copy(Mi.π),copy(Mi.π)
            if Mi.estimatePi == true
                for key in keys(Mi.mean_pi)
                  Mi.mean_pi[key]  =0.0
                  Mi.mean_pi2[key] =0.0
                end
            end
        end
    end

    if mme.M != 0
        for Mi in mme.M
            Mi.mΦΦArray=get_mΦΦarray(Φ,T,Mi.mArray)
        end
    end

    ############################################################################
    #  SET UP OUTPUT MCMC samples
    ############################################################################
    if output_samples_frequency != 0
        outfile=output_MCMC_samples_setup(mme,nIter-burnin,output_samples_frequency,output_folder*"/MCMC_samples")
    end
    ############################################################################
    # MCMC (starting values for sol (zeros);  mme.RNew; G0 are used)
    ############################################################################
    @showprogress "running MCMC for" for iter=1:nIter
        ########################################################################
        # 1.1. Non-Marker Location Parameters
        ########################################################################
        ycorr[:] = ycorr + mme.X*mme.sol
        mme.mmeRhs = mme.X'ycorr
        Gibbs(mme.mmeLhs,mme.sol,mme.mmeRhs,mme.R)
        ycorr[:] = ycorr - mme.X*mme.sol

        ########################################################################
        # 1.2 Marker Effects
        ########################################################################
        yfull[:]   = T'ycorr

        if mme.M != 0
            for Mi in mme.M
                if Mi.method in ["BayesC","BayesB","BayesA"]
                    locus_effect_variances = (Mi.method == "BayesC" ? fill(Mi.G,Mi.nMarkers) : Mi.G)
                    BayesABCRRM!(Mi.mArray,Mi.mpRinvm,wArray,yfull,
                        Mi.β,
                        Mi.δ,
                        Mi.α,
                        mme.R,locus_effect_variances,Mi.π,
                        Φ, whichzeros, Mi.mΦΦArray)
                    end
                end
            end


        # MTBayesABC!(mArray,mpm,wArray,
        #           betaArray,
        #           deltaArray,
        #           alphaArray,
        #           [mme.RNew 0;0 mme.RNew],locus_effect_variances,BigPi)
        ycorr[:]     = T*yfull

        #sample Pi
        if mme.M != 0
            for Mi in mme.M
                if Mi.estimatePi == true
                    samplePi(Mi.δ,Mi.π)
                end
            end
        end

        ########################################################################
        # 2.1 Genetic Covariance Matrix (Polygenic Effects) (variance.jl)
        # 2.2 varainces for (iid) random effects;not required(empty)=>jump out
        ########################################################################
        sampleVCs(mme,mme.sol)
        addVinv(mme)
        ########################################################################
        # 2.3 Residual Variance
        ########################################################################
        mme.ROld = mme.R
        mme.R = sample_variance(ycorr, length(ycorr), mme.df.residual, false, false)

        ########################################################################
        # 2.4 Marker Effects Variance
        ########################################################################
        if mme.M != 0
           for Mi in mme.M
               if Mi.method in ["RR-BLUP","BayesC","GBLUP"]
                   data = (Mi.method == "BayesC" ? Mi.β : Mi.α)
                   Mi.G =sample_variance(data, Mi.nMarkers, Mi.df, Mi.scale, false, false)
               end
           end
       end
        ########################################################################
        # 2.6 sample Scale parameter in prior for marker effect variances
        ########################################################################
        if mme.M != 0
            for Mi in mme.M
                if Mi.estimateScale == true
                    a = size(Mi.G,1)*Mi.df/2  + 1
                    b = sum(Mi.df ./ (2*Mi.G )) + 1
                    Mi.scale = rand(Gamma(a,1/b))
                end
            end
        end
        ########################################################################
        # 3.1 Save MCMC samples
        ########################################################################
        if iter>burnin && (iter-burnin)%output_samples_frequency == 0
            output_MCMC_samples(mme,mme.R,(mme.pedTrmVec!=0 ? inv(mme.Gi) : false),outfile)
            nsamples = (iter-burnin)/output_samples_frequency
            mme.solMean   += (mme.sol - mme.solMean)/nsamples
            mme.solMean2  += (mme.sol .^2 - mme.solMean2)/nsamples
            mme.meanVare  += (mme.R - mme.meanVare)/nsamples
            mme.meanVare2 += (mme.R .^2 - mme.meanVare2)/nsamples

            if mme.pedTrmVec != 0
                mme.G0Mean  += (inv(mme.Gi)  - mme.G0Mean )/nsamples
            end
            if mme.M != 0
                for Mi in mme.M
                    for i in 1:ncoeff
                        Mi.meanAlpha[i] += (Mi.α[i] - Mi.meanAlpha[i])/nsamples
                        Mi.meanAlpha2[i]+= (Mi.α[i].^2 - Mi.meanAlpha2[i])/nsamples
                        Mi.meanDelta[i] += (Mi.δ[i] - Mi.meanDelta[i])/nsamples
                    end
                    if Mi.estimatePi == true
                        for i in keys(Mi.π)
                              Mi.mean_pi[i] += (Mi.π[i]-Mi.mean_pi[i])/nsamples
                              Mi.mean_pi2[i] += (Mi.π[i].^2-Mi.mean_pi2[i])/nsamples
                        end
                    end
                    if Mi.method != "BayesB"
                        Mi.meanVara += (Mi.G - Mi.meanVara)/nsamples
                        Mi.meanVara2 += (Mi.G.^2 - Mi.meanVara2)/nsamples
                    end

                    if Mi.estimateScale == true
                        Mi.meanScaleVara += (Mi.scale - Mi.meanScaleVara)/nsamples
                        Mi.meanScaleVara2 += (Mi.scale .^2 - Mi.meanScaleVara2)/nsamples
                    end
                end
            end
        end
        ########################################################################
        # 3.2 Printout
        ########################################################################
        if iter%outFreq==0 && iter>burnin
            println("\nPosterior means at iteration: ",iter)
            println("Residual variance: ",round(mme.meanVare,digits=6))
        end
    end

    ############################################################################
    # After MCMC
    ############################################################################
    if output_samples_frequency != 0
      for (key,value) in outfile
        close(value)
      end
    end
    output=output_result(mme,output_folder,
                          mme.solMean,mme.meanVare,
                          mme.pedTrmVec!=0 ? mme.G0Mean : false,
                          mme.solMean2,mme.meanVare2,
                          mme.pedTrmVec!=0 ? mme.G0Mean2 : false)

    return output
end
