################################################################################
#MCMC for RR-BLUP, BayesC, BayesCπ and "conventional (no markers).
#
################################################################################
function MCMC_BayesianAlphabet_RRM(mme,df;
                     Φ                          = false,
                     burnin                     = 0,
                     Pi                         = 0.0,
                     starting_value             = false,
                     outFreq                    = 1000,
                     methods                    = "BayesC",
                     output_samples_frequency   = 0,
                     nIter                      = 100,
                     output_folder                = "MCMC_samples")
    ############################################################################
    # Working Variables
    # 1) samples at current iteration (starting values at the beginning, defaulting to zeros)
    # 2) posterior mean and variance at current iteration (zeros at the beginning)
    # 3) ycorr, phenotypes corrected for all effects
    ############################################################################
    mme.lhsVec = Symbol.(string.(1:size(mme.MCMCinfo.RRM,2)))
    mme.solMean, mme.solMean2  = zero(mme.sol),zero(mme.sol)
    if mme.M != 0
        for Mi in mme.M
            Mi.meanVara  = zero(Mi.G) #variable to save variance for marker effect
            Mi.meanVara2 = zero(Mi.G) #variable to save variance for marker effect
        end
    end
    #transition matrix from yfull to yobs
    T,whichzeros = matrix_yfull_to_yobs(df[!,1],mme.ySparse,df[!,:time])
    Z  = mkmat_incidence_factor(unique(df[!,1]) ,mme.obsID)
    mme.obsID     = unique(df[!,1])

    if mme.M != 0
        for Mi in mme.M
            Mi.genotypes = Z * Mi.genotypes
            Mi.nObs      = length(mme.obsID)
        end
    end
    #size of Φ
    ntime, ncoeff      = size(Φ)
    #location parameters
    sol,α       = starting_value[1:size(mme.mmeLhs,1)],starting_value[(size(mme.mmeLhs,1)+1):end]
    solMean, solMean2  = zero(sol),zero(sol)


    meanVare = meanVare2 = 0.0
    if mme.M != 0
        for Mi in mme.M
            #residual variance
            meanVara  = zero(Mi.G) #variable to save variance for marker effect
            meanVara2 = zero(Mi.G) #variable to save variance for marker effect
            meanScaleVara  =  zero(Mi.G) #variable to save Scale parameter for prior of marker effect variance
            meanScaleVara2 =  zero(Mi.G) #variable to save Scale parameter for prior of marker effect variance

            #marker effects
            mGibbs                        = GibbsMats(Mi.genotypes, ones(size(Mi.genotypes,1)))
            Mi.nObs,Mi.nMarkers, M              = mGibbs.nrows,mGibbs.ncols,mGibbs.X
            Mi.mArray,Mi.mRinvArray,Mi.mpRinvm  = mGibbs.xArray,mGibbs.xRinvArray,mGibbs.xpRinvx

            Mi.β  = [copy(Mi.α[1]) for coeffi = 1:ncoeff] #partial marker effeccts used in BayesB
            Mi.δ   = [ones(typeof(Mi.α[1][1]),Mi.nMarkers) for coeffi = 1:ncoeff] #inclusion indicator for marker effects
            Mi.meanDelta  = [zero(Mi.β[coeffi]) for coeffi = 1:ncoeff]
            Mi.α  = [copy(Mi.α[1]) for coeffi = 1:ncoeff]
            Mi.meanAlpha  = [zero(Mi.α[coeffi]) for coeffi = 1:ncoeff]
            Mi.meanAlpha2 = [zero(Mi.α[coeffi]) for coeffi = 1:ncoeff] #marker effects
        end
    end
    #polygenic effects (A), e.g, Animal+ Maternal
    if mme.pedTrmVec != 0
        G0Mean,G0Mean2  = zero(mme.Gi),zero(mme.Gi)
    end

    #phenotypes corrected for all effects #assumming starting values zeros
    ycorr    = vec(Matrix(mme.ySparse)-mme.X*mme.sol)
    yfull    = T'ycorr
    wArray         = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,ntime)#wArray is list reference of ycor
    if mme.M != 0
        for Mi in mme.M
            for timei = 1:ntime
                startPosi              = (timei-1)*Mi.nObs  + 1
                ptr                    = pointer(yfull,startPosi)
                wArray[timei]         = unsafe_wrap(Array,ptr,Mi.nObs)
            end

            #Pi value
            Mi.π,Mi.mean_pi,Mi.mean_pi2 = copy(Mi.π),copy(Mi.π),copy(Mi.π)
            if Mi.estimatePi == true
                for key in keys(Mi.mean_pi)
                        Mi.mean_pi[key]=0.0
                        Mi.mean_pi2[key]=0.0
                end
            end
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
    @showprogress "running MCMC ..." for iter=1:nIter
        ########################################################################
        # 1.1. Non-Marker Location Parameters
        ########################################################################
        ycorr[:] = ycorr + mme.X*sol
        mme.mmeRhs =  mme.X'ycorr
        Gibbs(mme.mmeLhs,mme.sol, mme.mmeRhs, mme.R)
        ycorr[:] = ycorr - mme.X*sol

        ########################################################################
        # 1.2 Marker Effects
        ########################################################################
        yfull[:]   = T'ycorr

        if mme.M != 0
            for Mi in mme.M
                if Mi.method in ["BayesC","BayesB","BayesA"]
                    locus_effect_variances = (methods=="BayesC" ? fill(Mi.G,Mi.nMarkers) : Mi.G)
                    mΦΦArray=get_mΦΦarray(Φ,T,Mi.mArray)
                    BayesABCRRM!(Mi,wArray,yfull,
                    mme.R,locus_effect_variances,
                    Φ, whichzeros, mΦΦArray)
                end
            end
        end

        # MTBayesABC!(mArray,mpm,wArray,
        #           betaArray,
        #           deltaArray,
        #           alphaArray,
        #           [mme.RNew 0;0 mme.RNew],locus_effect_variances,BigPi)
        ycorr     = T*yfull

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
        mme.R = sample_variance(ycorr, length(ycorr), mme.df.residual, mme.scaleR, mme.invweights)

        ########################################################################
        # 2.4 Marker Effects Variance
        ########################################################################
        if mme.M != 0
            for Mi in mme.M
                if Mi.method in ["RR-BLUP","BayesC","GBLUP"]
                    data = (Mi.method == "BayesC" ? Mi.β : Mi.α)
                    Mi.G =sample_variance(data, Mi.nMarkers, Mi.df, Mi.scale,false, false)
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
            nsamples = (iter-burnin)/output_samples_frequency
            mme.solMean   += (mme.sol - mme.solMean)/nsamples
            mme.solMean2  += (mme.sol .^2 - mme.solMean2)/nsamples
            mme.meanVare  += (mme.R - mme.meanVare)/nsamples
            mme.meanVare2 += (mme.R .^2 - mme.meanVare2)/nsamples

            if mme.pedTrmVec != 0
                mme.G0Mean  += (inv(mme.Gi)  - mme.G0Mean )/nsamples
                mme.G0Mean2 += (inv(mme.Gi) .^2  - mme.G0Mean2 )/nsamples
            end

            if mme.M != 0
                for Mi in mme.M
                    for coeffi in 1:ncoeff
                        Mi.meanAlpha[coeffi] += (Mi.α[coeffi] - Mi.meanAlpha[coeffi])/nsamples
                        Mi.meanAlpha2[coeffi]+= (Mi.α[coeffi].^2 - Mi.meanAlpha2[coeffi])/nsamples
                        Mi.meanDelta[coeffi] += (Mi.δ[coeffi] - Mi.meanDelta[coeffi])/nsamples
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
                    if Mi.estimateScale == true
                        Mi.meanScaleVara += (Mi.scale - Mi.meanScaleVara)/nsamples
                        Mi.meanScaleVara2 += (Mi.scale .^2 - Mi.meanScaleVara2)/nsamples
                    end
                end
            end

            #mean and variance of posterior distribution
            output_MCMC_samples(mme,mme.R,(mme.pedTrmVec!=0 ? inv(mme.Gi) : false),outfile)
        end

        ########################################################################
        # 3.2 Printout
        ########################################################################
        if iter%mme.MCMCinfo.printout_frequency==0 && iter>burnin
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
                          mme.pedTrmVec!=0 ? mme.G0Mean2 : false
                          #=
                          mme.M != 0 ? meanalphaArray : false,
                          mme.M != 0 ? meandeltaArray : false,
                          mme.M != 0 ? meanVara : false,
                          mme.M != 0 ? estimatePi : false,
                          mme.M != 0 ? BigPiMean : false,
                          mme.M != 0 ? estimateScale : false,
                          mme.M != 0 ? meanScaleVara : false,
                          solMean2,meanVare2,
                          mme.pedTrmVec!=0 ? G0Mean2 : false,
                          mme.M != 0 ? meanalphaArray2 : false,
                          mme.M != 0 ? meanVara2 : false,
                          mme.M != 0 ? BigPiMean2 : false,
                          mme.M != 0 ? meanScaleVara2 : false
                          =#)

    return output
end
