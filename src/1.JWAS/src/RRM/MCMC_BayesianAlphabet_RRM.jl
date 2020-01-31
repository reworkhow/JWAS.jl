################################################################################
#MCMC for RR-BLUP, BayesC, BayesCπ and "conventional (no markers).
#
################################################################################
function MCMC_BayesianAlphabet_RRM(nIter,mme,df;
                     Φ                          = false,
                     burnin                     = 0,
                     Pi                         = 0.0,
                     estimatePi                 = false,
                     estimateScale              = false,
                     starting_value             = false,
                     outFreq                    = 1000,
                     methods                    = "BayesC",
                     output_samples_frequency   = 0,
                     update_priors_frequency    = 0,
                     output_file                = "MCMC_samples")
    ############################################################################
    # Working Variables
    # 1) samples at current iteration (starting values at the beginning, defaulting to zeros)
    # 2) posterior mean and variance at current iteration (zeros at the beginning)
    # 3) ycorr, phenotypes corrected for all effects
    ############################################################################
    mme.lhsVec = Symbol.(string.(1:size(mme.MCMCinfo.RRM,2)))
    #transition matrix from yfull to yobs
    T,whichzeros = matrix_yfull_to_yobs(df[!,1],mme.ySparse,df[!,:time])
    Z  = mkmat_incidence_factor(unique(df[!,1]) ,mme.M.obsID)
    mme.M.genotypes = Z*mme.M.genotypes
    mme.M.obsID     = unique(df[!,1])
    mme.M.nObs      = length(mme.M.obsID)
    #size of Φ
    ntime, ncoeff      = size(Φ)
    #location parameters
    sol,α       = starting_value[1:size(mme.mmeLhs,1)],starting_value[(size(mme.mmeLhs,1)+1):end]
    solMean, solMean2  = zero(sol),zero(sol)
    #residual variance
    meanVare = meanVare2 = 0.0
    meanVara  = zero(mme.M.G) #variable to save variance for marker effect
    meanVara2 = zero(mme.M.G) #variable to save variance for marker effect
    meanScaleVara  =  zero(mme.M.G) #variable to save Scale parameter for prior of marker effect variance
    meanScaleVara2 =  zero(mme.M.G) #variable to save Scale parameter for prior of marker effect variance
    #polygenic effects (A), e.g, Animal+ Maternal
    if mme.pedTrmVec != 0
       G0Mean,G0Mean2  = zero(mme.Gi),zero(mme.Gi)
    end
    #marker effects
    mGibbs                        = GibbsMats(mme.M.genotypes)
    nObs,nMarkers, M              = mGibbs.nrows,mGibbs.ncols,mGibbs.X
    mArray,mpm,mRinvArray,mpRinvm = mGibbs.xArray,mGibbs.xpx,mGibbs.xRinvArray,mGibbs.xpRinvx

    #starting values for marker effects(zeros) and location parameters (sol)
    betaArray       = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,ncoeff) #BayesC
    deltaArray      = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,ncoeff) #BayesC
    meandeltaArray  = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,ncoeff) #BayesC
    alphaArray      = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,ncoeff) #BayesC,BayesC0
    meanalphaArray  = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,ncoeff) #BayesC
    meanalphaArray2 = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,ncoeff) #BayesC
    for coeffi = 1:ncoeff
        betaArray[coeffi]      = copy(α[(coeffi-1)*nMarkers+1:coeffi*nMarkers])
        deltaArray[coeffi]     = ones(typeof(α[1]),nMarkers)
        meandeltaArray[coeffi] = zero(betaArray[coeffi])
        alphaArray[coeffi]     = copy(α[(coeffi-1)*nMarkers+1:coeffi*nMarkers])
        meanalphaArray[coeffi] = zero(alphaArray[coeffi])
        meanalphaArray2[coeffi]= zero(alphaArray[coeffi])
    end

    #phenotypes corrected for all effects #assumming starting values zeros
    ycorr    = vec(Matrix(mme.ySparse))
    yfull    = T'ycorr
    wArray         = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,ntime)#wArray is list reference of ycor
    for timei = 1:ntime
        startPosi              = (timei-1)*nObs  + 1
        #ptr                    = pointer(yfull,startPosi)
        ptr                    = pointer(ycorr,startPosi)
        wArray[timei]         = unsafe_wrap(Array,ptr,nObs)
    end

    #Pi value
    BigPi,BigPiMean,BigPiMean2 = copy(Pi),copy(Pi),copy(Pi)
    for key in keys(BigPiMean)
      BigPiMean[key]=0.0
      BigPiMean2[key]=0.0
    end

    mΦΦArray=get_mΦΦarray(Φ,T,mArray)
    ############################################################################
    #  SET UP OUTPUT MCMC samples
    ############################################################################
    if output_samples_frequency != 0
        outfile=output_MCMC_samples_setup(mme,nIter-burnin,output_samples_frequency,output_file)
    end
    ############################################################################
    # MCMC (starting values for sol (zeros);  mme.RNew; G0 are used)
    ############################################################################
    @showprogress "running MCMC for "*methods*"..." for iter=1:nIter
        ########################################################################
        # 1.1. Non-Marker Location Parameters
        ########################################################################
        ycorr = ycorr + mme.X*sol
        rhs = mme.X'ycorr
        Gibbs(mme.mmeLhs,sol,rhs,mme.RNew)
        ycorr = ycorr - mme.X*sol

        ########################################################################
        # 1.2 Marker Effects
        ########################################################################
        #yfull[:]   = T'ycorr
        locus_effect_variances = (methods=="BayesC" ? fill(mme.M.G,nMarkers) : mme.M.G)
        # BayesABCRRM!(mArray,mpm,wArray,yfull,
        #             betaArray,
        #             deltaArray,
        #             alphaArray,
        #             mme.RNew,locus_effect_variances,BigPi,
        #             Φ, whichzeros, mΦΦArray)
        MTBayesABC!(mArray,mpm,wArray,
                  betaArray,
                  deltaArray,
                  alphaArray,
                  [mme.RNew 0;0 mme.RNew],locus_effect_variances,BigPi)
        #ycorr     = T*yfull

        #sample Pi
        if estimatePi == true
            samplePi(deltaArray,BigPi,BigPiMean,iter)
        end
        ########################################################################
        # 2.1 Genetic Covariance Matrix (Polygenic Effects) (variance.jl)
        # 2.2 varainces for (iid) random effects;not required(empty)=>jump out
        ########################################################################
        sampleVCs(mme,sol)
        addVinv(mme)
        ########################################################################
        # 2.3 Residual Variance
        ########################################################################
        mme.ROld = mme.RNew
        mme.RNew = sample_variance(ycorr, length(ycorr), mme.df.residual, mme.scaleRes)

        ########################################################################
        # 2.4 Marker Effects Variance
        ########################################################################
        SM    = zero(mme.M.scale)
        if methods == "BayesC"
            for coeffi = 1:ncoeff
                for coeffj = coeffi:ncoeff
                    SM[coeffi,coeffj]   = (betaArray[coeffi]'betaArray[coeffj])
                    SM[coeffj,coeffi]   = SM[coeffi,coeffj]
                end
            end
            mme.M.G = rand(InverseWishart(mme.df.marker + nMarkers, convert(Array,Symmetric(mme.M.scale + SM))))
        end
        ########################################################################
        # 2.6 sample Scale parameter in prior for marker effect variances
        ########################################################################
        if estimateScale == true
            a = size(mme.M.G,1)*mme.df.marker/2   + 1
            b = sum(mme.df.marker ./ (2*mme.M.G )) + 1
            mme.M.scale = rand(Gamma(a,1/b))
        end
        ########################################################################
        # 3.1 Save MCMC samples
        ########################################################################
        if iter>burnin && (iter-burnin)%output_samples_frequency == 0
            output_MCMC_samples(mme,sol,mme.RNew,(mme.pedTrmVec!=0 ? inv(mme.Gi) : false),BigPi,alphaArray,vec(mme.M.G),outfile)
            nsamples = (iter-burnin)/output_samples_frequency
            solMean   += (sol - solMean)/nsamples
            solMean2  += (sol .^2 - solMean2)/nsamples
            meanVare  += (mme.RNew - meanVare)/nsamples
            meanVare2 += (mme.RNew .^2 - meanVare2)/nsamples

            if mme.pedTrmVec != 0
                G0Mean  += (inv(mme.Gi)  - G0Mean )/nsamples
                G0Mean2 += (inv(mme.Gi) .^2  - G0Mean2 )/nsamples
            end
            for i in 1:ncoeff
                meanalphaArray[i] += (alphaArray[i] - meanalphaArray[i])/nsamples
                meanalphaArray2[i]+= (alphaArray[i].^2 - meanalphaArray2[i])/nsamples
                meandeltaArray[i] += (deltaArray[i] - meandeltaArray[i])/nsamples
            end
            if estimatePi == true
                for i in keys(BigPi)
                  BigPiMean[i] += (BigPi[i]-BigPiMean[i])/nsamples
                  BigPiMean2[i] += (BigPi[i].^2-BigPiMean2[i])/nsamples
                end
            end
            if methods != "BayesB"
                meanVara += (mme.M.G - meanVara)/nsamples
                meanVara2 += (mme.M.G .^2 - meanVara2)/nsamples
            end
            if estimateScale == true
                meanScaleVara += (mme.M.scale - meanScaleVara)/nsamples
                meanScaleVara2 += (mme.M.scale .^2 - meanScaleVara2)/nsamples
            end
        end
        ########################################################################
        # 3.2 Printout
        ########################################################################
        if iter%outFreq==0 && iter>burnin
            println("\nPosterior means at iteration: ",iter)
            println("Residual variance: ",round(meanVare,digits=6))
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
     output=output_result(mme,output_file,
                          solMean,meanVare,
                          mme.pedTrmVec!=0 ? G0Mean : false,
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
                          mme.M != 0 ? meanScaleVara2 : false)

    return output
end
