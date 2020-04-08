################################################################################
#MCMC for RR-BLUP, GBLUP, BayesABC, and conventional (no markers) methods
#
#GBLUP: pseudo markers are used.
#y = μ + a + e with mean(a)=0,var(a)=Gσ²=MM'σ² and G = LDL' <==>
#y = μ + Lα +e with mean(α)=0,var(α)=D*σ² : L orthogonal
#y2hat = cov(a2,a)*inv(var(a))*L*αhat =>
#      = (M2*M')*inv(MM')*L*αhat = (M2*M')*(L(1./D)L')(Lαhat)=(M2*M'*L(1./D))αhat
#
#Threshold trait:
#1)Sorensen and Gianola, Likelihood, Bayesian, and MCMC Methods in Quantitative
#Genetics
#2)Wang et al.(2013). Bayesian methods for estimating GEBVs of threshold traits.
#Heredity, 110(3), 213–219.
################################################################################
function MCMC_BayesianAlphabet(nIter,mme,df;
                     Rinv                       = false,
                     burnin                     = 0,
                     π                          = 0.0,
                     estimatePi                 = false,
                     estimate_variance          = false,
                     estimateScale              = false,
                     starting_value             = false,
                     outFreq                    = 1000,
                     methods                    = "BayesC",
                     output_samples_frequency   = 0,
                     update_priors_frequency    = 0,
                     output_file                = "MCMC_samples",
                     categorical_trait          = false)

    ############################################################################
    # Categorical Traits (starting values for maker effects defaulting to 0s)
    ############################################################################
    if categorical_trait == true
        category_obs,threshold = categorical_trait_setup(mme,starting_value)
    end
    ############################################################################
    # Working Variables
    # 1) samples at current iteration (starting values at the beginning, defaulting to zeros)
    # 2) posterior mean and variance at current iteration (zeros at the beginning)
    # 3) ycorr, phenotypes corrected for all effects
    ############################################################################
    #location parameters
    sol                = starting_value[1:size(mme.mmeLhs,1)]
    solMean, solMean2  = zero(sol),zero(sol)
    #residual variance
    if categorical_trait == false
        meanVare = meanVare2 = 0.0
    else
        mme.RNew = meanVare = meanVare2 = 1.0
    end
    #polygenic effects (A), e.g, Animal+ Maternal
    if mme.pedTrmVec != 0
       G0Mean,G0Mean2  = zero(mme.Gi),zero(mme.Gi)
    end
    #marker effects
    if mme.M != 0
        mGibbs                     = GibbsMats(mme.M.genotypes,Rinv)
        mme.M.mArray,mme.M.mRinvArray,mme.M.mpRinvm  = mGibbs.xArray,mGibbs.xRinvArray,mGibbs.xpRinvx

        if methods=="BayesB" #α=β.*δ
            mme.M.G        = fill(mme.M.G,mme.M.nMarkers) #a scalar in BayesC but a vector in BayeB
        end
        if methods=="BayesL"         #in the MTBayesLasso paper
            mme.M.G   /= 8           #mme.M.G is the scale Matrix, Sigma
            mme.M.scale /= 8
            gammaDist  = Gamma(1, 8) #8 is the scale parameter of the Gamma distribution (1/8 is the rate parameter)
            mme.M.gammaArray = rand(gammaDist,mme.M.nMarkers)
        end
        if methods=="GBLUP"
            mme.M.genotypes  = mme.M.genotypes ./ sqrt.(2*mme.M.alleleFreq.*(1 .- mme.M.alleleFreq))
            G       = (mme.M.genotypes*mme.M.genotypes'+ I*0.00001)/mme.M.nMarkers
            eigenG  = eigen(G)
            L       = eigenG.vectors
            D       = eigenG.values
            # α is pseudo marker effects of length nobs (starting values = L'(starting value for BV)
            mme.M.nMarkers= mme.M.nObs
            #reset parameters in output
            M2   = mme.output_genotypes ./ sqrt.(2*mme.M.alleleFreq.*(1 .- mme.M.alleleFreq))
            M2Mt = M2*mme.M.genotypes'/mme.M.nMarkers
            mme.output_genotypes = M2Mt*L*Diagonal(1 ./D)
            #reset parameter in mme.M
            mme.M.G         = mme.M.genetic_variance
            mme.M.scale     = mme.M.G*(mme.df.marker-2)/mme.df.marker
            mme.M.markerID  = string.(1:mme.M.nObs) #pseudo markers of length=nObs
            mme.M.genotypes = L
        end
        mme.M.α                            = starting_value[(size(mme.mmeLhs,1)+1):end]
        if methods == "GBLUP"
            mme.M.α  = L'mme.M.α
        end
        mme.M.β                                  = copy(α)          #partial marker effeccts used in BayesB
        mme.M.δ                                  = ones(typeof(α[1]),mme.M.nMarkers)   #inclusion indicator for marker effects
        mme.M.meanAlpha,mme.M.meanAlpha2         = zero(α),zero(α)  #marker effects
        mme.M.meanDelta                          = zero(δ)          #inclusion indicator
        mme.M.mean_pi,mme.M.mean_pi2             = 0.0,0.0          #inclusion probability
        mme.M.meanVara,mme.M.meanVara2           = 0.0,0.0          #marker effect variances
        mme.M.meanScaleVara,mme.M.meanScaleVara2 = 0.0,0.0          #scale parameter for prior of marker effect variance
    end
    #phenotypes corrected for all effects
    ycorr       = vec(Matrix(mme.ySparse)-mme.X*sol)
    if mme.M != 0 && α!=zero(α)
        ycorr = ycorr - mme.M.genotypes*α
    end
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
        # 0. Categorical traits (liabilities)
        ########################################################################
        if categorical_trait == true
            ycorr = categorical_trait_sample_liabilities(mme,ycorr,category_obs,threshold)
        end
        ########################################################################
        # 1.1. Non-Marker Location Parameters
        ########################################################################
        ycorr = ycorr + mme.X*sol
        if Rinv == false
            rhs = mme.X'ycorr
        else
            rhs = mme.X'Diagonal(Rinv)*ycorr
        end

        Gibbs(mme.mmeLhs,sol,rhs,mme.RNew)

        ycorr = ycorr - mme.X*sol
        ########################################################################
        # 1.2 Marker Effects
        ########################################################################
        if mme.M !=0
            if methods in ["BayesC","BayesB","BayesA"]
                locus_effect_variances = (methods=="BayesC" ? fill(mme.M.G,mme.M.nMarkers) : mme.M.G)
                nLoci = BayesABC!(mme.M,ycorr,mme.RNew,locus_effect_variances)
            elseif methods=="RR-BLUP"
                BayesC0!(mme.M,ycorr,mme.RNew)
                nLoci = mme.M.nMarkers
            elseif methods == "BayesL"
                BayesL!(mme.M,ycorr,mme.RNew)
                nLoci = mme.M.nMarkers
            elseif methods == "GBLUP"
                ycorr       = ycorr + mme.M.genotypes*mme.M.α
                lhs         = Rinv .+ mme.RNew./(mme.M.G*D)
                mean1       = mme.M.genotypes'*(Rinv.*ycorr)./lhs
                mme.M.α     = mean1 + randn(mme.M.nObs).*sqrt.(mme.RNew./lhs)
                ycorr       = ycorr - mme.M.genotypes*mme.M.α
            end
            #sample Pi
            if estimatePi == true
                mme.M.π = samplePi(nLoci, mme.M.nMarkers)
            end
        end
        ########################################################################
        # 2.1 Genetic Covariance Matrix (Polygenic Effects) (variance.jl)
        # 2.2 varainces for (iid) random effects;not required(empty)=>jump out
        ########################################################################
        if mme.MCMCinfo.estimate_variance == true
            sampleVCs(mme,sol)
            addVinv(mme)
        end
        ########################################################################
        # 2.3 Residual Variance
        ########################################################################
        if categorical_trait == false && mme.MCMCinfo.estimate_variance == true
            mme.ROld = mme.RNew
            mme.RNew = sample_variance(ycorr.* (Rinv!=false ? sqrt.(Rinv) : 1.0), length(ycorr), mme.df.residual, mme.scaleRes)
        end
        ########################################################################
        # 2.4 Marker Effects Variance
        ########################################################################
        if mme.M != 0 && mme.MCMCinfo.estimate_variance == true
            if methods in ["BayesC","RR-BLUP"]
                mme.M.G  = sample_variance(α, nLoci, mme.df.marker, mme.M.scale)
            elseif methods == "BayesB"
                for j=1:mme.M.nMarkers
                    mme.M.G[j] = sample_variance(β[j],1,mme.df.marker, mme.M.scale)
                end
            elseif methods == "BayesL"
                ssq = 0.0
                for i=1:size(α,1)
                    ssq += α[i]^2/mme.M.gammaArray[i]
                end
                mme.M.G = (ssq + mme.df.marker*mme.M.scale)/rand(Chisq(nLoci+mme.df.marker))
               # MH sampler of gammaArray (Appendix C in paper)
                sampleGammaArray!(mme.M.gammaArray,α,mme.M.G)
            elseif methods == "GBLUP"
                mme.M.G  = sample_variance(α./sqrt.(D), mme.M.nObs,mme.df.marker, mme.M.scale)
            else
                error("Sampling of marker effect variances is not available")
            end
        end
        ########################################################################
        # 2.5 Update priors using posteriors (empirical)
        ########################################################################
        if update_priors_frequency !=0 && iter%update_priors_frequency==0
            if mme.M!=0 && methods != "BayesB"
                mme.M.scale   = meanVara*(mme.df.marker-2)/mme.df.marker
            end
            if mme.pedTrmVec != 0
                mme.scalePed  = G0Mean*(mme.df.polygenic - size(mme.pedTrmVec,1) - 1)
            end
            mme.scaleRes  =  meanVare*(mme.df.residual-2)/mme.df.residual
            println("\n Update priors from posteriors.")
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
            if mme.M != 0
                output_MCMC_samples(mme,sol,mme.RNew,(mme.pedTrmVec!=0 ? inv(mme.Gi) : false),mme.M.π,mme.M.α,mme.M.G,outfile)
            else
                output_MCMC_samples(mme,sol,mme.RNew,(mme.pedTrmVec!=0 ? inv(mme.Gi) : false),false,false,false,outfile)
            end

            nsamples = (iter-burnin)/output_samples_frequency
            solMean   += (sol - solMean)/nsamples
            solMean2  += (sol .^2 - solMean2)/nsamples
            meanVare  += (mme.RNew - meanVare)/nsamples
            meanVare2 += (mme.RNew .^2 - meanVare2)/nsamples

            if mme.pedTrmVec != 0
                G0Mean  += (inv(mme.Gi)  - G0Mean )/nsamples
                G0Mean2 += (inv(mme.Gi) .^2  - G0Mean2 )/nsamples
            end
            if mme.M != 0
                mme.M.meanAlpha  += (mme.M.α - mme.M.meanAlpha)/nsamples
                mme.M.meanAlpha2 += (mme.M.α .^2 - mme.M.meanAlpha2)/nsamples
                mme.M.meanDelta  += (mme.M.δ - mme.M.meanDelta)/nsamples

                if estimatePi == true
                    mme.M.mean_pi += (mme.M.π-mme.M.mean_pi)/nsamples
                    mme.M.mean_pi2 += (mme.M.π .^2-mme.M.mean_pi2)/nsamples
                end
                if methods != "BayesB"
                    mme.M.meanVara += (mme.M.G - mme.M.meanVara)/nsamples
                    mme.M.meanVara2 += (mme.M.G .^2 - mme.M.meanVara2)/nsamples
                end
                if estimateScale == true
                    mme.M.meanScaleVara += (mme.M.scale - mme.M.meanScaleVara)/nsamples
                    mme.M.meanScaleVara2 += (mme.M.scale .^2 - mme.M.meanScaleVara2)/nsamples
                end
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
    if methods == "GBLUP"
        mv(output_file*"_marker_effects_variances.txt",output_file*"_genetic_variance(REML).txt")
    end
    output=output_result(mme,output_file,
                         solMean,meanVare,
                         mme.pedTrmVec!=0 ? G0Mean : false,
                         mme.M != 0 ? mme.M.meanAlpha : false,
                         mme.M != 0 ? mme.M.meanDelta : false,
                         mme.M != 0 ? mme.M.meanVara : false,
                         mme.M != 0 ? estimatePi : false,
                         mme.M != 0 ? mme.M.mean_pi : false,
                         mme.M != 0 ? estimateScale : false,
                         mme.M != 0 ? mme.M.meanScaleVara : false,
                         solMean2,meanVare2,
                         mme.pedTrmVec!=0 ? mme.M.G0Mean2 : false,
                         mme.M != 0 ? mme.M.meanAlpha2 : false,
                         mme.M != 0 ? mme.M.meanVara2 : false,
                         mme.M != 0 ? mme.M.mean_pi2 : false,
                         mme.M != 0 ? mme.M.meanScaleVara2 : false)
    return output
end
