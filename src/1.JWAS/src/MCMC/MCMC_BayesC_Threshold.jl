################################################################################
#MCMC for RR-BLUP, BayesC, BayesCπ and "conventional (no markers).
################################################################################
function MCMC_BayesC_threshold(nIter,mme,df;
                     burnin                     = 0,
                     π                          = 0.0,
                     estimatePi                 = false,
                     estimateScale              = false,
                     sol                        = false,
                     outFreq                    = 1000,
                     methods                    = "BayesC",
                     output_samples_frequency   = 0,
                     update_priors_frequency    = 0,
                     output_file                = "MCMC_samples")

    ############################################################################
    # Pre-Check
    ############################################################################
    #starting values for location parameters(no marker) are sol
    #starting values for thresholds  -∞ < t1=0 < t2 < ... < t_{#category-1} < +∞
    # where t1=0 (must be fixed to center the distribution) and t_{#category-1}<1.
    #Then liabilities are sampled and stored in mme.ySparse
    categories = copy(mme.ySparse) # categories (1,2,2,3,1...)
    thresholds = collect(range(0, length=length(unique(mme.categories)),
                         stop=1))[1:(end-1)]  #t1,t2,...,t_{#c-1}

    cmean      = mme.X*sol #maker effects all zeros
    for i in 1:length(categories) #1,2,2,3,1...
        cmeani = cmean[i]
        while true
            liability = cmeani + randn()
            if sum(liability .> thresholds)+1 == categories[i]
                mme.ySparse[i] = liability
                break
            end
        end
    end
    solMean     = zero(sol)
    ############################################################################
    # PRIORS
    ############################################################################
    #prior for residual variance
    # vRes        = mme.RNew
    # nuRes       = mme.df.residual
    # scaleRes    = vRes*(nuRes-2)/nuRes
    # meanVare    = 0.0
    vRes = 1.0

    #priors for variances explained by polygenic effects (A) e.g Animal+ Maternal
    if mme.pedTrmVec != 0
       ν         = mme.df.polygenic
       pedTrmVec = mme.pedTrmVec
       k         = size(pedTrmVec,1)  #2
       νG0       = ν + k #final df for this inverse wisahrt
       G0        = inv(mme.GiNew)
       P         = G0*(νG0 - k - 1)
       S         = zeros(Float64,k,k)
       G0Mean    = zeros(Float64,k,k)
    end

    #priors for marker effect variance
    if mme.M != 0
        mGibbs      = GibbsMats(mme.M.genotypes)
        nObs,nMarkers,mArray,mpm,M = mGibbs.nrows,mGibbs.ncols,mGibbs.xArray,mGibbs.xpx,mGibbs.X
        dfEffectVar = mme.df.marker
        vEff        = mme.M.G
        scaleVar    = vEff*(dfEffectVar-2)/dfEffectVar #scale factor for locus effects
        #println("Init scaleVar to ",scaleVar)
        meanVara     = 0.0 #variable to save variance for marker effect
        meanScaleVar = 0.0 #variable to save Scale parameter for prior of marker effect variance
        #vectors to save solutions for marker effects
        α           = zeros(nMarkers)#starting values for marker effeccts are zeros
        δ           = zeros(nMarkers)#inclusion indicator for marker effects
        meanAlpha   = zeros(nMarkers)#vectors to save solutions for marker effects
        mean_pi     = 0.0
        if methods=="BayesB" #α=β.*δ
            β              = zeros(nMarkers) ##partial marker effeccts
            locusEffectVar = fill(vEff,nMarkers)
            vEff           = locusEffectVar #vEff is scalar for BayesC but a vector for BayeB
        end
    end

    ############################################################################
    #  WORKING VECTORS (ycor, saving values)
    ############################################################################
    #adjust y for starting values
    ycorr       = vec(Matrix(mme.ySparse)-mme.X*sol)

    ############################################################################
    #  SET UP OUTPUT MCMC samples
    ############################################################################
    if output_samples_frequency != 0
        outfile=output_MCMC_samples_setup(mme,nIter-burnin,output_samples_frequency,output_file)
    end

    ############################################################################
    # MCMC (starting values for sol (zeros);  vRes; G0 are used)
    ############################################################################
    @showprogress "running MCMC for "*methods*"..." for iter=1:nIter

        ########################################################################
        # 0. Categorical traits
        # 0.1 liabilities
        ########################################################################
        cmean = mme.ySparse - ycorr #liabilities - residuals
        for i in 1:length(categories) #1,2,2,3,1...
            cmeani = cmean[i]
            while true
                liability = cmeani + randn()
                if sum(liability .> thresholds)+1 == categories[i]
                    mme.ySparse[i] = liability
                    break
                end
            end
        end

        ########################################################################
        # 0. Categorical traits
        # 0.2 Thresholds
        ########################################################################
        #the first threshold t1=0 must be fixed to center the distribution
        for i in 2:length(thresholds) #e.g., t1, t2 for categories 1,2,3
            lowerboundry  = max(mme.ySparse[mme.categories .== i])
            upperboundry  = min(mme.ySparse[mme.categories .== (i+1)])
            thresholds[i] = rand(Uniform(lowerboundry,upperboundry))
        end
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
        if mme.M !=0
            if methods=="BayesC"
                nLoci = sampleEffectsBayesC!(mArray, mpm, ycorr, α, δ,vRes, vEff, π)
            elseif methods=="RR-BLUP"
                sampleEffectsBayesC0!(mArray,mpm,ycorr,α,vRes,vEff)
                nLoci = nMarkers
            elseif methods=="BayesB"
                nLoci = sampleEffectsBayesB!(mArray,mpm,ycorr,α,β,δ,vRes,locusEffectVar,π)
            end
            if iter > burnin
                meanAlpha += (α - meanAlpha)/(iter-burnin)
            end

            #sample Pi
            if estimatePi == true
                π = samplePi(nLoci, nMarkers)
                if iter > burnin
                    mean_pi += (π-mean_pi)/(iter-burnin)
                end
            end
        end
        ########################################################################
        # 2.1 Genetic Covariance Matrix (Polygenic Effects) (variance.jl)
        ########################################################################
        if mme.pedTrmVec != 0
            G0=sample_variance_pedigree(mme,pedTrmVec,sol,P,S,νG0)
            addA(mme)
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
        # mme.ROld = mme.RNew
        # vRes     = sample_variance(ycorr, length(ycorr), nuRes, scaleRes)
        # mme.RNew = vRes
        # if iter > burnin
        #     meanVare += (vRes - meanVare)/(iter-burnin)
        # end
        ########################################################################
        # 2.4 Marker Effects Variance
        ########################################################################
        if mme.M != 0
            if methods != "BayesB" && mme.MCMCinfo.estimate_variance
                vEff  = sample_variance(α, nLoci, dfEffectVar, scaleVar)
            elseif methods == "BayesB"
                for j=1:nMarkers
                    vEff[j] = sample_variance(β[j],1,dfEffectVar, scaleVar)
                end
            end
            if iter > burnin && methods != "BayesB"
                meanVara += (vEff - meanVara)/(iter-burnin)
            end
        end

        ########################################################################
        # 2.5 Update priors using posteriors (empirical)
        ########################################################################
        if update_priors_frequency !=0 && iter%update_priors_frequency==0
            if mme.M!=0 && methods != "BayesB"
                scaleVar    = meanVara*(dfEffectVar-2)/dfEffectVar #scale factor for locus effects
            end
            if mme.pedTrmVec != 0
                P  = G0Mean*(νG0 - k - 1)
            end
            scaleRes    =  meanVare*(nuRes-2)/nuRes
            println("\n Update priors from posteriors.")
        end

        ########################################################################
        # 2.6 Estimate Scale parameter in prior for variance of marker effects
        ########################################################################
        if estimateScale == true
            a = size(vEff,1)*dfEffectVar/2   + 1
            b = sum(dfEffectVar ./ (2vEff )) + 1
            scaleVar = rand(Gamma(a,1/b))
            #println("scaleVar = ",scaleVar)
            if iter > burnin
                meanScaleVar += (scaleVar - meanScaleVar)/(iter-burnin)
            end
        end
        ########################################################################
        # 3.1 Save MCMC samples
        ########################################################################
        if output_samples_frequency != 0 && (iter-burnin)%output_samples_frequency==0 && iter>burnin
            if mme.M != 0
                output_MCMC_samples(mme,sol,vRes,(mme.pedTrmVec!=0 ? G0 : false),π,α,vEff,outfile)
            else
                output_MCMC_samples(mme,sol,vRes,(mme.pedTrmVec!=0 ? G0 : false),false,false,false,outfile)
            end
        end

        ########################################################################
        # 3.2 Printout
        ########################################################################
        if iter%outFreq==0 && iter>burnin
            println("\nPosterior means at iteration: ",iter)
            println("Residual variance: ",round(meanVare,digits=6))
            if mme.pedTrmVec !=0
                println("Polygenic effects covariance matrix \n",round.(G0Mean,digits=3))
            end
            if mme.M != 0
                if methods!="BayesB"
                    println("Marker effects variance: ",round(meanVara,digits=6))
                end
                if estimatePi == true
                    println("π: ", round(mean_pi,digits=3))
                end
            end
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
    if mme.M != 0
        output=output_result(mme,solMean,meanVare,(mme.pedTrmVec!=0 ? G0Mean : false),output_samples_frequency,
                             meanAlpha,meanVara,estimatePi,mean_pi,estimateScale,meanScaleVar,output_file)
    else
        output=output_result(mme,solMean,meanVare,(mme.pedTrmVec!=0 ? G0Mean : false),output_samples_frequency,
                             false,false,false,false,output_file)
    end
    return output
end
