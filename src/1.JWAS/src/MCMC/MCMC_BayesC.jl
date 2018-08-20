#MCMC for RR-BLUP, BayesC, BayesCpi and "conventional (no markers).
function MCMC_BayesC(nIter,mme,df;
                     burnin                     = 0,
                     π                          = 0.0,
                     estimatePi                 = false,
                     sol                        = false,
                     outFreq                    = 1000,
                     methods                    = "BayesC",
                     output_samples_frequency   = 0,
                     output_file                = "MCMC_samples")

    ############################################################################
    # Pre-Check
    ############################################################################
    #starting values for location parameters(no marker) are sol
    solMean     = zero(sol)

    if methods == "conventional (no markers)"
        if mme.M!=0
            error("Conventional analysis runs without genotypes!")
        elseif estimatePi == true
            error("conventional (no markers) analysis runs with estimatePi = false.")
        end
    elseif methods=="RR-BLUP"
        if mme.M == 0
            error("RR-BLUP runs with genotypes")
        elseif π != 0.0
            error("RR-BLUP runs with π=0.")
        elseif estimatePi == true
            error("RR-BLUP runs with estimatePi = false.")
        end
    elseif methods=="BayesC"
        if mme.M == 0
            error("BayesC runs with genotypes.")
        end
    elseif methods=="BayesB"
        if mme.M==0
            error("BayesB runs with genotypes.")
        end
    end
    ############################################################################
    # PRIORS
    ############################################################################
    #prior for residual variance
    vRes        = mme.RNew
    nuRes       = mme.df.residual
    scaleRes    = vRes*(nuRes-2)/nuRes
    meanVare    = 0.0

    #priors for variances explained by polygenic effects (A) e.g Animal+ Maternal
    if mme.pedTrmVec != 0
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
        nObs,nMarkers,mArray,mpm,M = mGibbs.nrows,mGibbs.ncols,mGibbs.xArray,mGibbs.xpx,mGibbs.X
        dfEffectVar = mme.df.marker
        vEff        = mme.M.G
        scaleVar    = vEff*(dfEffectVar-2)/dfEffectVar #scale factor for locus effects
        meanVara    = 0.0 #variable to save variance for marker effect
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
        out_i,outfile=output_MCMC_samples_setup(mme,nIter-burnin,output_samples_frequency,output_file)
    end

    ############################################################################
    # MCMC (starting values for sol (zeros);  vRes; G0 are used)
    ############################################################################
    @showprogress "running MCMC for "*methods*"..." for iter=1:nIter

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
            G0=sample_variance_pedigree(mme,pedTrmVec,sol,P,S,νG0) #better add A outside
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
        vRes     = sample_variance(ycorr, length(ycorr), nuRes, scaleRes)
        mme.RNew = vRes
        if iter > burnin
            meanVare += (vRes - meanVare)/(iter-burnin)
        end
        ########################################################################
        # 2.4 Marker Effects Variance
        ########################################################################
        if mme.M != 0
            if methods != "BayesB"
                vEff  = sample_variance(α, nLoci, dfEffectVar, scaleVar)
            else
                for j=1:nMarkers
                    vEff[j] = sample_variance(α[j],1,dfEffectVar, scaleVar)
                end
            end

            if iter > burnin
                meanVara += (vEff - meanVara)/(iter-burnin)
            end
        end

        ########################################################################
        # 3.1 Save MCMC samples
        ########################################################################
        if output_samples_frequency != 0 && (iter-burnin)%output_samples_frequency==0 && iter>burnin
            if mme.M != 0
                out_i=output_MCMC_samples(mme,out_i,sol,vRes,(mme.pedTrmVec!=0 ? G0 : false),π,α,vEff,outfile)
            else
                out_i=output_MCMC_samples(mme,out_i,sol,vRes,(mme.pedTrmVec!=0 ? G0 : false),false,false,false,outfile)
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
                             meanAlpha,meanVara,estimatePi,mean_pi,output_file)
    else
        output=output_result(mme,solMean,meanVare,(mme.pedTrmVec!=0 ? G0Mean : false),output_samples_frequency,
                             false,false,false,false,output_file)
    end
    return output
end
