function MT_MCMC_BayesianAlphabet(nIter,mme,df;
                        Rinv                       = false,
                        burnin                     = 0,
                        Pi                         = 0.0,
                        estimatePi                 = false,
                        estimate_variance          = true,
                        estimateScale              = false,
                        sol                        = false,
                        outFreq                    = 1000,
                        output_samples_frequency   = 0,
                        methods                    = "conventional (no markers)",
                        missing_phenotypes         = false,
                        constraint                 = false,
                        update_priors_frequency    = 0,
                        output_file                = "MCMC_samples",
                        causal_structure           = false)

    ############################################################################
    # Pre-Check
    ############################################################################
    #starting values for location parameters(no marker) are sol
    sol,α       = sol[1:size(mme.mmeLhs,1)],sol[(size(mme.mmeLhs,1)+1):end]
    solMean, solMean2  = zero(sol),zero(sol)

    if mme.M != 0
        BigPi,BigPiMean,BigPiMean2 = copy(Pi),copy(Pi),copy(Pi)
        if estimatePi == true
          for key in keys(BigPiMean)
            BigPiMean[key]=0.0
            BigPiMean2[key]=0.0
          end
        end
    end
    #if methods == "BayesCC"  labels,BigPi,BigPiMean=setPi(Pi)  end
    ############################################################################
    # PRIORS
    ############################################################################
    ##Phenotypes adjusted for all effects
    ycorr       = vec(Matrix(mme.ySparse)-mme.X*sol)
    if mme.M != 0 && α!=zero(α)
        ycorr      = ycorr - mme.M.genotypes*α
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

    nTraits        = mme.nModels
    nObs           = div(length(ycorr),nTraits)

    #save posterior mean for residual variance
    meanVare  = zero(mme.R)
    meanVare2 = zero(mme.R)
    #save poseterior mean for variances explained by polygenic effects (A) e.g Animal+ Maternal
    if mme.pedTrmVec != 0
        G0Mean,G0Mean2  = zero(mme.Gi),zero(mme.Gi)
    end
    #priors for marker effect variance
    if mme.M != 0
        ########################################################################
        #Priors for marker covaraince matrix
        ########################################################################
        mGibbs                        = GibbsMats(mme.M.genotypes,Rinv)
        nObs,nMarkers, M              = mGibbs.nrows,mGibbs.ncols,mGibbs.X
        mArray,mRinvArray,mpRinvm     = mGibbs.xArray,mGibbs.xRinvArray,mGibbs.xpRinvx
        if methods=="BayesL"# in BayesL mme.M.G is the scale Matrix, Sigma, in MTBayesLasso paper
            mme.M.G /= 4*(nTraits+1)
            mme.M.scale /= 4*(nTraits+1)
        end
        vEff      = mme.M.G
        meanVara  = zero(mme.R) #variable to save variance for marker effect
        meanVara2 = zero(mme.R) #variable to save variance for marker effect
        meanScaleVara  =  zero(mme.R) #variable to save Scale parameter for prior of marker effect variance
        meanScaleVara2 =  zero(mme.R) #variable to save Scale parameter for prior of marker effect variance
        ########################################################################
        ##WORKING VECTORS
        wArray         = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,nTraits)#wArray is list reference of ycor
        ########################################################################
        #Arrays to save solutions for marker effects
        #In RR-BLUP and BayesL, only alphaArray is used.
        #In BayesA, B and C, alphaArray,deltaArray, betaArray are used.
        ########################################################################
        if methods=="BayesB"
            mme.M.G = fill(mme.M.G,nMarkers)#locus specific covaraince matrix
            nLoci   = zeros(nTraits)
        end
        if methods=="BayesL"
           gammaDist  = Gamma((nTraits+1)/2, 8) #8 (1/8): the scale (rate) parameter
           gammaArray = rand(gammaDist,nMarkers)
        end
        if methods=="GBLUP"
            mme.M.genotypes  = mme.M.genotypes ./ sqrt.(2*mme.M.alleleFreq.*(1 .- mme.M.alleleFreq))
            G       = (mme.M.genotypes*mme.M.genotypes'+ I*0.00001)/nMarkers
            eigenG  = eigen(G)
            L       = eigenG.vectors
            D       = eigenG.values
            # α is pseudo marker effects of length #genotyped inds (starting values = L'(starting value for BV)
            nMarkers= size(mme.M.genotypes,1)
            α       = ((α != zero(α)) ? L'α : zeros(nTraits*nMarkers))  #starting values for pseudo marker effect
            #reset parameters in output
            M2   = mme.output_genotypes ./ sqrt.(2*mme.M.alleleFreq.*(1 .- mme.M.alleleFreq))
            M2Mt = M2*mme.M.genotypes'/nMarkers
            mme.output_genotypes = M2Mt*L*Diagonal(1 ./D)
            #reset parameter in mme.M
            mme.M.G         = mme.M.genetic_variance
            mme.M.scale     = mme.M.G*(mme.df.marker-nTraits-1)
            mme.M.markerID  = string.(1:nMarkers)
            mme.M.genotypes = L
        end
        #starting values for marker effects(zeros) and location parameters (sol)
        betaArray       = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,nTraits) #BayesC
        deltaArray      = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,nTraits) #BayesC
        meandeltaArray  = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,nTraits) #BayesC
        alphaArray      = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,nTraits) #BayesC,BayesC0
        meanalphaArray  = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,nTraits) #BayesC
        meanalphaArray2 = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,nTraits) #BayesC
        for traiti = 1:nTraits
            startPosi              = (traiti-1)*nObs  + 1
            ptr                    = pointer(ycorr,startPosi)
            wArray[traiti]         = unsafe_wrap(Array,ptr,nObs)
            betaArray[traiti]      = copy(α[(traiti-1)*nMarkers+1:traiti*nMarkers])
            deltaArray[traiti]     = ones(typeof(α[1]),nMarkers)
            meandeltaArray[traiti] = zero(betaArray[traiti])
            alphaArray[traiti]         = copy(α[(traiti-1)*nMarkers+1:traiti*nMarkers])
            meanalphaArray[traiti]     = zero(alphaArray[traiti])
            meanalphaArray2[traiti]    = zero(alphaArray[traiti])
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
    @showprogress "running MCMC for "*methods*"..." for iter=1:nIter
        ########################################################################
        # 1.1. Non-Marker Location Parameters
        ########################################################################
        Gibbs(mme.mmeLhs,sol,mme.mmeRhs)
        ycorr[:] = ycorr[:] - mme.X*sol
        ########################################################################
        # 1.2 Marker Effects
        ########################################################################
        if mme.M != 0
          #WILL ADD BURIN INSIDE
          if methods in ["BayesC","BayesB","BayesA"]
            locus_effect_variances = (methods=="BayesC" ? fill(mme.M.G,nMarkers) : mme.M.G)
            MTBayesABC!(mArray,mRinvArray,mpRinvm,
                        wArray,betaArray,deltaArray,alphaArray,
                        mme.R,locus_effect_variances,BigPi)
          elseif methods == "RR-BLUP"
            MTBayesC0!(mArray,mRinvArray,mpRinvm,
                       wArray,alphaArray,
                       mme.R,mme.M.G)
          elseif methods == "BayesL"
            MTBayesL!(mArray,mRinvArray,mpRinvm,
                      wArray,alphaArray,gammaArray,
                      mme.R,mme.M.G)
          elseif methods == "GBLUP"
            iR0      = inv(mme.R)
            iGM      = inv(mme.M.G)
            for trait = 1:nTraits
                wArray[trait][:] = wArray[trait] + mme.M.genotypes*alphaArray[trait]
            end
            lhs    = [iR0*Rinv[i] + iGM/D[i] for i=1:length(D)]
            RHS    = (mme.M.genotypes'Diagonal(Rinv))*reshape(ycorr,nObs,nTraits)*iR0 #size nmarkers (=nObs) * nTraits
            rhs    = [RHS[i,:] for i in 1:size(RHS,1)]  #not column major
            Σα     = Symmetric.(inv.(lhs))              #nTrait X nTrait
            μα     = Σα.*rhs
            αs     = rand.(MvNormal.(μα,Σα))
            for markeri = 1:nMarkers
                for traiti = 1:nTraits
                    alphaArray[traiti][markeri] = αs[markeri][traiti]
                end
            end
            for trait = 1:nTraits
                wArray[trait][:] = wArray[trait] - mme.M.genotypes*alphaArray[trait]
            end
          end
          if estimatePi == true
            if methods in ["BayesC","BayesB"]
              samplePi(deltaArray,BigPi,BigPiMean,iter)
            elseif methods == "BayesCC"
              samplePi(deltaArray,BigPi,BigPiMean,iter,labels)
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

        if mme.M != 0 && estimate_variance == true
            SM    = zero(mme.M.scale)
            if methods == "BayesC"
                for traiti = 1:nTraits
                    for traitj = traiti:nTraits
                        SM[traiti,traitj]   = (betaArray[traiti]'betaArray[traitj])
                        SM[traitj,traiti]   = SM[traiti,traitj]
                    end
                end
                mme.M.G = rand(InverseWishart(mme.df.marker + nMarkers, convert(Array,Symmetric(mme.M.scale + SM))))
            elseif methods == "RR-BLUP"
                for traiti = 1:nTraits
                    for traitj = traiti:nTraits
                        SM[traiti,traitj]   = (alphaArray[traiti]'alphaArray[traitj])
                        SM[traitj,traiti]   = SM[traiti,traitj]
                    end
                end
                mme.M.G = rand(InverseWishart(mme.df.marker + nMarkers, convert(Array,Symmetric(mme.M.scale + SM))))
            elseif methods == "BayesL"
                for traiti = 1:nTraits
                    alphai = alphaArray[traiti]./gammaArray
                    for traitj = traiti:nTraits
                        SM[traiti,traitj]   = (alphai'alphaArray[traitj])
                        SM[traitj,traiti]   = SM[traiti,traitj]
                    end
                end
                mme.M.G = rand(InverseWishart(mme.df.marker + nMarkers, convert(Array,Symmetric(mme.M.scale + SM))))
                sampleGammaArray!(gammaArray,alphaArray,mme.M.G)# MH sampler of gammaArray (Appendix C in paper)
            elseif methods == "BayesB"
                marker_effects_matrix = betaArray[1]
                for traiti = 2:nTraits
                    marker_effects_matrix = [marker_effects_matrix betaArray[traiti]]
                end
                marker_effects_matrix = marker_effects_matrix'

                beta2 = [marker_effects_matrix[:,i]*marker_effects_matrix[:,i]'
                         for i=1:size(marker_effects_matrix,2)]
                for markeri = 1:nMarkers
                  mme.M.G[markeri] = rand(InverseWishart(mme.df.marker + 1, convert(Array,Symmetric(mme.M.scale + beta2[markeri]))))
                end
            elseif methods == "GBLUP"
                for traiti = 1:nTraits
                    for traitj = traiti:nTraits
                        alphaArrayi         = alphaArray[traiti]
                        alphaArrayj         = alphaArray[traitj]
                        SM[traiti,traitj]   = alphaArrayi'*Diagonal(1 ./D)*alphaArrayj
                        SM[traitj,traiti]   = SM[traiti,traitj]
                    end
                end
                mme.M.G = rand(InverseWishart(mme.df.marker + nMarkers, convert(Array,Symmetric(mme.M.scale + SM))))
            else
                error("Sampling of marker effect vairances is not available")
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
            if mme.M != 0
                if methods in ["BayesC","RR-BLUP","BayesL","GBLUP"]
                    output_MCMC_samples(mme,sol,mme.R,(mme.pedTrmVec!=0 ? inv(mme.Gi) : false),BigPi,alphaArray,vec(mme.M.G),outfile)
                elseif methods == "BayesB"
                    output_MCMC_samples(mme,sol,mme.R,(mme.pedTrmVec!=0 ? inv(mme.Gi) : false),BigPi,alphaArray,hcat([x for x in mme.M.G]...),outfile)
                end
            else
                output_MCMC_samples(mme,sol,mme.R,(mme.pedTrmVec!=0 ? inv(mme.Gi) : false),false,fill(false,nTraits),false,outfile)
            end
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
                for trait in 1:nTraits
                    meanalphaArray[trait] += (alphaArray[trait] - meanalphaArray[trait])/nsamples
                    meanalphaArray2[trait]+= (alphaArray[trait].^2 - meanalphaArray2[trait])/nsamples
                    meandeltaArray[trait] += (deltaArray[trait] - meandeltaArray[trait])/nsamples
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
