function MT_MCMC_BayesianAlphabet(nIter,mme,df;
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
        for key in keys(BigPiMean)
          BigPiMean[key]=0.0
          BigPiMean2[key]=0.0
        end
    end
    #if methods == "BayesCC"  labels,BigPi,BigPiMean=setPi(Pi)  end
    ############################################################################
    # PRIORS
    ############################################################################
    ##Phenotypes adjusted for all effects (assuming zeros now)
    ycorr          = vec(Matrix(mme.ySparse))
    nTraits        = mme.nModels
    nObs           = div(length(ycorr),nTraits)
    if missing_phenotypes==true
        RiNotUsing   = mkRi(mme,df) #fill up missing phenotypes patterns
    end

    #save posterior mean for residual variance
    meanVare  = zeros(Float64,nTraits,nTraits)
    meanVare2 = zeros(Float64,nTraits,nTraits)
    #save poseterior mean for variances explained by polygenic effects (A) e.g Animal+ Maternal
    if mme.pedTrmVec != 0
        G0Mean,G0Mean2  = zero(mme.Gi),zero(mme.Gi)
    end
    #priors for marker effect variance
    if mme.M != 0
        ########################################################################
        #Priors for marker covaraince matrix
        ########################################################################
        mGibbs      = GibbsMats(mme.M.genotypes)
        nObs,nMarkers,mArray,mpm,M = mGibbs.nrows,mGibbs.ncols,mGibbs.xArray,mGibbs.xpx,mGibbs.X
        dfEffectVar = mme.df.marker
        if methods=="BayesL"# in BayesL mme.M.G is the scale Matrix, Sigma, in MTBayesLasso paper
            mme.M.G /= 4*(nTraits+1)
            mme.M.scale /= 4*(nTraits+1)
        end
        vEff        = mme.M.G
        meanVara  = zeros(Float64,nTraits,nTraits) #variable to save variance for marker effect
        meanVara2 = zeros(Float64,nTraits,nTraits) #variable to save variance for marker effect
        meanScaleVara  =  zeros(Float64,nTraits,nTraits) #variable to save Scale parameter for prior of marker effect variance
        meanScaleVara2 =  zeros(Float64,nTraits,nTraits) #variable to save Scale parameter for prior of marker effect variance
        ########################################################################
        ##WORKING VECTORS
        wArray         = Array{Array{Float64,1}}(undef,nTraits)#wArray is list reference of ycor
        ########################################################################
        #Arrays to save solutions for marker effects
        #In RR-BLUP and BayesL, only alphaArray is used.
        #In BayesA, B and C, alphaArray,deltaArray, betaArray are used.
        ########################################################################
        #starting values for marker effects(zeros) and location parameters (sol)
        betaArray     = Array{Array{Float64,1}}(undef,nTraits) #BayesC
        deltaArray     = Array{Array{Float64,1}}(undef,nTraits) #BayesC
        meandeltaArray = Array{Array{Float64,1}}(undef,nTraits) #BayesC
        alphaArray         = Array{Array{Float64,1}}(undef,nTraits) #BayesC,BayesC0
        meanalphaArray     = Array{Array{Float64,1}}(undef,nTraits) #BayesC
        meanalphaArray2     = Array{Array{Float64,1}}(undef,nTraits) #BayesC

        for traiti = 1:nTraits
            startPosi              = (traiti-1)*nObs  + 1
            ptr                    = pointer(ycorr,startPosi)
            wArray[traiti]         = unsafe_wrap(Array,ptr,nObs)
            betaArray[traiti]      = copy(α[(traiti-1)*nMarkers+1:traiti*nMarkers])
            deltaArray[traiti]     = ones(nMarkers)
            meandeltaArray[traiti] = zeros(nMarkers)
            alphaArray[traiti]         = copy(α[(traiti-1)*nMarkers+1:traiti*nMarkers])
            meanalphaArray[traiti]     = zeros(nMarkers)
            meanalphaArray2[traiti]    = zeros(nMarkers)
        end

        if methods=="BayesB"
            arrayG  = fill(mme.M.G,nMarkers)#locus specific covaraince matrix
            mme.M.G = arrayG
            nLoci   = zeros(nTraits)
        end
        if methods=="BayesL"
           gammaDist  = Gamma((nTraits+1)/2, 8) # 8 is the scale parameter of the Gamma distribution (1/8 is the rate parameter)
           gammaArray = rand(gammaDist,nMarkers)
        end
    end

    ############################################################################
    # Starting values for SEM
    ############################################################################
    if causal_structure != false #starting values
        #starting values for 1) structural coefficient λij (i≂̸j) is zero,thus
        #now Λycorr = ycorr, later variable "ycorr" actually denotes Λycorr for
        #coding convinience 2)no missing phenotypes
        Y  = get_sparse_Y_FRM(wArray,causal_structure) #here wArray is for phenotypes (before corrected)
        Λ  = Matrix{Float64}(I,nTraits,nTraits) #structural coefficient λij (i≂̸j) is zero
        Λy = kron(Λ,sparse(1.0I,nObs,nObs))*mme.ySparse

        causal_structure_filename = "strcuture_coefficient_MCMC_samples.txt"
        causal_structure_outfile  = open(causal_structure_filename,"w")   #write MCMC samples for Λ to a txt file
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
          iR0      = inv(mme.R)
          iGM      = (methods!="BayesB" ? inv(mme.M.G) : inv.(mme.M.G)) #an array in BayesB
          #WILL ADD BURIN INSIDE
          if methods == "BayesC"
            sampleMarkerEffectsBayesC!(mArray,mpm,wArray,
                                      betaArray,
                                      deltaArray,
                                      alphaArray,
                                      iR0,iGM,BigPi)
          elseif methods == "RR-BLUP"
            sampleMarkerEffectsMTBayesC0!(mArray,mpm,wArray,
                                          alphaArray,
                                          iR0,iGM)
          elseif methods == "BayesB"
            sampleMarkerEffectsBayesB!(mArray,mpm,wArray,
                                       betaArray,
                                       deltaArray,
                                       alphaArray,
                                       iR0,iGM,BigPi)
          elseif methods == "BayesL"
            sampleMarkerEffectsMTBayesL!(mArray,mpm,wArray,
                                         alphaArray,
                                         gammaArray,
                                         iR0,iGM)
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
            sample_variance(mme,resVec,constraint=constraint)
        end
        Ri = kron(inv(mme.R),SparseMatrixCSC{Float64}(I, nObs, nObs))
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
        ########################################################################
        if mme.pedTrmVec != 0 && estimate_variance == true
            sample_variance_pedigree(mme,sol)
            addA(mme)
        end
        ########################################################################
        # 2.2 varainces for (iid) random effects;not required(empty)=>jump out
        ########################################################################
        if estimate_variance == true
            sampleVCs(mme,sol)
            addLambdas(mme)
        end
        ########################################################################
        # 2.3 Marker Covariance Matrix
        ########################################################################

        if mme.M != 0 && estimate_variance == true
            SM    = zero(mme.M.scale)
            if !(methods in ["BayesB","BayesL","RR-BLUP"])
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
            else
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
                if methods in ["BayesC","RR-BLUP","BayesL"]
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

            if mme.pedTrmVec !=0
              println("Polygenic effects covariance matrix \n",round.(G0Mean,digits=6))
            end

            if mme.M != 0
              if methods != "BayesB"
                  println("Marker effects covariance matrix: \n",round.(meanVara,digits=6))
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


function sampleGammaArray!(gammaArray,alphaArray,mmeMG)
    Gi = inv(mmeMG)
    nMarkers = size(gammaArray,1)
    nTraits  = length(alphaArray[1])==1 ? 1 : length(alphaArray)

    Q  = zeros(nMarkers)
    nTraits > 1 ? calcMTQ!(Q,nMarkers,nTraits,alphaArray,Gi) : calcSTQ!(Q,nMarkers,alphaArray,Gi)
    gammaDist = Gamma(0.5,4) # 4 is the scale parameter, which corresponds to a rate parameter of 1/4
    candidateArray = 1 ./ rand(gammaDist,nMarkers)
    uniformArray = rand(nMarkers)
    acceptProbArray = exp.(Q ./4 .*(2 ./ gammaArray - candidateArray))
    replace = uniformArray .< acceptProbArray
    gammaArray[replace] = 2 ./ candidateArray[replace]
end

function calcMTQ!(Q,nMarkers,nTraits,alphaArray,Gi)
    for locus = 1:nMarkers
        for traiti = 1:nTraits
            for traitj = 1:nTraits
                Q[locus] += alphaArray[traiti][locus]*alphaArray[traitj][locus]*Gi[traiti,traitj]
            end
        end
    end
end

function calcSTQ!(Q,nMarkers,alphaArray,Gi)
        Q .= alphaArray.^2 ./Gi
end
