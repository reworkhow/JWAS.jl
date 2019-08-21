function MT_MCMC_BayesC(nIter,mme,df;
                        burnin                     = 0,
                        Pi                         = 0.0,
                        estimatePi                 = false,
                        estimate_variance          = true,
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
    solMean     = zero(sol)

    if methods in ["RR-BLUP","BayesL"]
        BigPi = copy(Pi) #temporary for output_MCMC_samples function
    elseif methods in ["BayesC","BayesB"]
        BigPi = copy(Pi)
        BigPiMean = copy(Pi)
        for key in keys(BigPiMean)
          BigPiMean[key]=0.0
        end
    elseif methods == "BayesCC"
        if mme.M == 0
            error("BayesCC runs with genotypes")
        end
      #label with the index of Array
      #labels[1]=[0.0,0.0],labels[2]=[1.0,1.0]
      #BigPi[1]=0.3,BigPi[2]=0.7
      nlabels= length(Pi)
      labels = Array(Array{Float64,1},nlabels)
      BigPi  = Array(Float64,nlabels)
      whichlabel=1
      for pair in sort(collect(Pi), by=x->x[2],rev=true)
        key=pair[1]
        labels[whichlabel]=copy(key)
        BigPi[whichlabel]=copy(pair[2])
        whichlabel = whichlabel+1
      end
      BigPiMean = zeros(BigPi)
    end
    ############################################################################
    # PRIORS
    ############################################################################
    #Priors for residual covariance matrix
    ν       = mme.df.residual
    nObs    = size(df,1)
    nTraits = size(mme.lhsVec,1)
    νR0     = ν + nTraits
    R0      = mme.R
    PRes    = R0*(νR0 - nTraits - 1)
    SRes    = zeros(Float64,nTraits,nTraits)
    R0Mean  = zeros(Float64,nTraits,nTraits)
    scaleRes= diag(mme.R)*(ν-2)/ν #only for chisq for constraint diagonal R

    if missing_phenotypes==true
        RiNotUsing   = mkRi(mme,df) #fill up missing phenotypes patterns
    end

    #Priors for polygenic effect covariance matrix
    if mme.pedTrmVec != 0
      ν         = mme.df.polygenic
      pedTrmVec = mme.pedTrmVec
      k         = size(pedTrmVec,1)
      νG0       = ν + k
      G0        = inv(mme.Gi)
      P         = G0*(νG0 - k - 1)
      S         = zeros(Float64,k,k)
      G0Mean    = zeros(Float64,k,k)
    end

    ##SET UP MARKER PART
    if mme.M != 0
        ########################################################################
        #Priors for marker covaraince matrix
        ########################################################################
        mGibbs      = GibbsMats(mme.M.genotypes)
        nObs,nMarkers,mArray,mpm,M = mGibbs.nrows,mGibbs.ncols,mGibbs.xArray,mGibbs.xpx,mGibbs.X
        dfEffectVar = mme.df.marker
        if methods=="BayesL"
            mme.M.G /= 4*(nTraits+1) # in BayesL mme.M.G is the scale Matrix, Sigma, in the MTBayesLasso paper
        end
        vEff        = mme.M.G
        νGM         = dfEffectVar + nTraits
        PM          = vEff*(νGM - nTraits - 1)
        SM          = zeros(Float64,nTraits,nTraits)
        GMMean      = zeros(Float64,nTraits,nTraits)
        ########################################################################
        ##WORKING VECTORS
        ########################################################################
        ycorr          = vec(Matrix(mme.ySparse)) #ycorr for different traits
        wArray         = Array{Array{Float64,1}}(undef,nTraits)#wArray is list reference of ycor
        ########################################################################
        #Arrays to save solutions for marker effects
        ########################################################################
        #starting values for marker effects(zeros) and location parameters (sol)
        alphaArray     = Array{Array{Float64,1}}(undef,nTraits) #BayesC,BayesC0
        meanAlphaArray = Array{Array{Float64,1}}(undef,nTraits) #BayesC,BayesC0
        deltaArray     = Array{Array{Float64,1}}(undef,nTraits) #BayesC
        meanDeltaArray = Array{Array{Float64,1}}(undef,nTraits) #BayesC
        uArray         = Array{Array{Float64,1}}(undef,nTraits) #BayesC
        meanuArray     = Array{Array{Float64,1}}(undef,nTraits) #BayesC

        for traiti = 1:nTraits
            startPosi              = (traiti-1)*nObs  + 1
            ptr                    = pointer(ycorr,startPosi)
            #wArray[traiti]         = pointer_to_array(ptr,nObs)
            wArray[traiti]         = unsafe_wrap(Array,ptr,nObs)

            alphaArray[traiti]     = zeros(nMarkers)
            meanAlphaArray[traiti] = zeros(nMarkers)
            deltaArray[traiti]     = ones(nMarkers)
            meanDeltaArray[traiti] = zeros(nMarkers)
            uArray[traiti]         = zeros(nMarkers)
            meanuArray[traiti]     = zeros(nMarkers)
        end

        if methods=="BayesB"
            arrayG  = Array{Array{Float64,2}}(undef,nMarkers)#locus specific covaraince matrix
            for markeri = 1:nMarkers
              arrayG[markeri]= vEff
            end
            nLoci = zeros(nTraits)
        end
        if methods=="BayesL"
           gammaDist  = Gamma((nTraits+1)/2, 8) # 8 is the scale paramenter of the Gamma distribution (1/8 is the rate parameter)
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

        if iter > burnin
            solMean += (sol - solMean)/(iter-burnin)
        end
        ########################################################################
        # 1.2 Marker Effects
        ########################################################################
        if mme.M != 0
          ycorr[:] = ycorr[:] - mme.X*sol
          iR0      = inv(mme.R)
          iGM      = (methods!="BayesB" ? inv(mme.M.G) : [inv(G) for G in arrayG])
          #WILL ADD BURIN INSIDE
          if methods == "BayesC"
            sampleMarkerEffectsBayesC!(mArray,mpm,wArray,
                                      alphaArray,meanAlphaArray,
                                      deltaArray,meanDeltaArray,
                                      uArray,meanuArray,
                                      iR0,iGM,iter,BigPi,burnin)
          elseif methods == "RR-BLUP"
            sampleMarkerEffectsMTBayesC0!(mArray,mpm,wArray,alphaArray,
                               meanAlphaArray,iR0,iGM,iter,burnin)
          elseif methods == "BayesCC"
            sampleMarkerEffectsBayesCC!(mArray,mpm,wArray,
                                        alphaArray,meanAlphaArray,
                                        deltaArray,meanDeltaArray,
                                        uArray,meanuArray,
                                        iR0,iGM,iter,BigPi,labels,burnin)
          elseif methods == "BayesB"
            sampleMarkerEffectsBayesB!(mArray,mpm,wArray,
                                       alphaArray,meanAlphaArray,
                                       deltaArray,meanDeltaArray,
                                       uArray,meanuArray,
                                       iR0,iGM,iter,BigPi,burnin)
          elseif methods == "BayesL"
            sampleMarkerEffectsMTBayesL!(mArray,mpm,wArray,alphaArray,
                               meanAlphaArray,gammaArray,iR0,iGM,iter,burnin)
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
        resVec = (mme.M==0 ? (mme.ySparse - mme.X*sol) : ycorr)
        #here resVec is alias for ycor ***

        if missing_phenotypes==true
          resVec[:]=sampleMissingResiduals(mme,resVec)
        end

        if estimate_variance == true
            for traiti = 1:nTraits
                startPosi = (traiti-1)*nObs + 1
                endPosi   = startPosi + nObs - 1
                for traitj = traiti:nTraits
                    startPosj = (traitj-1)*nObs + 1
                    endPosj   = startPosj + nObs - 1
                    SRes[traiti,traitj] = (resVec[startPosi:endPosi]'resVec[startPosj:endPosj])[1,1]
                    SRes[traitj,traiti] = SRes[traiti,traitj]
                end
            end

            R0      = rand(InverseWishart(νR0 + nObs, convert(Array,Symmetric(PRes + SRes))))

            #for constraint R, chisq
            if constraint == true
                R0 = zeros(nTraits,nTraits)
                for traiti = 1:nTraits
                    R0[traiti,traiti]= (SRes[traiti,traiti]+ν*scaleRes[traiti])/rand(Chisq(nObs+ν))
                end
            end
        end

        mme.R = R0
        Ri    = kron(inv(R0),SparseMatrixCSC{Float64}(I, nObs, nObs))

        if iter > burnin
            R0Mean  += (R0  - R0Mean )/(iter-burnin)
        end
        ########################################################################
        # -- LHS and RHS for conventional MME (No Markers)
        # -- Position: between new Ri and new Ai
        ########################################################################
        X          = mme.X
        mme.mmeLhs = X'Ri*X #LHS for normal equation (no random effects)
        if mme.M != 0
          ycorr[:]   = ycorr[:] + X*sol
          #same to ycorr[:]=resVec+X*sol
        end
        mme.mmeRhs = (mme.M == 0 ? (X'Ri*mme.ySparse) : (X'Ri*ycorr))
        ########################################################################
        # 2.2 Genetic Covariance Matrix (Polygenic Effects)
        ########################################################################
        if mme.pedTrmVec != 0 && estimate_variance == true #will optimize taking advantage of symmetry
            for (i,trmi) = enumerate(pedTrmVec)
                pedTrmi   = mme.modelTermDict[trmi]
                startPosi = pedTrmi.startPos
                endPosi   = startPosi + pedTrmi.nLevels - 1
                for (j,trmj) = enumerate(pedTrmVec)
                    pedTrmj    = mme.modelTermDict[trmj]
                    startPosj  = pedTrmj.startPos
                    endPosj    = startPosj + pedTrmj.nLevels - 1
                    S[i,j]     = (sol[startPosi:endPosi]'*mme.Ai*sol[startPosj:endPosj])[1,1]
                end
            end
            pedTrm1 = mme.modelTermDict[pedTrmVec[1]]
            q = pedTrm1.nLevels

            G0 = rand(InverseWishart(νG0 + q, convert(Array,Symmetric(P + S))))
            mme.Gi = inv(G0)
            addA(mme)
            if iter > burnin
                G0Mean  += (G0  - G0Mean )/(iter-burnin)
            end
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
            if !(methods in ["BayesB","BayesL"])
                for traiti = 1:nTraits
                    for traitj = traiti:nTraits
                        SM[traiti,traitj]   = (alphaArray[traiti]'alphaArray[traitj])[1]
                        SM[traitj,traiti]   = SM[traiti,traitj]
                    end
                end
                mme.M.G = rand(InverseWishart(νGM + nMarkers, convert(Array,Symmetric(PM + SM))))
                if iter > burnin
                    GMMean  += (mme.M.G  - GMMean)/(iter-burnin)
                end
            elseif methods == "BayesL"
                for traiti = 1:nTraits
                    alphai = alphaArray[traiti]./gammaArray
                    for traitj = traiti:nTraits
                        SM[traiti,traitj]   = (alphai'alphaArray[traitj])[1]
                        SM[traitj,traiti]   = SM[traiti,traitj]
                    end
                end
                mme.M.G = rand(InverseWishart(νGM + nMarkers, convert(Array,Symmetric(PM + SM))))
                if iter > burnin
                    GMMean  += (mme.M.G  - GMMean)/(iter-burnin)
                end
                # MH sampler of gammaArray (Appendix C in paper)
                sampleGammaArray!(gammaArray,alphaArray,mme.M.G)
            else
                marker_effects_matrix = alphaArray[1]
                for traiti = 2:nTraits
                    marker_effects_matrix = [marker_effects_matrix alphaArray[traiti]]
                end
                marker_effects_matrix = marker_effects_matrix'

                alpha2 = [marker_effects_matrix[:,i]*marker_effects_matrix[:,i]'
                         for i=1:size(marker_effects_matrix,2)]

                for markeri = 1:nMarkers
                  arrayG[markeri] = rand(InverseWishart(νGM + 1, convert(Array,Symmetric(PM + alpha2[markeri]))))
                end
            end
        end

        ########################################################################
        # 3.1 Causal Relationships among phenotypes (Structure Equation Model)
        ########################################################################
        if causal_structure != false
            sample4λ = get_Λ(Y,mme.R,ycorr,Λy,mme.ySparse,causal_structure) #no missing phenotypes
        end
        ########################################################################
        # 3.2 Update priors using posteriors (empirical)
        ########################################################################
        if update_priors_frequency !=0 && iter%update_priors_frequency==0
            if mme.M!=0
                PM = GMMean*(νGM - nTraits - 1)
            end
            if mme.pedTrmVec != 0
                P  = G0Mean*(νG0 - k - 1)
            end
            PRes    = R0Mean*(νR0 - nTraits - 1)
            println("\n Update priors from posteriors.")
        end
        ########################################################################
        # OUTPUT
        ########################################################################
        if output_samples_frequency != 0 && (iter-burnin)%output_samples_frequency==0 && iter>burnin
            if mme.M != 0
                if methods in ["BayesC","BayesCC"]
                    output_MCMC_samples(mme,sol,R0,(mme.pedTrmVec!=0 ? G0 : false),BigPi,uArray,vec(mme.M.G),outfile)
                elseif methods in ["RR-BLUP","BayesL"]
                    output_MCMC_samples(mme,sol,R0,(mme.pedTrmVec!=0 ? G0 : false),BigPi,alphaArray,vec(mme.M.G),outfile)
                elseif methods == "BayesB"
                    #res=hcat([x for x in arrayG]...)
                    output_MCMC_samples(mme,sol,R0,(mme.pedTrmVec!=0 ? G0 : false),BigPi,uArray,hcat([x for x in arrayG]...),outfile)
                end
            else
                output_MCMC_samples(mme,sol,R0,(mme.pedTrmVec!=0 ? G0 : false),false,fill(false,nTraits),false,outfile)
            end
            if causal_structure != false
                writedlm(causal_structure_outfile,sample4λ',',')
            end
        end

        if iter%outFreq==0 && iter>burnin
            println("\nPosterior means at iteration: ",iter)
            println("Residual covariance matrix: \n",round.(R0Mean,digits=6))

            if mme.pedTrmVec !=0
              println("Polygenic effects covariance matrix \n",round.(G0Mean,digits=6))
            end

            if mme.M != 0
              if methods != "BayesB"
                  println("Marker effects covariance matrix: \n",round.(GMMean,digits=6))
              end
              if methods in ["BayesC","BayesB"] && estimatePi == true
                println("π: \n",BigPiMean)
              elseif methods=="BayesCC" && estimatePi == true
                println("π for ", labels)
                println(round.(BigPiMean,3))
              end
            end
            println()
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
    if mme.M != 0 && methods in ["RR-BLUP","BayesL"]
        output=output_result(mme,solMean,R0Mean,(mme.pedTrmVec!=0 ? G0Mean : false),output_samples_frequency,
                             meanAlphaArray,GMMean,estimatePi,false,false,false,output_file)
    elseif mme.M != 0
        output=output_result(mme,solMean,R0Mean,(mme.pedTrmVec!=0 ? G0Mean : false),output_samples_frequency,
                             meanuArray,GMMean,estimatePi,BigPiMean,false,false,output_file)
    else
        output=output_result(mme,solMean,R0Mean,(mme.pedTrmVec!=0 ? G0Mean : false),output_samples_frequency,
                             false,false,false,false,false,false,output_file)
    end

    #OUTPUT Marker Effects
    if mme.M != 0
      if methods != "BayesB"
          output["Posterior mean of marker effects covariance matrix"] = GMMean
      end

      if methods=="BayesC"||methods=="BayesCC"||methods=="BayesB"
        output["Model frequency"] = meanDeltaArray
      end
    end

    return output
end

function sampleGammaArray!(gammaArray,alphaArray,mmeMG)
    Gi = inv(mmeMG)
    nMarkers = size(gammaArray,1)
    nTraits  = length(alphaArray)
    Q  = zeros(nMarkers)
    for locus = 1:nMarkers
        for traiti = 1:nTraits
            for traitj = 1:nTraits
                Q[locus] += alphaArray[traiti][locus]*alphaArray[traitj][locus]*Gi[traiti,traitj]
            end
        end
    end
    gammaDist = Gamma(0.5,4) # 4 is the scale parameter, which corresponds to a rate parameter of 1/4
    candidateArray = 1 ./ rand(gammaDist,nMarkers)
    uniformArray = rand(nMarkers)
    acceptProbArray = exp.(Q ./4 .*(2 ./ gammaArray - candidateArray))
    replace = uniformArray .< acceptProbArray
    gammaArray[replace] = 2 ./ candidateArray[replace]
end
