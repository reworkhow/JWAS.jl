function MCMC_BayesC(nIter,mme,df,Pi;
                      sol=false,outFreq=100,
                      missing_phenotypes=false,
                      constraint=nothing,
                      output_marker_effects_frequency=0)

    #Pi is of length nTrait^2

    if size(mme.mmeRhs)==()
       getMME(mme,df)
    end

    if sol == false
       sol=zeros(size(mme.mmeLhs,1))
    end #starting value for sol can be provided

    p = size(mme.mmeLhs,1)
    solMean = fill(0.0,p)

    #priors for residual covariance matrix
    ν       = 4
    nObs    = size(df,1)
    nTraits = size(mme.lhsVec,1)
    νR0     = ν + nTraits
    R0      = mme.R
    prior_R = copy(mme.R)
    PRes    = R0*(νR0 - nTraits - 1)
    SRes    = zeros(Float64,nTraits,nTraits)
    R0Mean  = zeros(Float64,nTraits,nTraits)

    #priors for genetic variance matrix
    if mme.ped != 0
        ν         = 4
        pedTrmVec = mme.pedTrmVec
        k         = size(pedTrmVec,1)
        νG0       = ν + k
        G0        = inv(mme.Gi)
        P         = G0*(νG0 - k - 1)
        S         = zeros(Float64,k,k)
        G0Mean    = zeros(Float64,k,k)
    end

    #priors for marker covaraince matrix
    nObs,nMarkers  = size(mme.M.genotypes)
    dfEffectVar    = 4.0
    mGibbs         = GibbsMats(mme.M.genotypes)
    nObs           = mGibbs.nrows
    nMarkers       = mGibbs.ncols
    mArray         = mGibbs.xArray
    mpm            = mGibbs.xpx
    M              = mGibbs.X
    vEff           = mme.M.G
    νGM            = dfEffectVar + nTraits
    PM             = vEff*(νGM - nTraits - 1)
    SM             = zeros(Float64,nTraits,nTraits)
    GMMean         = zeros(Float64,nTraits,nTraits)

    #starting values for marker effects are all zeros
    #starting values for other location parameters are sol

    ycorr          = vec(full(mme.ySparse))
    wArray         = Array(Array{Float64,1},nTraits)
    alphaArray     = Array(Array{Float64,1},nTraits)
    meanAlphaArray = Array(Array{Float64,1},nTraits)
    deltaArray     = Array(Array{Float64,1},nTraits)
    meanDeltaArray = Array(Array{Float64,1},nTraits)
    uArray         = Array(Array{Float64,1},nTraits)
    meanuArray     = Array(Array{Float64,1},nTraits)

    for traiti = 1:nTraits
        startPosi              = (traiti-1)*nObs  + 1
        ptr                    = pointer(ycorr,startPosi)
        wArray[traiti]         = pointer_to_array(ptr,nObs) #ycorr for different traits
                                                            #wArray is list version reference of ycor
        alphaArray[traiti]     = zeros(nMarkers)
        meanAlphaArray[traiti] = zeros(nMarkers)
        #deltaArray[traiti]     = zeros(nMarkers) #starting values for deltaArray should follow staring vlaus for pi
        deltaArray[traiti]     = ones(nMarkers) #starting values for deltaArray should follow staring vlaus for pi
        meanDeltaArray[traiti] = zeros(nMarkers)
        uArray[traiti]         = zeros(nMarkers)
        meanuArray[traiti]     = zeros(nMarkers)
    end

    nLoci = zeros(nTraits)

    BigPi = copy(Pi)
    BigPiMean = copy(Pi)
    for key in keys(BigPiMean)
        BigPiMean[key]=0.0
    end

    if output_marker_effects_frequency != 0  #write samples for marker effects to a txt file
        outfile = Array{IOStream}(nTraits)
        for traiti in 1:nTraits
           outfile[traiti]=open("marker_effects"*"_"*string(mme.lhsVec[traiti])*"_$(now()).txt","w")
        end

        if mme.M.markerID[1]!="NA"
           for traiti in 1:nTraits
              writedlm(outfile[traiti],mme.M.markerID')
           end
        end
     end










    for iter=1:nIter
        #####################################
        #sample non-marker location parameter
        #####################################
        Gibbs(mme.mmeLhs,sol,mme.mmeRhs)
        ycorr[:] = ycorr[:] - mme.X*sol

        #####################################
        #sample marker effects
        #####################################
        iR0,iGM = inv(mme.R),inv(mme.M.G)
        sampleMarkerEffectsBayesC!(mArray,mpm,wArray,
                                   alphaArray,meanAlphaArray,
                                   deltaArray,meanDeltaArray,
                                   uArray,meanuArray,
                                   iR0,iGM,iter,BigPi)

        temp = deltaArray[1]
        for traiti = 2:nTraits
          temp = [temp deltaArray[traiti]]
        end

        iloci = 1
        nLoci_array=zeros(2^nTraits)
        for i in keys(BigPi) #assume order of key won't change
          temp2 = broadcast(-,temp,i')
          nLoci =  sum(mean(abs(temp2),2).==0.0)
          nLoci_array[iloci] = nLoci +1
          iloci = iloci +1
        end
        tempPi = rand(Dirichlet(nLoci_array))

        iloci = 1
        for i in keys(BigPi)
          BigPi[i] = tempPi[iloci]
          BigPiMean[i] += (tempPi[iloci]-BigPiMean[i])/iter
          iloci = iloci +1
        end
        #####################################
        #sample residual covariance matrix
        #AND marker covariance matrix
        #####################################
        resVec = ycorr
        if missing_phenotypes==true
          sampleMissingResiduals(mme,resVec)
        end

        for traiti = 1:nTraits
            startPosi = (traiti-1)*nObs + 1
            endPosi   = startPosi + nObs - 1
            for traitj = traiti:nTraits
                startPosj = (traitj-1)*nObs + 1
                endPosj   = startPosj + nObs - 1
                SRes[traiti,traitj] = (resVec[startPosi:endPosi]'resVec[startPosj:endPosj])[1,1]
                SRes[traitj,traiti] = SRes[traiti,traitj]
                SM[traiti,traitj]   = (alphaArray[traiti]'alphaArray[traitj])[1]
                #  SM[traiti,traitj]   = (uArray[traiti]'uArray[traitj])[1]
                SM[traitj,traiti]   = SM[traiti,traitj]
            end
        end

        mme.M.G = rand(InverseWishart(νGM + nMarkers, PM + SM))
        R0      = rand(InverseWishart(νR0 + nObs, PRes + SRes))
        if contraint != nothing
          R0 = R0.*!constraint+ prior_R.*constraint
        end

        mme.R = R0
        if missing_phenotypes==true
          RiNotUsing   = mkRi(mme,df) #for missing value;updata mme.ResVar
        end

        R0    = mme.R
        Ri    = kron(inv(R0),speye(nObs))

        X          = mme.X
        mme.mmeLhs = X'Ri*X
        ycorr[:]   = ycorr[:] + mme.X*sol
        mme.mmeRhs = mme.X'Ri*ycorr

        #####################################
        #sample genetic variance matrix (polygenic effects)
        #can make this more efficient by taking advantage of symmetry
        #####################################
        if mme.ped != 0
            for (i,trmi) = enumerate(pedTrmVec)
                pedTrmi  = mme.modelTermDict[trmi]
                startPosi  = pedTrmi.startPos
                endPosi    = startPosi + pedTrmi.nLevels - 1
                for (j,trmj) = enumerate(pedTrmVec)
                    pedTrmj  = mme.modelTermDict[trmj]
                    startPosj  = pedTrmj.startPos
                    endPosj    = startPosj + pedTrmj.nLevels - 1
                    S[i,j] = (sol[startPosi:endPosi]'*mme.Ai*sol[startPosj:endPosj])[1,1]
                end
            end
            pedTrm1 = mme.modelTermDict[pedTrmVec[1]]
            q = pedTrm1.nLevels
            G0 = rand(InverseWishart(νG0 + q, P + S))
            mme.Gi = inv(G0)
            addA(mme)

            G0Mean  += (G0  - G0Mean )/iter
        end

        solMean += (sol - solMean)/iter
        R0Mean  += (R0  - R0Mean )/iter
        GMMean  += (mme.M.G  - GMMean)/iter

        if iter%outFreq==0
            println("posterior means at sample: ",iter)
            println("Residual covariance matrix: \n",R0Mean)
            println("Marker effects covariance matrix: \n",GMMean,"\n")
            println("π: \n",BigPiMean)
        end

        if output_marker_effects_frequency != 0  #write samples for marker effects to a txt file
          if iter%output_marker_effects_frequency==0
            for traiti in 1:nTraits
              writedlm(outfile[traiti],uArray[traiti]')
            end
          end
        end
    end

    if output_marker_effects_frequency != 0  #write samples for marker effects to a txt file
        for traiti in 1:nTraits
           close(outfile[traiti])
        end
    end

    output = Dict()
    output["posterior mean of location parameters"]    = [getNames(mme) solMean]
    output["posterior mean of marker effects covariance matrix"]    = GMMean
    output["posterior mean of residual covaraince matrix"]          = R0Mean

    if mme.ped != 0
      output["posterior mean of polygenic effects covariance matrix"] = G0Mean
    end


    if mme.M.markerID[1]!="NA"
      markerout        = []
      for markerArray in meanuArray
        push!(markerout,[mme.M.markerID markerArray])
      end
    else
      markerout        = meanuArray
    end
    output["posterior mean of marker effects"] = markerout
    output["posterior mean of Pi"] = BigPiMean

    return output
end
