#deprected
function make_valpha(G,pi;option="theory") #from marker effect variance with π=0
    if option=="theory"
        offdiag(a,b)=a*b
    elseif option=="max_qtl" #guanrantee posdef
        offdiag(a,b)=max(a,b)
    elseif option=="min_qtl"
        offdiag(a,b)=min(a,b)
    elseif option=="sqrt"  #guarante posdef
        offdiag(a,b)=sqrt(a*b)
    elseif option=="zero" #make off-dignal of maker vaiance matrix zero
        offdiag(a,b)= Inf
    else
        error("unknown option")
    end

    pi_c        = 1 - pi
    mat_pi_c    = diagm(vec(pi_c))
    for i=1:size(mat_pi_c,2)
        for j=1:size(mat_pi_c,2)
            if i!=j
                mat_pi_c[i,j]=offdiag(mat_pi_c[i,i],mat_pi_c[j,j])
            end
        end
    end
    Valpha    = G./mat_pi_c
    println("marker effect variance with π=$pi is positive definite??: ", isposdef(Valpha))
    Valpha
end


############################
##Bugs not found
##BayesC0
##rewrote in another way
############################


function sample_effects_ycorr!(xArray,xpx,yCorr,alphaArray,meanAlpha,invR0,invG0,
                               iIter,traiti)
    # function to sample effects for trait i
    # invR0 is inverse of R0[traiti,traiti]
    # invG0 is inverse of G0

    nEffects = length(xArray)
    α = alphaArray[traiti] # vector of effects for traiti across all loci
    nTraits = length(alphaArray)
    β = zeros(nTraits)     # vector of effects for locus j across all traits
                           #2x1 vector
    for j=1:nEffects
        for k=1:nTraits
            β[k] = alphaArray[k][j]
        end
        x        = xArray[j]
        lhs      = xpx[j]*invR0 + invG0[traiti,traiti]
        rhs      = dot(x,yCorr) - (invG0[traiti,:]*β)[1]
        invLhs   = 1.0/lhs
        mean     = invLhs*rhs + α[j]
        oldAlpha = α[j]
        α[j]     = mean + randn()*sqrt(invLhs)
        BLAS.axpy!((oldAlpha-α[j,1])*invR0,x,yCorr)
        meanAlpha[j] += (α[j] - meanAlpha[j])/iIter
    end
end

function initVec(a::Float64,x::Array{Float64,1})
    for v in x
        v = a
    end
end

function sampleMCMC(nIter,mme,df;sol=false,outFreq=100,thin=100)
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
    nObs,nMarkers  = size(mme.M.X)

    if mme.M != 0
        dfEffectVar = 10.0 #10000000.0
        nObs,nMarkers  = size(mme.M.X)

        vEff    = mme.M.G/mme.M.mean2pq
        νGM     = dfEffectVar + nTraits
        PM      = vEff*(νGM - nTraits - 1)
        SM      = zeros(Float64,nTraits,nTraits)
        GMMean  = zeros(Float64,nTraits,nTraits)
        mme.M.G = vEff

        mArray = mme.M.xArray
        mpm    = [dot(mme.M.X[:,i],mme.M.X[:,i]) for i=1:size(mme.M.X,2)]
        M      = mme.M.X
    end


    #staring values for
    #starting values for marker effects are all zeros
    #starting values for other location parameters are sol

    ycorr          = vec(full(mme.ySparse))
    wArray         = Array(Array{Float64,1},nTraits)
    alphaArray     = Array(Array{Float64,1},nTraits)
    meanAlphaArray = Array(Array{Float64,1},nTraits)
    for traiti = 1:nTraits
        startPosi              = (traiti-1)*nObs  + 1
        ptr                    = pointer(ycorr,startPosi)
        wArray[traiti]         = pointer_to_array(ptr,nObs) #ycorr for different traits
                                                            #wArray is list version reference of ycor
        alphaArray[traiti]     = zeros(nMarkers)
        meanAlphaArray[traiti] = zeros(nMarkers)
    end
    yCorr  = zeros(nObs)

    #file1=open("marker_effects_samples_trait1.txt","w")
    #file2=open("marker_effects_samples_trait2.txt","w")

    for iter=1:nIter
        #####################################
        #sample non-marker location parameter
        #####################################
        Gibbs(mme.mmeLhs,sol,mme.mmeRhs)
        ycorr[:] = ycorr[:] - mme.X*sol

        #####################################
        #sample marker effects
        #####################################
        iR0 = inv(mme.R)
        iGM = inv(mme.M.G)
        for traiti = 1:nTraits # sample alpha_i for each trait
            initVec(0.0,yCorr)
            for traitj =1:nTraits # making yCorr for trait i
                yCorr[:] = yCorr[:] + iR0[traiti,traitj]*wArray[traitj]
            end
            #sample marker effetcs for trait
            oldAlpha = copy(alphaArray[traiti])
            sample_effects_ycorr!(mArray,mpm,yCorr,alphaArray,
                                  meanAlphaArray[traiti],iR0[traiti,traiti],
                                  iGM,iter,traiti)
            # update w_i
            wArray[traiti][:] = wArray[traiti][:] + mme.M.X*(oldAlpha - alphaArray[traiti])
        end



        #####################################
        #sample residual covariance matrix
        #AND marker covariance matrix
        #####################################
        resVec = ycorr
        #sampleMissingResiduals(mme,resVec)

        for traiti = 1:nTraits
            startPosi = (traiti-1)*nObs + 1
            endPosi   = startPosi + nObs - 1
            for traitj = traiti:nTraits
                startPosj = (traitj-1)*nObs + 1
                endPosj   = startPosj + nObs - 1
                SRes[traiti,traitj] = (resVec[startPosi:endPosi]'resVec[startPosj:endPosj])[1,1]
                SRes[traitj,traiti] = SRes[traiti,traitj]
                SM[traiti,traitj]   = (alphaArray[traiti]'alphaArray[traitj])[1]
                SM[traitj,traiti]   = SM[traiti,traitj]
            end
        end

        #SRes[1,2]=SRes[2,1]=0.0
        #SM[1,2]=SM[2,1]=0.0

        #mme.M.G = rand(InverseWishart(νGM + nMarkers, PM + SM))
        #temp  = rand(InverseWishart(νGM + nMarkers, PM + SM))

        #R0    = rand(InverseWishart(νR0 + nObs, PRes + SRes))
        #mme.R = R0
        #Ri   = mkRi(mme,df) #for missing value
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
        #GMMean  += (temp  - GMMean)/iter


        if iter%outFreq==0
            println("at sample: ",iter)
            println(PRes+SRes)
            println(PM+SM)
            #println(R0Mean)
            #println("SM=",SM/nMarkers)
            #println("temp=",temp)
        end

        #if iter%thin==0
        #    writedlm(file1,meanAlphaArray[1]')
        #    writedlm(file2,meanAlphaArray[2]')
        #end
    end

    #close(file1)
    #close(file2)


    output = Dict()
    output["posteriorMeanLocationParms"] = solMean
    output["posteriorMeanAlpha"]         = meanAlphaArray
    #output["posteriorMeanG0"]            = G0Mean
    output["posteriorMeanR0"]            = R0Mean
    output["posteriorMeanGM"]            = GMMean
    return output
end
