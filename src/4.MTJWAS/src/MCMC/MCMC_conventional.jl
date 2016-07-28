function MCMC_conventional(nIter,mme,df;
                           sol=false,outFreq=100,
                           missing_phenotypes=false,
                          constraint=false)
    if size(mme.mmeRhs)==()
       getMME(mme,df)
    end

    p = size(mme.mmeLhs,1)
    if sol == false
       sol=zeros(p)
    end #starting value for sol can be provided

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
    scaleRes= diag(mme.R)*(ν-2)/ν #for chisq for constrant diagonal R

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

    for iter=1:nIter

        #sample non-marker location parameter
        Gibbs(mme.mmeLhs,sol,mme.mmeRhs)

        #sample residual covaraince matrix
        resVec = mme.ySparse - mme.X*sol
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
            end
        end
        R0    = rand(InverseWishart(νR0 + nObs, PRes + SRes))

        #for constraint R, chisq
        if constraint == true
          R0 = zeros(nTraits,nTraits)
          for traiti = 1:nTraits
            R0[i,i]= (SRes[traiti,traiti]+ν*scaleRes[traiti])/rand(Chisq(nObs+ν))
          end
        end

        mme.R = R0
        if missing_phenotypes==true
          RiNotUsing   = mkRi(mme,df) #for missing value;updata mme.ResVar
        end
        R0         = mme.R
        Ri         = kron(inv(R0),speye(nObs))
        X          = mme.X
        mme.mmeLhs = X'Ri*X
        mme.mmeRhs = X'Ri*mme.ySparse

        #sample genetic variance matrix (polygenic effects)
        # can make this more efficient by taking advantage of symmetry
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

        if iter%outFreq==0
            println("\nPosterior means at iteration: ",iter)
            println("Residual covariance matrix: \n",round(R0Mean,3),"\n")
        end
    end
    output = Dict()
    output["posterior mean of location parameters"] = [getNames(mme) solMean]
    output["posterior mean of residual covaraince matrix"] = R0Mean
    if mme.ped != 0
      output["posterior mean of genetic covariance matrix"] = G0Mean
    end

    return output
end
