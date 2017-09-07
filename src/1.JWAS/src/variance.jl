################################################################################
#sample variance of Marker or Residual effects
################################################################################
function sample_variance(x, n, df, scale)
    return (dot(x,x) + df*scale)/rand(Chisq(n+df))
end

################################################################################
#sample variances for iid Random Location effects
################################################################################
function sampleVCs(mme::MME,sol::Array{Float64,1})
    for effect in  mme.rndTrmVec
        trmi       = effect.term
        startPosi  = trmi.startPos
        endPosi    = startPosi + trmi.nLevels - 1
        x          = sol[startPosi:endPosi]
        effect.vcOld  = effect.vcNew
        effect.vcNew  = sample_variance(x,trmi.nLevels, effect.df, effect.scale)
    end
end

################################################################################
# sample Genetic Covariance Matrix (Polygenic Effects)
################################################################################
function sample_variance_pedigree(mme,pedTrmVec,sol,P,S,νG0)
    for (i,trmi) = enumerate(pedTrmVec)
        pedTrmi   = mme.modelTermDict[trmi]
        startPosi = pedTrmi.startPos
        endPosi   = startPosi + pedTrmi.nLevels - 1
        for (j,trmj) = enumerate(pedTrmVec)
            pedTrmj   = mme.modelTermDict[trmj]
            startPosj = pedTrmj.startPos
            endPosj   = startPosj + pedTrmj.nLevels - 1
            S[i,j]    = (sol[startPosi:endPosi]'*mme.Ai*sol[startPosj:endPosj])[1,1]
        end
    end

    pedTrm1 = mme.modelTermDict[pedTrmVec[1]]
    q  = pedTrm1.nLevels
    G0 = rand(InverseWishart(νG0 + q, convert(Array,Symmetric(P + S)))) #better invchi when ST-PBLUP (and no maternal...)

    mme.GiOld = copy(mme.GiNew)
    mme.GiNew = inv(G0)
    addA(mme) #add Ainverse*lambda

    return G0
end
