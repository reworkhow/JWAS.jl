################################################################################
#  SAMPLE VARIANCES FOR MARKER OR RESIDUAL EFFECTS (GIBBS SAMPLER)             #
#*******************************************************************************
#sample from scaled-Inv-χ2 (n+df, (dot(x,x) + df*scale)/(n+df))                *
#Given X ~ Inv-⁠χ2(df) <=> Y = (df*scale)*X ~ Scale-Inv-⁠χ2(df,scale)            *
#Given X ~ Inv-⁠χ2(df) <=> Y = 1/X ~ χ2(df)                                     *
#These transformations are indicated by names already                          *
#*******************************************************************************
function sample_variance(x, n, df, scale)
    return (dot(x,x) + df*scale)/rand(Chisq(n+df))
end

################################################################################
#  SAMPLE VARIANCES FOR OTHER RANDOM EFFECTS (GIBBS SAMPLER)                   #
#*******************************************************************************
#Given Scale-Inverse-Wishart(ν,S)=Inverse-Gamma(ν/2,S/2)=Scale-Inv-chi2(ν,S/ν) *
#variances for random effects (non-marker) are sampled from                    *
#Scale-Inverse-Wishart(ν,S) for coding simplicity under all situations         *
#even when the prior is a scalar (e.g., scaled-Inv-χ2 in single-trait PBLUP    *
#and no maternal) by treating it as a 1x1 matrix.                              *
#*******************************************************************************

################################################################################
# sample variances for IID random location effects
################################################################################
function sampleVCs(mme::MME,sol::Array{Float64,1})
    for random_term in mme.rndTrmVec
      term_array = random_term.term_array
      nLevels    = mme.modelTermDict[term_array[1]].nLevels
      Vi         = (random_term.Vinv!=0) ? random_term.Vinv : SparseMatrixCSC{Float64}(I,nLevels,nLevels)
      S          = zeros(length(term_array),length(term_array))
      for (i,termi) = enumerate(term_array)
          randTrmi   = mme.modelTermDict[termi]
          startPosi  = randTrmi.startPos
          endPosi    = startPosi + randTrmi.nLevels - 1
          for (j,termj) in enumerate(term_array)
            randTrmj    = mme.modelTermDict[termj]
            startPosj   = randTrmj.startPos
            endPosj     = startPosj + randTrmj.nLevels - 1
            S[i,j]      = sol[startPosi:endPosi]'*Vi*sol[startPosj:endPosj]
          end
       end
       q  = mme.modelTermDict[term_array[1]].nLevels
       G0 = rand(InverseWishart(random_term.df + q, convert(Array,Symmetric(random_term.scale + S))))

       random_term.GiOld = copy(random_term.GiNew)
       random_term.GiNew = copy(inv(G0))
       random_term.Gi    = copy(inv(G0))
    end
end



################################################################################
# sample variances for Polygenic Effects (Genetic Covariance Matrix)           #
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
            S[i,j]    = sol[startPosi:endPosi]'*mme.Ai*sol[startPosj:endPosj]
        end
    end

    q  = mme.modelTermDict[pedTrmVec[1]].nLevels
    G0 = rand(InverseWishart(νG0 + q, convert(Array,Symmetric(P + S))))

    mme.GiOld = copy(mme.GiNew)
    mme.GiNew = copy(inv(G0))
    mme.Gi    = copy(inv(G0))
    return G0
end
