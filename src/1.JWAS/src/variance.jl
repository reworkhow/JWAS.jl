################################################################################
#sample variance of Marker or Residual effects
################################################################################
function sample_variance(x, n, df, scale)
    return (dot(x,x) + df*scale)/rand(Chisq(n+df))
end

################################################################################
#sample variances for iid Random Location effects
################################################################################
# function sampleVCs(mme::MME,sol::Array{Float64,1})
#     for effect in  mme.rndTrmVec
#         trmi       = effect.term
#         startPosi  = trmi.startPos
#         endPosi    = startPosi + trmi.nLevels - 1
#         x          = sol[startPosi:endPosi]
#         effect.vcOld  = effect.vcNew
#         effect.vcNew  = sample_variance(x,trmi.nLevels, effect.df, effect.scale)
#     end
# end

function sampleVCs(mme::MME,sol::Array{Float64,1})
    for random_term in mme.rndTrmVec
      term_array = random_term.term_array
      Vi         = (random_term.Vinv!=0) ? random_term.Vinv : speye(mme.modelTermDict[term_array[1]].nLevels)
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
            S[i,j]    = sol[startPosi:endPosi]'*mme.Ai*sol[startPosj:endPosj]
        end
    end

    q  = mme.modelTermDict[pedTrmVec[1]].nLevels
    G0 = rand(InverseWishart(νG0 + q, convert(Array,Symmetric(P + S)))) #better invchi when ST-PBLUP (and no maternal...)

    mme.GiOld = copy(mme.GiNew)
    mme.GiNew = copy(inv(G0))
    mme.Gi    = copy(inv(G0))

    addA(mme)

    return G0
end
