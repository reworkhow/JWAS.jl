#  Distributions.jl package always returns Float64
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

function sample_variance(mme,resVec,Rinv;constraint=false)
    νR0     = mme.df.residual
    PRes    = mme.scaleRes
    SRes    = zero(PRes)

    nTraits = mme.nModels
    nObs    = div(length(resVec),nTraits)
    Ri      = Diagonal(Rinv)
    for traiti = 1:nTraits
        startPosi = (traiti-1)*nObs + 1
        endPosi   = startPosi + nObs - 1
        for traitj = traiti:nTraits
            startPosj = (traitj-1)*nObs + 1
            endPosj   = startPosj + nObs - 1
            SRes[traiti,traitj] = resVec[startPosi:endPosi]'*Ri*resVec[startPosj:endPosj]
            SRes[traitj,traiti] = SRes[traiti,traitj]
        end
    end
    if constraint == false
        mme.R  = rand(InverseWishart(νR0 + nObs, convert(Array,Symmetric(PRes + SRes))))
    else     #for constraint R, chisq
        ν        = νR0 - nTraits
        scaleRes = diag(mme.scaleRes/(νR0 - nTraits - 1))*(ν-2)/ν #diag(R_prior)*(ν-2)/ν
        for traiti = 1:nTraits
            mme.R[traiti,traiti]= (SRes[traiti,traiti]+ν*scaleRes[traiti])/rand(Chisq(nObs+ν))
        end
    end
    if mme.MCMCinfo.double_precision == false
        mme.R = Float32.(mme.R)
    end
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
# sample variances for random location effects
################################################################################
function sampleVCs(mme::MME,sol::Union{Array{Float64,1},Array{Float32,1}})
    for random_term in mme.rndTrmVec
      term_array = random_term.term_array
      myI        = SparseMatrixCSC{mme.MCMCinfo.double_precision ? Float64 : Float32}(I, mme.modelTermDict[term_array[1]].nLevels, mme.modelTermDict[term_array[1]].nLevels)
      Vi         = (random_term.Vinv!=0) ? random_term.Vinv : myI
      S          = zeros(length(term_array),length(term_array))
      for i = 1:length(term_array)
          termi      = term_array[i]
          randTrmi   = mme.modelTermDict[termi]
          startPosi  = randTrmi.startPos
          endPosi    = startPosi + randTrmi.nLevels - 1
          for j = i:length(term_array)
              termj     = term_array[j]
              randTrmj  = mme.modelTermDict[termj]
              startPosj = randTrmj.startPos
              endPosj   = startPosj + randTrmj.nLevels - 1
              S[i,j]    = sol[startPosi:endPosi]'*Vi*sol[startPosj:endPosj]
              S[j,i]    = S[i,j]
          end
      end
       q  = mme.modelTermDict[term_array[1]].nLevels
       G0 = rand(InverseWishart(random_term.df + q, convert(Array,Symmetric(random_term.scale + S))))
       if mme.MCMCinfo.double_precision == false
           G0 = Float32.(G0)
       end
       random_term.GiOld = copy(random_term.GiNew)
       random_term.GiNew = copy(inv(G0))
       random_term.Gi    = copy(inv(G0))
       if random_term.randomType == "A"
           mme.Gi    = random_term.Gi
       end
    end
end
