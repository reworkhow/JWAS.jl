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

    ntraits = mme.nModels
    nObs    = div(length(resVec),ntraits)
    Ri      = Diagonal(Rinv)
    for traiti = 1:ntraits
        startPosi = (traiti-1)*nObs + 1
        endPosi   = startPosi + nObs - 1
        for traitj = traiti:ntraits
            startPosj = (traitj-1)*nObs + 1
            endPosj   = startPosj + nObs - 1
            SRes[traiti,traitj] = resVec[startPosi:endPosi]'*Ri*resVec[startPosj:endPosj]
            SRes[traitj,traiti] = SRes[traiti,traitj]
        end
    end
    if constraint == false
        mme.R  = rand(InverseWishart(νR0 + nObs, convert(Array,Symmetric(PRes + SRes))))
    else     #for constraint R, chisq
        ν        = νR0 - ntraits
        scaleRes = diag(mme.scaleRes/(νR0 - ntraits - 1))*(ν-2)/ν #diag(R_prior)*(ν-2)/ν
        for traiti = 1:ntraits
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


########################################################################
# Marker Effects Variance
########################################################################
function sample_marker_effect_variance(Mi)
    if Mi.ntraits == 1
        if Mi.method in ["BayesC","RR-BLUP"]
            Mi.G  = sample_variance(Mi.α[1], Mi.nLoci, Mi.df, Mi.scale)
        elseif Mi.method == "BayesB"
            for j=1:Mi.nMarkers
                Mi.G[j] = sample_variance(Mi.β[1][j],1,Mi.df, Mi.scale)
            end
        elseif Mi.method == "BayesL"
            ssq = 0.0
            for i=1:size(Mi.α[1],1)
                ssq += Mi.α[1][i]^2/Mi.gammaArray[i]
            end
            Mi.G = (ssq + Mi.df*Mi.scale)/rand(Chisq(Mi.nLoci+Mi.df))
            # MH sampler of gammaArray (Appendix C in paper)
            sampleGammaArray!(Mi.gammaArray,Mi.α[1],Mi.G)
        elseif Mi.method == "GBLUP"
            Mi.G  = sample_variance(Mi.α[1]./sqrt.(Mi.D), Mi.nObs,Mi.df, Mi.scale)
        end
    else
        SM    = zero(Mi.scale)
        if Mi.method == "BayesC"
            for traiti = 1:Mi.ntraits
                for traitj = traiti:Mi.ntraits
                    SM[traiti,traitj]   = (Mi.β[traiti]'Mi.β[traitj])
                    SM[traitj,traiti]   = SM[traiti,traitj]
                end
            end
            Mi.G = rand(InverseWishart(Mi.df + Mi.nMarkers, convert(Array,Symmetric(Mi.scale + SM))))
        elseif Mi.method == "RR-BLUP"
            for traiti = 1:Mi.ntraits
                for traitj = traiti:Mi.ntraits
                    SM[traiti,traitj]   = (Mi.α[traiti]'Mi.α[traitj])
                    SM[traitj,traiti]   = SM[traiti,traitj]
                end
            end
            Mi.G = rand(InverseWishart(Mi.df + Mi.nMarkers, convert(Array,Symmetric(Mi.scale + SM))))
        elseif Mi.method == "BayesL"
            for traiti = 1:Mi.ntraits
                alphai = Mi.α[traiti]./Mi.gammaArray
                for traitj = traiti:Mi.ntraits
                    SM[traiti,traitj]   = (alphai'Mi.α[traitj])
                    SM[traitj,traiti]   = SM[traiti,traitj]
                end
            end
            Mi.G = rand(InverseWishart(Mi.df + Mi.nMarkers, convert(Array,Symmetric(Mi.scale + SM))))
            sampleGammaArray!(Mi.gammaArray,Mi.α,Mi.G)# MH sampler of gammaArray (Appendix C in paper)
        elseif Mi.method == "BayesB"
            marker_effects_matrix = Mi.β[1]
            for traiti = 2:Mi.ntraits
                marker_effects_matrix = [marker_effects_matrix Mi.β[traiti]]
            end
            marker_effects_matrix = marker_effects_matrix'
            beta2 = [marker_effects_matrix[:,i]*marker_effects_matrix[:,i]' for i=1:size(marker_effects_matrix,2)]
            for markeri = 1:Mi.nMarkers
                Mi.G[markeri] = rand(InverseWishart(Mi.df + 1, convert(Array,Symmetric(Mi.scale + beta2[markeri]))))
            end
        elseif Mi.method == "GBLUP"
            for traiti = 1:Mi.ntraits
                for traitj = traiti:Mi.ntraits
                    alphaArrayi         = Mi.α[traiti]
                    alphaArrayj         = Mi.α[traitj]
                    SM[traiti,traitj]   = alphaArrayi'*Diagonal(1 ./Mi.D)*alphaArrayj
                    SM[traitj,traiti]   = SM[traiti,traitj]
                end
            end
            Mi.G = rand(InverseWishart(Mi.df + Mi.nMarkers, convert(Array,Symmetric(Mi.scale + SM))))
        end
    end
end
