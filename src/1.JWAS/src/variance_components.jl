################################################################################
# SAMPLE VARIANCES of MARKER, RESIDUAL and other Random EFFECTS (GIBBS SAMPLER)#
#*******************************************************************************
#singel-trait i.i.d (residual variance, marker effect variance)
#*******************************************************************************
#sample from Scale-Inv-chi2(ν,SSE/ν)  which is a conjugate prior,
#1st parameter: degree of freedom, 2nd parameter: scale,
#ν is (posterior) degree of freedom, and SSE is (posterior) sum of square (SSE).
#
# ν=n+df, SSE = dot(x,x) + df*scale (in conditional posterior)
# x    : data, e.g., residuals (ycorr) or marker effects (α or β)
# n    : length(x)
# df   : degree of freedom for the prior
# scale: scale parameters for the prior (different from scale in InverseWishart)
#
#How to sample from scale-inv-⁠χ2(ν,scale) with scale = SSE/ν
# (1) X ~ Scale-Inv-⁠χ2(ν,SSE/ν)
# <=> X = (df*(SSE/ν))*Y = SSE*Y and Y ~ Inv-⁠χ2(ν)
# (2) Y ~ Inv-⁠χ2(ν)
# <=> Y = 1/Z and Z ~ χ2(ν)
# Thus, X = SSE/Z and Z ~ χ2(ν)
#
#*******************************************************************************
#multi-trait i.i.d (residual variance)
#*******************************************************************************
#sample from InverseWishart(ν,SSE), which is a conjugate prior,
#1st parameter: degree of freedom, 2nd parameter: scale,
#ν is (posterior) degree of freedom, and SSE is (posterior) sum of square (SSE).
#
# ν = n+df, SSE = scale + SSE* (in conditional posterior)
# SSE* : ∑residuals_i^2 on diagonal for variance of trait i,
#        ∑(residuals_i*residuals_j) on off-diagonal for covariance between trait i and j
# n    : number of observations
# df   : degree of freedom for the prior
# scale: scale parameters for the prior (different from scale in Scale-Inv-chi2)
#
#*******************************************************************************
# single/multiple trait i.i.d/non-i.i.d random effects (others)
# (e.g.,polygenic effects in Pedigree-based BLUP)
#*******************************************************************************
#Given Scale-Inverse-Wishart(ν,SSE)=Inverse-Gamma(ν/2,SSE/2)=Scale-Inv-chi2(ν,SSE/ν),
#1st parameter: df or shape, 2nd parameter: scale,
#ν is (posterior) df, and SSE is the (posterior) sum of square.
#(wikipedia:Inverse-Wishart distribution;Scaled inverse chi-squared distribution)
#
#covariance matrix should be sampled from Scale-Inverse-Wishart(ν,SSE) and scalar
#variance from Scale-Inv-chi2(ν,SSE/ν), for coding simplicity, scalar variance,
#e.g., scaled-Inv-χ2 in single-trait PBLUP without correlated maternal effect,
#is treated as a 1x1 matrix and all are sampled from Scale-Inverse-Wishart(ν,SSE).
#*******************************************************************************

#*******************************************************************************
#Note:
#1. Scale parameters (2nd parameter) is SSE/df in scale-inv-⁠χ2(df,scale)
#   but sum of square (SSE) in scale-Inv-Wishart(df,scale)
#   (from distributions.jl/wikipedia)
#2. Distributions.jl package always returns values of type Float64
#*******************************************************************************
#single-trait i.i.d (residuzl, marker effect variance)
function sample_variance(x, n, df, scale)
    return (dot(x,x) + df*scale)/rand(Chisq(n+df))
end

function sample_variance(x, n, df, scale, invweights)
    sample_variance(x.* (invweights!=false ? sqrt.(invweights) : 1.0),n,df,scale)
end

#multi-trait i.i.d  #?reduce(hcat,array of array)' may be used to replace loops with matrix multiplication
function sample_variance(ycorr_array, nobs, df, scale, invweights, constraint; binary_trait_index=false)
    if invweights != false
        invweights  = Diagonal(invweights)
    end
    ntraits = length(ycorr_array)
    SSE   = zeros(ntraits,ntraits)
    for traiti = 1:ntraits
        ycorri    = ycorr_array[traiti]
        for traitj = traiti:ntraits
            ycorrj    = ycorr_array[traitj]
            SSE[traiti,traitj] = (invweights == false) ? dot(ycorri,ycorrj) : ycorri'*invweights*ycorrj
            if constraint == true #diagonal elements only
                break
            end
            SSE[traitj,traiti] = SSE[traiti,traitj]
        end
    end
    if constraint == false
        if binary_trait_index==false
            R  = rand(InverseWishart(df + nobs, convert(Array,Symmetric(scale + SSE))))
        else
            R  = sample_from_conditional_inverse_Wishart(df + nobs, convert(Array,Symmetric(inv(scale + SSE))), binary_trait_index)
        end
    else  #diagonal elements only, from scale-inv-⁠χ2
        R  = Diagonal(zeros(ntraits))
        for traiti = 1:ntraits
            R[traiti,traiti]= (SSE[traiti,traiti]+df*scale[traiti,traiti])/rand(Chisq(nobs+df)) 
        end
    end
    return R
end

#Futher update is needed
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
       G0 = rand(InverseWishart(random_term.Gi.df + q, convert(Array,Symmetric(random_term.Gi.scale + S))))
       if mme.MCMCinfo.double_precision == false
           G0 = Float32.(G0)
       end
       random_term.GiOld.val = copy(random_term.GiNew.val)
       random_term.GiNew.val = copy(inv(G0))
       random_term.Gi.val    = copy(inv(G0))
       if random_term.randomType == "A"
           mme.Gi    = random_term.Gi
       end
    end
end
########################################################################
# Marker Effects Variance
########################################################################
function sample_marker_effect_variance(Mi)
    if Mi.method == "BayesL"
        invweights = 1 ./ Mi.gammaArray
    elseif Mi.method == "GBLUP"
        invweights = 1 ./ Mi.D
    else
        invweights = false
    end
    if Mi.ntraits == 1
        if Mi.method in ["BayesC","BayesL","RR-BLUP","GBLUP"]
            nloci = Mi.method == "BayesC" ? sum(Mi.δ[1]) : Mi.nMarkers
            Mi.G.val  = sample_variance(Mi.α[1], nloci, Mi.G.df, Mi.G.scale, invweights)
            if Mi.method == "BayesL"
                sampleGammaArray!(Mi.gammaArray,Mi.α,Mi.G.val) # MH sampler of gammaArray (Appendix C in paper)
            end
        elseif Mi.method == "BayesB"
            for j=1:Mi.nMarkers
                Mi.G.val[j] = sample_variance(Mi.β[1][j],1,Mi.G.df, Mi.G.scale)
            end
        end
    else
        if Mi.method in ["RR-BLUP","BayesC","BayesL","GBLUP"]
            data = (Mi.method == "BayesC" ? Mi.β : Mi.α)
            Mi.G.val =sample_variance(data, Mi.nMarkers, Mi.G.df, Mi.G.scale, invweights, Mi.G.constraint)
            if Mi.method == "BayesL"
                sampleGammaArray!(Mi.gammaArray,Mi.α,Mi.G.val) #MH sampler of gammaArray (Appendix C in paper)
            end
        elseif Mi.method == "BayesB" #potential slowdown (scalar multiplication is used instead of matrices)
            marker_effects_matrix = reduce(hcat,Mi.β)'
            for i = 1:Mi.nMarkers
                data    = marker_effects_matrix[:,i]
                Mi.G.val[i] = sample_variance(data, 1, Mi.G.df, Mi.G.scale, false, Mi.G.constraint)
            end
        end
    end
end

function sampleGammaArray!(gammaArray,alphaArray,mmeMG)
    Gi = inv(mmeMG)
    nMarkers = size(gammaArray,1)
    ntraits  = length(alphaArray)

    Q  = zeros(nMarkers)
    ntraits > 1 ? calcMTQ!(Q,nMarkers,ntraits,alphaArray,Gi) : calcSTQ!(Q,nMarkers,alphaArray[1],Gi)
    gammaDist = Gamma(0.5,4) # 4 is the scale parameter, which corresponds to a rate parameter of 1/4
    candidateArray = 1 ./ rand(gammaDist,nMarkers)
    uniformArray = rand(nMarkers)
    acceptProbArray = exp.(Q ./4 .*(2 ./ gammaArray - candidateArray))
    replace = uniformArray .< acceptProbArray
    gammaArray[replace] = 2 ./ candidateArray[replace]
end

function calcMTQ!(Q,nMarkers,ntraits,alphaArray,Gi)
    for locus = 1:nMarkers
      for traiti = 1:ntraits
          for traitj = 1:ntraits
              Q[locus] += alphaArray[traiti][locus]*alphaArray[traitj][locus]*Gi[traiti,traitj]
          end
      end
    end
end

function calcSTQ!(Q,nMarkers,alphaArray,Gi)
    Q .= alphaArray.^2 ./Gi
end
