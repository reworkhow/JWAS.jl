function sampleMarkerEffectsMTBayesL!(xArray,xpx,wArray,alphaArray,meanAlpha,gammaArray,invR0,invG0,
                               iIter,burnin)
    # function to sample effects
    # invR0 is inverse of R0
    # invG0 is inverse of G0

    getGi(x,j,Gi) = size(gammaArray,1)>1 ? Gi ./ x[j] : Gi
    nEffects = length(xArray)
    nTraits  = length(alphaArray)
    Lhs = zeros(nTraits,nTraits)
    Rhs = zeros(nTraits)
    β = zeros(nTraits)             # vector of effects for locus j across all traits
    for j=1:nEffects
        for k=1:nTraits
            β[k] = alphaArray[k][j]
        end
        oldβ = copy(β) #add for optimization
        x = xArray[j]
        # unadjust for locus j
        for trait = 1:nTraits
            #wArray[trait][:] = wArray[trait][:] + x*alphaArray[trait][j]
            #Rhs[trait] = dot(x,wArray[trait])
            Rhs[trait] = dot(x,wArray[trait])+xpx[j]*alphaArray[trait][j]
        end
        #println(wArray[1])
        Rhs = invR0*Rhs
        Lhs = dot(x,x)*invR0 + getGi(gammaArray,j,invG0)
        #println("Rhs",Rhs)
        #println("Lhs",Lhs)
        for trait = 1:nTraits
            lhs  = Lhs[trait,trait]
            ilhs = 1/lhs
            rhs  = Rhs[trait] - (Lhs[trait,:]'β)[1]
            mu   = ilhs*rhs + β[trait]
            alphaArray[trait][j] = mu + randn()*sqrt(ilhs)
            β[trait] = alphaArray[trait][j]
        end
        # adjust for locus j
        for trait = 1:nTraits
            #wArray[trait][:] = wArray[trait][:] - x*alphaArray[trait][j]
            BLAS.axpy!(oldβ[trait]-β[trait],x,wArray[trait])
            if iIter>burnin
                meanAlpha[trait][j] += (β[trait] - meanAlpha[trait][j])/(iIter-burnin)
            end
        end

    end
end

sampleMarkerEffectsMTBayesC0!(xArray,xpx,wArray,alphaArray,meanAlpha,invR0,invG0,
                               iIter,burnin) =
  sampleMarkerEffectsMTBayesL!(xArray,xpx,wArray,alphaArray,meanAlpha,[1.0],invR0,invG0,
                                                              iIter,burnin)
