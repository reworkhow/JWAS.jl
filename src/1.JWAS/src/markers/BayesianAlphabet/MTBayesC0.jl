function sampleMarkerEffects!(xArray,xpx,wArray,alphaArray,meanAlpha,invR0,invG0,
                               iIter)
    # function to sample effects
    # invR0 is inverse of R0
    # invG0 is inverse of G0

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
        Lhs = dot(x,x)*invR0 + invG0
        #println("Rhs",Rhs)
        #println("Lhs",Lhs)
        for trait = 1:nTraits
            lhs  = Lhs[trait,trait]
            ilhs = 1/lhs
            rhs  = (Rhs[trait] - Lhs[trait,:]*β)[1,1]
            mu   = ilhs*rhs + β[trait]
            alphaArray[trait][j] = mu + randn()*sqrt(ilhs)
            β[trait] = alphaArray[trait][j]
        end
        # adjust for locus j
        for trait = 1:nTraits
            #wArray[trait][:] = wArray[trait][:] - x*alphaArray[trait][j]
            BLAS.axpy!(oldβ[trait]-β[trait],x,wArray[trait])
            meanAlpha[trait][j] += (β[trait] - meanAlpha[trait][j])/iIter
        end

    end
end
