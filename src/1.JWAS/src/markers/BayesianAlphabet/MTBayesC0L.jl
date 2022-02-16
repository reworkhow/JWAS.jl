function MTBayesL!(genotypes,ycorr_array,vare)
    MTBayesL!(genotypes.mArray,genotypes.mRinvArray,genotypes.mpRinvm,
            ycorr_array,genotypes.α,genotypes.gammaArray,vare,genotypes.G)
end

function MTBayesC0!(genotypes,ycorr_array,vare)
    MTBayesL!(genotypes.mArray,genotypes.mRinvArray,genotypes.mpRinvm,
            ycorr_array,genotypes.α,[1.0],vare,genotypes.G)
end

function MTBayesL!(xArray,xRinvArray,xpRinvx,
                   wArray,alphaArray,gammaArray,
                   vare,varEffect)
    invR0    = inv(vare)
    invG0    = inv(varEffect)
    nMarkers = length(xArray)
    ntraits  = length(alphaArray)
    Lhs      = zeros(ntraits,ntraits)
    Rhs      = zeros(ntraits)
    newα     = zeros(typeof(alphaArray[1][1]),ntraits)
    oldα     = zeros(typeof(alphaArray[1][1]),ntraits)


    function getGi_function(gammaArray)
        f1(x,j,Gi)  = Gi ./ x[j]
        f2(x,j,Gi)  = Gi
        size(gammaArray,1)>1 ? f1 : f2
    end
    getGi = getGi_function(gammaArray)

    for marker=1:nMarkers
        x, xRinv = xArray[marker], xRinvArray[marker]
        for trait=1:ntraits #unadjust for locus j
            oldα[trait] = newα[trait] = alphaArray[trait][marker]
            Rhs[trait]  = dot(xRinv,wArray[trait])+xpRinvx[marker]*oldα[trait]
        end
        Rhs = invR0*Rhs
        Lhs = xpRinvx[marker]*invR0 + getGi(gammaArray,marker,invG0)
        for trait = 1:ntraits
            lhs  = Lhs[trait,trait]
            ilhs = 1/lhs
            rhs  = Rhs[trait] - (Lhs[trait,:]'newα)[1]
            mu   = ilhs*rhs + newα[trait]
            alphaArray[trait][marker] = mu + randn()*sqrt(ilhs)
            newα[trait] = alphaArray[trait][marker]
        end
        for trait = 1:ntraits
            BLAS.axpy!(oldα[trait]-newα[trait],x,wArray[trait])
        end
    end
end
