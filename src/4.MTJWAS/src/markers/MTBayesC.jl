
function sampleMarkerEffectsBayesC!(xArray,xpx,wArray,alphaArray,meanAlphaArray,
                                    deltaArray,meanDeltaArray,
                                    uArray,meanuArray,
                                    invR0,invG0,iIter,BigPi)
    nMarkers = length(xArray)
    nTraits  = length(alphaArray)
    Ginv     = invG0
    Rinv     = invR0

    α        = zeros(nTraits)
    newu     = zeros(nTraits)
    oldu     = zeros(nTraits)
    δ        = zeros(nTraits)
    w        = zeros(nTraits) #for rhs

    for marker=1:nMarkers

        x    = xArray[marker]

        for trait = 1:nTraits
            α[trait]  = alphaArray[trait][marker]
         oldu[trait]  = newu[trait] = uArray[trait][marker]
            δ[trait]  = deltaArray[trait][marker]
            w[trait]  = dot(x,wArray[trait])+xpx[marker]*oldu[trait]
        end

        for k=1:nTraits
            Ginv11 = Ginv[k,k]
            nok    = deleteat!(collect(1:nTraits),k)
            Ginv12 = Ginv[k,nok]
            C11    = Ginv11+Rinv[k,k]*xpx[marker]
            C12    = Ginv12+xpx[marker]*Rinv[k,nok]*diagm(δ[nok])
            #C12    = Ginv12+xpx[marker]*Rinv[k,nok].*δ[nok]' #δ[:,nok] : row vector,

            invLhs0  = 1/Ginv11
            rhs0     = - Ginv12*α[nok]
            gHat0    = (rhs0*invLhs0)[1,1]

            invLhs1  = 1/C11
            rhs1     = w'*Rinv[:,k]-C12*α[nok] #w transpose
            gHat1    = (rhs1*invLhs1)[1,1]

            d0 = copy(δ)
            d1 = copy(δ)
            d0[k] = 0.0
            d1[k] = 1.0

            logDelta0  = -0.5*(log(Ginv11)- gHat0^2*Ginv11) + log(BigPi[d0]) #logPi
            logDelta1  = -0.5*(log(C11)-gHat1^2*C11) + log(BigPi[d1]) #logPiComp

            probDelta1 =  1.0/(1.0+exp(logDelta0-logDelta1))
            if(rand()<probDelta1)
                δ[k] = 1
                α[k] = newu[k] = gHat1 + randn()*sqrt(invLhs1)
            else
                α[k] = gHat0 + randn()*sqrt(invLhs0)
                δ[k] = 0
                newu[k] = 0
            end
        end

        # adjust for locus j
        for trait = 1:nTraits
            BLAS.axpy!(oldu[trait]-newu[trait],x,wArray[trait])
            meanAlphaArray[trait][marker] += (α[trait] - meanAlphaArray[trait][marker])/iIter
            meanDeltaArray[trait][marker] += (δ[trait] - meanDeltaArray[trait][marker])/iIter
            meanuArray[trait][marker]     += (newu[trait] - meanuArray[trait][marker])/iIter

            alphaArray[trait][marker]      = α[trait]
            deltaArray[trait][marker]      = δ[trait]
            uArray[trait][marker]          = newu[trait]
        end
    end
end
