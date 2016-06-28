
function sampleMarkerEffectsBayesC!(xArray,xpx,wArray,alphaArray,meanAlphaArray,
                                    deltaArray,meanDeltaArray,
                                    uArray,meanuArray,
                                    invR0,invG0,iIter,logPi)

    nMarkers = length(xArray)
    nTraits = length(alphaArray)
    Ginv     = invG0
    Rinv     = invR0

       α = zeros(nTraits)
    newu = zeros(nTraits)
       δ = zeros(nTraits)

       w = zeros(nTraits) #for rhs

    for j=1:nMarkers

        x    = xArray[j]

        for k=1:nTraits
               α[k] = alphaArray[k][j]
            newu[k] = uArray[k][j]
               δ[k] = deltaArray[k][j]
        end

        oldu = copy(newu)

        #α    = alpha[j,:]
        #δ    = delta[j,:]
        #newu = u[j,:]
        #oldu = u[j,:]

        for trait = 1:nTraits
            w[trait] = dot(x,wArray[trait])+xpx[j]*oldu[trait]
            #w    = x'yCor+xpx[j]*oldu #adjust y except locus j
        end

        for k=1:nTraits
            Ginv11 = Ginv[k,k]
            nok    = deleteat!(collect(1:nTraits),k)
            Ginv12 = Ginv[k,nok]
            C11    = Ginv11+Rinv[k,k]*xpx[j]
            #C12    = Ginv12+mpm[j]*Rinv[k,nok]*diagm(δ[nok])#work for Rinv_sub
            C12    = Ginv12+xpx[j]*Rinv[k,nok].*δ[nok]' #δ[:,nok] : row vector

            invLhs0  = 1/Ginv11
            rhs0     = - Ginv12*α[nok]
            gHat0    = (rhs0*invLhs0)[1,1]

            invLhs1  = 1/C11
            rhs1     = w'*Rinv[:,k]-C12*α[nok] #w transpose
            gHat1    = (rhs1*invLhs1)[1,1]

            logDelta0  = -0.5*(log(Ginv11)- gHat0^2*Ginv11) + logPi[1,k] #logPi
            logDelta1  = -0.5*(log(C11)-gHat1^2*C11) + logPi[2,k] #logPiComp
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
        #yCor[:,:] = yCor + m*(oldu-newu)

        # adjust for locus j
        for trait = 1:nTraits
            #wArray[trait][:] = wArray[trait][:] - x*alphaArray[trait][j]
            BLAS.axpy!(oldu[trait]-newu[trait],x,wArray[trait])
            meanAlphaArray[trait][j] += (α[trait] - meanAlphaArray[trait][j])/iIter
            meanDeltaArray[trait][j] += (δ[trait] - meanDeltaArray[trait][j])/iIter
            meanuArray[trait][j]     += (newu[trait] - meanuArray[trait][j])/iIter

            alphaArray[trait][j]      = α[trait]
            deltaArray[trait][j]      = δ[trait]
            uArray[trait][j]          = newu[trait]
        end
    end
end
