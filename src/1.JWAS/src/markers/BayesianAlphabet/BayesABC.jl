function BayesABC!(xArray,xRinvArray,xpRinvx,
                   yCorr,
                   α,β,δ,
                   vare,varEffects,π)

    logPi         = log(π)
    logPiComp     = log(1-π)
    logDelta0     = logPi
    invVarRes     = 1/vare
    invVarEffects = 1 ./  varEffects
    logVarEffects = log.(varEffects)
    nLoci         = 0
    nMarkers      = length(α)

    for j=1:nMarkers
        x, xRinv = xArray[j], xRinvArray[j]
        rhs = (dot(xRinv,yCorr) + xpRinvx[j]*α[j])*invVarRes
        lhs = xpRinvx[j]*invVarRes + invVarEffects[j]
        invLhs = 1/lhs
        gHat   = rhs*invLhs
        logDelta1  = -0.5*(log(lhs) + logVarEffects[j] - gHat*rhs) + logPiComp
        probDelta1 = 1/(1+ exp(logDelta0 - logDelta1))
        oldAlpha = α[j]

        if(rand()<probDelta1)
            δ[j] = 1
            β[j] = gHat + randn()*sqrt(invLhs)
            α[j] = β[j]
            BLAS.axpy!(oldAlpha-α[j],x,yCorr)
            nLoci = nLoci + 1
        else
            if (oldAlpha!=0)
                BLAS.axpy!(oldAlpha,x,yCorr)
            end
            δ[j] = 0
            β[j] = randn()*sqrt(varEffects[j])
            α[j] = 0
        end
    end
    return nLoci
end
