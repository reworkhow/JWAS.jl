function sampleEffectsBayesC!(xArray,
                              xpx,
                              yCorr,
                              α,
                              δ,
                              vare,
                              varEffects,
                              π)

    logPi         = log(π)
    logPiComp     = log(1.0-π)
    logVarEffects = log(varEffects)
    logDelta0     = logPi
    invVarRes     = 1.0/vare
    invVarEffects = 1.0/varEffects
    nLoci         = 0

    nMarkers      = length(α)

    for j=1:nMarkers
        x = xArray[j]
        rhs = (dot(x,yCorr) + xpx[j]*α[j])*invVarRes
        lhs = xpx[j]*invVarRes + invVarEffects
        invLhs = 1.0/lhs
        gHat   = rhs*invLhs
        logDelta1  = -0.5*(log(lhs) + logVarEffects - gHat*rhs) + logPiComp
        probDelta1 = 1.0/(1.0 + exp(logDelta0 - logDelta1))
        oldAlpha = α[j]

        if(rand()<probDelta1)
            δ[j] = 1
            α[j] = gHat + randn()*sqrt(invLhs)
            BLAS.axpy!(oldAlpha-α[j],x,yCorr)
            nLoci = nLoci + 1
        else
            if (oldAlpha!=0)
                BLAS.axpy!(oldAlpha,x,yCorr)
            end
            δ[j] = 0
            α[j] = 0
        end
    end
    return nLoci
end
