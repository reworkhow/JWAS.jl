function sampleEffectsBayesB!(xArray,
                              xpx,
                              yCorr,
                              u,
                              α,
                              δ,
                              vare,
                              locusEffectVar,
                              π)

    logPi             = log(π)
    logPiComp         = log(1.0-π)
    logDelta0         = logPi
    invVarRes         = 1.0/vare
    invlocusEffectVar = 1.0./locusEffectVar
    loglocusEffectVar = log.(locusEffectVar)
    nLoci             = 0

    nMarkers      = length(α)

    for j=1:nMarkers
        x = xArray[j]
        rhs = (dot(x,yCorr) + xpx[j]*u[j])*invVarRes
        lhs = xpx[j]*invVarRes + invlocusEffectVar[j]
        invLhs = 1.0/lhs
        gHat   = rhs*invLhs
        logDelta1  = -0.5*(log(lhs) + loglocusEffectVar[j] - gHat*rhs) + logPiComp
        probDelta1 = 1.0/(1.0 + exp(logDelta0 - logDelta1))
        oldu = u[j]

        if(rand()<probDelta1)
            δ[j] = 1
            α[j] = gHat + randn()*sqrt(invLhs)
            u[j] = α[j]
            BLAS.axpy!(oldu-u[j],x,yCorr)
            nLoci = nLoci + 1
        else
            if (oldu!=0)
                BLAS.axpy!(oldu,x,yCorr)
            end
            δ[j] = 0
            α[j] = randn()*sqrt(locusEffectVar[j])
            u[j] = 0
        end
    end
    return nLoci
end
