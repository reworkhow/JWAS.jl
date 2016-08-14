function sampleEffectsBayesC!(mats::GibbsMats,current::Current,out::Output)

    logPi         = log(current.π)
    logPiComp     = log(1.0-current.π)
    logDelta0     = logPi
    invVarRes     = 1.0/current.varResidual
    logVarEffect  = log(current.varEffect)
    invVarEffect  = 1.0/current.varEffect

    xArray        = mats.xArray
    xpx           = mats.xpx
    yCorr         = current.yCorr
    α             = current.α
    δ             = current.δ
    u             = current.u
    varRes        = current.varResidual
    λ             = current.varResidual/current.varEffect
    iIter         = 1/current.iter
    nLoci         = 0

    nEffects      = length(u)
    for j=1:nEffects
        x = xArray[j]
        rhs = (dot(x,yCorr) + xpx[j]*α[j])*invVarRes
        lhs = xpx[j]*invVarRes + invVarEffect
        invLhs = 1.0/lhs
        gHat   = rhs*invLhs
        logDelta1  = -0.5*(log(lhs) + logVarEffect - gHat*rhs) + logPiComp
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

    current.nLoci          =  nLoci
    out.meanMarkerEffects +=  (α - out.meanMarkerEffects)*iIter #the way in sampe_ycor may be better
    out.modelFreq         +=  (δ - out.modelFreq)*iIter
end

export sampleEffectsBayesC!
