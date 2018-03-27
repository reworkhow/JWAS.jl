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

###############################################################################
#	Heterogeneous residuals, "weight variables"
# an argument "Rinv" is required
###############################################################################



###############################################################################
#1, comput Xj'Xj
###############################################################################
#OLD
##xpx       = [dot(X[:,i],X[:,i]) for i=1:size(X,2)]

#NEW
#function getXpRinvX(X, Rinv)
#    ncol = size(X)[2]
#    XpRinvX = [((X[:,i].*Rinv)'X[:,i])[1]::Float64 for i=1:ncol]
#    return XpRinvX
#end
#XpRinvX = getXpRinvX(X, Rinv)

###############################################################################
#2 sample marker effcts
###############################################################################
#OLD
##above

#NEW
#
# function sampleEffectsBayesC!(nMarkers,
#                               xArray,
#                               XpRinvX,
#                               yCorr,
#                               α,
#                               δ,
#                               vare,
#                               varEffects,
#                               π,
#                               Rinv)
#
#     logPi         = log(π)
#     logPiComp     = log(1.0-π)
#     logVarEffects = log(varEffects)
#     logDelta0     = logPi
#     invVarRes     = 1.0/vare
#     invVarEffects = 1.0/varEffects
#     nLoci = 0
#
#     for j=1:nMarkers
#         x = xArray[j]
#         rhs = (dot(x.*Rinv,yCorr) + XpRinvX[j]*α[j])*invVarRes
#         lhs = XpRinvX[j]*invVarRes + invVarEffects
#         invLhs = 1.0/lhs
#         gHat   = rhs*invLhs
#         logDelta1  = -0.5*(log(lhs) + logVarEffects - gHat*rhs) + logPiComp
#         probDelta1 = 1.0/(1.0 + exp(logDelta0 - logDelta1))
#         oldAlpha = α[j]
#
#         if(rand()<probDelta1)
#             δ[j] = 1
#             α[j] = gHat + randn()*sqrt(invLhs)
#             BLAS.axpy!(oldAlpha-α[j],x,yCorr)
#             nLoci = nLoci + 1
#         else
#             if (oldAlpha!=0)
#                 BLAS.axpy!(oldAlpha,x,yCorr)
#             end
#             δ[j] = 0
#             α[j] = 0
#         end
#     end
#     return nLoci
# end

###############################################################################
#3, sample residual variances
###############################################################################

#OLD
#vare = sampleVariance(ycorr,nInd,nuRes,scaleRes)

#NEW
#vare = sampleVariance(yCorr.*RinvSqrt, nObs, nuRes, scaleRes)
