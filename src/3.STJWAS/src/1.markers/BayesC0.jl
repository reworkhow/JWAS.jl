function sampleEffectsBayesC0!(xArray,xpx,yCorr,α,vRes,vEff)#sample vare and vara
    nMarkers      = length(α)
    λ    = vRes/vEff
    for j=1:nMarkers
        x        = xArray[j]
        rhs      = dot(x,yCorr) + xpx[j]*α[j,1]
        lhs      = xpx[j] + λ
        invLhs   = 1.0/lhs
        mean     = invLhs*rhs
        oldAlpha = α[j,1]
        α[j]     = mean + randn()*sqrt(invLhs*vRes)
        BLAS.axpy!(oldAlpha-α[j,1],x,yCorr)
    end
end
