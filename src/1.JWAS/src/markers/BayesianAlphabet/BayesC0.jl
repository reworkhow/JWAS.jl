function sampleEffectsBayesL!(xArray,xpx,yCorr,α,gammaArray,vRes,vEff,Rinv) # sample vare and vara
    nMarkers      = length(α)
    λ    = vRes/vEff

   function getLambdaFunction(x)
        f1(x,j,λ)  = λ/x[j]
        f2(x,j,λ)  = λ
        length(x)>1 ? f1 : f2
    end
    getLambda = getLambdaFunction(gammaArray)
    for j=1:nMarkers
        x        = xArray[j]
        rhs      = dot(x.*Rinv,yCorr) + xpx[j]*α[j,1]
        lhs      = xpx[j] + getLambda(gammaArray,j,λ)
        invLhs   = 1.0/lhs
        mean     = invLhs*rhs
        oldAlpha = α[j,1]
        α[j]     = mean + randn()*sqrt(invLhs*vRes)
        BLAS.axpy!(oldAlpha-α[j,1],x,yCorr)
    end
end

function sampleEffectsBayesC0!(xArray,xpx,yCorr,α,vRes,vEff)
    sampleEffectsBayesL!(xArray,xpx,yCorr,α,[1.0],vRes,vEff)
end
