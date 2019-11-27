function BayesL!(xArray,xpx,xRinvArray,xpRinvx, #Heterogeneous residuals
                 yCorr,α,gammaArray,vRes,vEff)
    nMarkers      = length(α)
    λ    = vRes/vEff

   function getLambdaFunction(x)
        f1(x,j,λ)  = λ/x[j]
        f2(x,j,λ)  = λ
        length(x)>1 ? f1 : f2
    end
    getLambda = getLambdaFunction(gammaArray)
    for j=1:nMarkers
        x, xRinv = xArray[j], xRinvArray[j]
        rhs      = dot(xRinv,yCorr) + xpx[j]*α[j]
        lhs      = xpx[j] + getLambda(gammaArray,j,λ)
        invLhs   = 1.0/lhs
        mean     = invLhs*rhs
        oldAlpha = α[j]
        α[j]     = mean + randn()*sqrt(invLhs*vRes)
        BLAS.axpy!(oldAlpha-α[j],x,yCorr)
    end
end

function BayesC0!(xArray,xpx,xRinvArray,xpRinvx,yCorr,α,vRes,vEff)
    BayesL!(xArray,xpx,xRinvArray,xpRinvx,yCorr,α,[1.0],vRes,vEff)
end
