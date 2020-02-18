function BayesL!(xArray,xRinvArray,xpRinvx,
                 yCorr,
                 α,gammaArray,
                 vRes,vEff)
    nMarkers = length(α)
    λ        = vRes/vEff
    function get_lambda_function(x)
        f1(x,j,λ)  = λ/x[j]
        f2(x,j,λ)  = λ
        length(x)>1 ? f1 : f2
    end
    getlambda = get_lambda_function(gammaArray)
    for j=1:nMarkers
        x, xRinv = xArray[j], xRinvArray[j]
        rhs      = dot(xRinv,yCorr) + xpRinvx[j]*α[j]
        lhs      = xpRinvx[j] + getlambda(gammaArray,j,λ)
        invLhs   = 1.0/lhs
        mean     = invLhs*rhs
        oldAlpha = α[j]
        α[j]     = mean + randn()*sqrt(invLhs*vRes)
        BLAS.axpy!(oldAlpha-α[j],x,yCorr)
    end
end

function BayesC0!(xArray,xRinvArray,xpRinvx,yCorr,α,vRes,vEff)
    BayesL!(xArray,xRinvArray,xpRinvx,yCorr,α,[1.0],vRes,vEff)
end
