function megaBayesL!(genotypes,wArray,vare)
    Threads.@threads for i in 1:length(wArray) #ntraits
        BayesL!(genotypes.mArray,genotypes.mRinvArray,genotypes.mpRinvm,
            wArray[i],genotypes.α[i],genotypes.gammaArray,vare[i,i],genotypes.G.val[i,i])
    end
end

function BayesL!(genotypes,ycorr,vare)
    BayesL!(genotypes.mArray,genotypes.mRinvArray,genotypes.mpRinvm,
            ycorr,genotypes.α[1],genotypes.gammaArray,vare,genotypes.G.val)
end

function megaBayesC0!(genotypes,wArray,vare)
    Threads.@threads for i in 1:length(wArray) #ntraits
        BayesL!(genotypes.mArray,genotypes.mRinvArray,genotypes.mpRinvm,
                wArray[i],genotypes.α[i],[1.0],vare[i,i],genotypes.G.val[i,i])
    end
end

function BayesC0!(genotypes,ycorr,vare)
    BayesL!(genotypes.mArray,genotypes.mRinvArray,genotypes.mpRinvm,
            ycorr,genotypes.α[1],[1.0],vare,genotypes.G.val)
end

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
