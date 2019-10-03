#Created by Rohan Fernando and Hao Cheng on 2/22/2015

function sampleEffectsBayesC0!(nMarkers,
                        nObs,
                        xArray::Array{Array{Float32,1},1},
                        xpx::Array{Float32,1},
                        yCorr::Array{Float32,1},
                        α::Array{Float32,1},
                        meanAlpha::Array{Float32,2},
                        vare::Float32,
                        varEffects::Float32,
                        x::Array{Float32,1})
    for j=1:nMarkers
        x = xArray[j]
        rhs::Float32 = dot(x,yCorr) + xpx[j]*α[j]
        lhs::Float32      = xpx[j] + vare/varEffects
        invLhs::Float32   = 1.0/lhs
        mean::Float32     = invLhs*rhs
        oldAlpha::Float32 = α[j]
        α[j] = mean + randn()*sqrt(invLhs*vare)
        BLAS.axpy!(oldAlpha-α[j],x,yCorr)
    end
    nothing
end

function BayesC0!(options,X,y)
    seed            =   options.seed    # set the seed for the random number generator
    chainLength     =   options.chainLength   # number of iterations
    probFixed       =   options.probFixed    # parameter "pi" the probability SNP effect is zero
    estimatePi      =   options.estimatePi   # "yes" or "no"
    dfEffectVar     =   options.dfEffectVar    # hyper parameter (degrees of freedom) for locus effect variance
    nuRes           =   options.nuRes    # hyper parameter (degrees of freedom) for residual variance
    varGenotypic    =   options.varGenotypic    # used to derive hyper parameter (scale) for locus effect variance
    varResidual     =   options.varResidual    # used to derive hyper parameter (scale) for locus effect variance
    scaleVar        =   varGenotypic*(dfEffectVar-2)/dfEffectVar    # scale factor for locus effects
    scaleRes        =   varResidual*(nuRes-2)/nuRes                 # scale factor for residual varianc
    numIter         =   chainLength

    nObs,nMarkers=size(X)

    xpx = [(X[:,i]'X[:,i])[1]::Float32 for i=1:nMarkers]
    x   = Array(Float32,nObs)
    xArray = Array(Array{Float32,1},nMarkers)
    for i=1:nMarkers
        xArray[i] = get_column(X,i)
    end

    #initial values
    logPi      =log(probFixed)
    logPiComp  =log(1-probFixed)
    vare       = float32(1.0)
    varEffects = float32(1.0)
    yCorr      = copy(y)
    α          = float32(fill(0.0,nMarkers))
    #return values
    meanMu     = 0
    meanAlpha  = float32(zeros(nMarkers,1))

    chi1=Chisq(nObs+nuRes)
    chi2=Chisq(dfEffectVar+nMarkers)

    for i=1:numIter
        # sample intercept
        yCorr  = float32(yCorr+mu)
        rhs    = float32(sum(yCorr))
        invLhs = float32(1.0/(nObs))
        mean   = rhs*invLhs
        mu     = mean + randn()*sqrt(invLhs*vare)
        yCorr  = yCorr - mu
        meanMu = meanMu + (mu - meanMu)/i

        # sample effects
        sampleEffectsBayesC0!(nMarkers,nObs,xArray,xpx,yCorr,α,meanAlpha,vare,varEffects,x)
        meanAlpha = meanAlpha + (α - meanAlpha)/i

        # sample residula variance
        vare = float32((dot(yCorr,yCorr)+nuRes*scaleRes)/rand(chi1))

        #sameple locus effect variance
        varEffects = float32((scaleVar*dfEffectVar + dot(α,α))/rand(chi2) )

        if (i%100)==0
            yhat = meanMu+X*meanAlpha
            resCorr = cor(yhat,yhat) #modify
            println ("Correlation of between true and predicted breeding value: ", resCorr[1])
        end
    end
end
