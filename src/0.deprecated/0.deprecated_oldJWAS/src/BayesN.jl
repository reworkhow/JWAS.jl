#Created by Jian Zeng on 3/6/2015

function make_windows(map, width)
    nMarkers = size(map)[1]
    windows  = Array(Array{Float64,1},nMarkers)
    chrStart = map[1,2]
    posStart = map[1,3]
    winId = 1
    for i=1:nMarkers
        chr = map[i,2]
        pos = map[i,3]
        if (chr > chrStart | pos-posStart > width)
            winId = winId + 1
        end
        windows[winId] = [windows[winId],i]
    end
    resize!(windows, winId)
end


function sampleEffectsBayesN!(yCorr, nObs, nWindows, windows, windowSizes, xArray, XpRinvX, Rinv,
                              u, α, Δ, δ, Π, varWindows, vare, k)
    logPi         = log(Π)
    logPiComp     = log(1.0-Π)
    invVarRes     = 1.0/vare
    nWins = 0
    nLoci = 0
    for i=1:nWindows
        window = windows[i]
        windowEffect = xArray[window]*α[window]
        diffQuadratic = 2*dot(yCorr.*Rinv,windowEffect)
        if (Δ[i]==1)
            diffQuadratic = diffQuadratic + dot(windowEffect.*Rinv,windowEffect)
        else
            diffQuadratic = diffQuadratic - dot(windowEffect.*Rinv,windowEffect)
        end
        logDelta0ToDelta1 = -0.f*diffQuadratic*invVarRes + logPi - logPiComp
        probDelta1 = 1.0/(1.0 + exp(logDelta0ToDelta1))
        if (rand() < probDelta1)
            if (Δ[i]==0)
                BLAS.axpy!(-1,windowEffect,yCorr)
            end
            Δ[i] = 1
            π = (windowSizes[i] - k)/windowSizes[i]
            if (k > windowSizes[i])
                π = 0
            end
            nMrks = sampleEffectsBayesCPi(windowSizes[i], xArray[window], XpRinvX[window], yCorr,
                                          α[window], δ[window], vare, varWindows[i], π, Rinv)
            u[window] = α[window]
            nLoci = nLoci + nMrks
            nWins = nWins + 1
        else
            Δ[i] = 0
            sdWindow = sqrt(varWindows[i])
            for j=1:windowSizes[i]
                if (rand() < 0.5)
                    δ[window[j]] = 1
                    α[window[j]] = randn()*sdWindow
                else
                    δ[window[j]] = 0
                    α[window[j]] = 0
                end
            end
            u[window] = 0
        end
    end

    return [nWins nLoci]
end


function BayesN!(options,X,y,C,Rinv)
    # input options
    seed            =   options.seed            # set the seed for the random number generator
    chainLength     =   options.chainLength     # number of iterations
    probFixed       =   options.probFixed       # parameter "pi" the probability window effect is zero
    estimatePi      =   options.estimatePi      # "yes" or "no"
    dfEffectVar     =   options.dfEffectVar     # hyper parameter (degrees of freedom) for locus effect variance
    nuRes           =   options.nuRes           # hyper parameter (degrees of freedom) for residual variance
    varGenotypic    =   options.varGenotypic    # used to derive hyper parameter (scale) for locus effect variance
    varResidual     =   options.varResidual     # used to derive hyper parameter (scale) for locus effect variance
    scaleRes        =   options.scaleRes        # scale factor for residual varianc
    k               =   options.fittedSNPperWindow
    windowWidth     =   options.windowWidth     # in mega-basepaire unit
    map             =   options.markerMap       # physical map positions of SNPs

    # prepare
    numIter         =   chainLength
    nObs,nMarkers   =   size(X)
    nFixedEffects   =   size(C)[2]

    windows = make_windows(map, windowWidth)
    nWindows = length(windows)
    windowSizes = [length(windows[i]) for i=1:nWindows]
    markerMeans = center!(X)
    XpRinvX = getXpRinvX(X, Rinv)

    β          =  zeros(nFixedEffects)  # sample of fixed effects
    α          =  zeros(nMarkers)       # sample of partial marker effects conditional only on δ
    u          =  zeros(nMarkers)       # sample of marker effects
    Δ          =  zeros(nWindows)       # inclusion indicator for windows
    δ          =  zeros(nMarkers)       # inclusion indicator for marker effects
    p          =  markerMeans/2.0
    mean2pq    =  (2*p*(1-p)')[1,1]
    varEffects =  varGenotypic/((1-probFixed)*mean2pq)
    scaleVar   =  varEffects*(dfEffectVar-2)/dfEffectVar
    vare       =  varResidual
    Π          =  probFixed
    yCorr      =  copy(y)
    RinvSqrt   =  sqrt(Rinv)
    varWindows =  fill(varEffects, nWindows)

    # output variables
    meanFxdEff   = zeros(nFixedEffects)
    meanMrkEff   = zeros(nMarkers)
    mdlFrqMrk    = zeros(nMarkers)
    mdlFrqWin    = zeros(nMarkers)
    Pi           = zeros(chainLength)
    genVar       = zeros(chainLength)
    resVar       = zeros(chainLength)

    # analysis start
    for i=1:numIter
        # sample fixed effects
        sampleFixedEffects!(yCorr, nFixedEffects, C, Rinv, β, vare)
        meanFxdEff = meanFxdEff + (β - meanFxdEff)/i

        # sample marker effects
        nWins, nLoci = sampleEffectsBayesN!(yCorr, nObs, nWindows, windows, windowSizes, xArray, XpRinvX, Rinv,
                                            u, α, Δ, δ, Π, varWindows, vare, k)
        meanMrkEff   = meanMrkEff + (u - meanMrkEff)/i
        mdlFrqMrk    = mdlFrqMrk  + (δ - mdlFrqMrk )/i
        mdlFrqWin    = mdlFrqWin  + (Δ - mdlFrqWin )/i
        genVar[i]    = var(X*α)

        # sample residula variance
        vare = sampleVariance(yCorr.*RinvSqrt, nObs, nuRes, scaleRes)
        resVar[i] = vare

        # sameple effect variance for every window
        for j=1:nWindows
            varWindows[j] = sampleVariance(α[windows[j]], windowSizes[j], dfEffectVar, scaleVar)
        end

        # sample π
        if (estimatePi == "yes")
            Π = samplePi(nWins, nWindows)[1,1]
        end
        Pi[i] = Π[1]

        # display progress
        if (i%100)==0
            yCorr = y - C*β - X*u  # remove rounding errors
            println ("Iter ", i, ", number of windows ", nWins, ", number of loci ", nLoci)
        end
    end

    # output list
    output = Dict()
    output["posterior mean of fixed effects"]         = meanFxdEff
    output["posterior mean of marker effects"]        = meanMrkEff
    output["window model frequency"]                  = mdlFrqWin
    output["marker model frequency"]                  = mdlFrqMrk
    output["posterior sample of Pi"]                  = Pi
    output["posterior sample of genotypic variance"]  = genVar
    output["posterior sample of residual variance"]   = resVar

    return output
end

