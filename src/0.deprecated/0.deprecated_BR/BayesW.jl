function window_genotype(geno::Genotypes,input::QTL.InputParameters)#assume one WBegin now
    XAll= geno.genotypes
    sel = fill(false,size(XAll,2))
    sel[WBegin:WEnd] = true
    X  = XAll[:,sel]
    X2 = XAll[:,!sel]

    geno4window = make_genotypes(X) #one window genotypes

    s    = sqrt(var(X2,1))
    sel  = s.!=0
    X2   = X2[:,sel]
    G    = X2*X2'/size(X2,2)
    v,U  = eig(G)
    vCum = cumsum(v)./sum(v)
    startAt = 1
    for i=1:length(vCum)
        if (vCum[i] > 0.0045)
            startAt = i
            break
        end
    end
    W = U[:,startAt:end]
    d = v[startAt:end]
    return W,d, geno4window
end


function samplePolygenicEffects!(yCorr,W,d,z,vare,vara)
    λ = vare/vara
    yCorr[:] += W*z
    rhs  = W'yCorr
    di   = λ./d
    zHat = rhs./(1 + di)
    z[:] = zHat + randn(length(zHat)).*sqrt(vare./(1+di))
    yCorr[:] -= W*z
end

function BayesW!(y::Array{Float64,1},
                 geno::Genotypes,
                 fixed::FixedMatrix,
                 input::QTL.InputParameters;outFreq=5000,WBegin,WEnd)

    ###INPUT
    windowSize      =   options.windowSize      # number of markers in window for computing window variance
    numIter         =   chainLength

    W,d,geno4window = window_genotype(geno,input)
    ###START
    y   = y
    X   = fixed.C
    M   = geno4window.genotypes

    current   = Current(input,geno4window,fixed,y)
    output    = QTL.Output(input,geno,fixed)

    mGibbs = GibbsMats(M)
    xGibbs = GibbsMats(X)

    meanVare  = 0.0
    meanVara  = 0.0
    meanPi    = 0.0

    #initial values
    z          = zeros(length(d))

    #return values
    polyVar    = zeros(chainLength)

    nWindows    = round(Int64,nMarkers/windowSize)
    winVarProps = zeros(chainLength,nWindows)

    # MCMC sampling
    for iter = 1:input.chainLength

        current.iter += 1
        # sample fixed effects
        sample_fixed!(xGibbs,current,output)
        # sample marker effects
        sampleEffectsBayesC!(mGibbs,current,output)
        #sample polygenic effect
        samplePolygenicEffects!(current.yCorr,W,d,z,current.varResidual,varPoly)
        # sample residual vairance
        current.varResidual=sample_variance(current.yCorr, geno.nObs, input.nuRes, current.scaleRes)
        # sample marker vairance
        current.varEffect = sample_variance(current.α, current.nLoci, input.dfEffectVar, current.scaleVar)
        # sample polygenic variance
        varPoly = sampleVariance(z./sqrt(d), length(d), dfPoly, scalePoly)
        polyVar[i] = varPoly
        genVar[i] = var(M*current.u + W*z)

        if (input.estimatePi == true)
            current.π = samplePi(current.nLoci, geno.nMarkers)
        end

        wEnd = 0
        for win=1:nWindows
            wStart = wEnd + 1
            wEnd  += windowSize
            wEnd   = (wEnd > nMarkers) ? nMarkers:wEnd
            winVarProps[i,win] = var(X[:,wStart:wEnd]*α[wStart:wEnd])/genVar[i]
        end


       if (iter%outFreq ==0)
       ##yCorr = y - C*β - X*α - W*z # remove rounding errors
        @printf("Iteration %d with %d loci included in the model, mean residual/marker effect %6.3f/%6.3f with mean π= %6.3f.\n",
               iter,current.nLoci, meanVare, meanVara, meanPi)
       end
    end

    ###OUTPUT
    output = Dict()
    output["posterior mean of fixed effects"]         = meanFxdEff
    output["posterior mean of marker effects"]        = meanAlpha
    output["model frequency"]                         = mdlFrq
    output["posterior sample of pi"]                  = pi
    output["posterior sample of scale"]               = scale
    output["posterior sample of genotypic variance"]  = genVar
    output["posterior sample of residual variance"]   = resVar
    output["Proportion of variance explained by genomic window"]   = winVarProps
    output["Posterior sample of polygenic variance"]  = polyVar
    return output

end
