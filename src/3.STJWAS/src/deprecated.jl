function getSolJ(mme::MME, df::DataFrame)
    if size(mme.mmeRhs)==()
        getMME(mme,df)
    end
    p = size(mme.mmeRhs,1)
    return [getNames(mme) Jacobi(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,0.3,tolerance=0.000001)]
end

function getSolG(mme::MME, df::DataFrame;outFreq=10)
    if size(mme.mmeRhs)==()
        getMME(mme,df)
    end
    p = size(mme.mmeRhs,1)
    return [getNames(mme) GaussSeidel(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,tolerance=0.000001,output=outFreq)]
end

function getSolGibbs(mme::MME, df::DataFrame;nIter=50000,outFreq=100)
    if size(mme.mmeRhs)==()
        getMME(mme,df)
    end
    p = size(mme.mmeRhs,1)
    return [getNames(mme) Gibbs(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,mme.RNew,nIter,outFreq=outFreq)]
end

function sampleMCMC_novar(nIter,mme,df;sol=zeros(size(mme.mmeRhs,1)),outFreq=100)
    if size(mme.mmeRhs)==()
        MMEModule.getMME(mme,df)
    end
    p = size(mme.mmeRhs,1)
    solMean = zeros(p)

    initSampleArrays(mme,nIter)

    #prior for residual variance
    vRes     = mme.RNew

    #priors for marker effect variance
    if mme.M != 0
        nObs,nLoci = size(mme.M.X)
        α          = zeros(Float64,nLoci)
        meanAlpha  = zeros(Float64,nLoci)
        mArray     = mme.M.xArray
        mpm        = [dot(mme.M.X[:,i],mme.M.X[:,i]) for i=1:size(mme.M.X,2)]
        M          = mme.M.X
        dfEffectVar= 4
        vEff       = mme.M.G/mme.M.mean2pq
        scaleVar   = vEff*(dfEffectVar-2)/dfEffectVar #scale factor for locus effects
    end

    #starting values for marker effects are all zeros
    #starting values for other location parameters are sol
    ycorr    = vec(full(mme.ySparse)-mme.X*sol)
    meanVare = 0.0
    meanVara = 0.0

    for iter=1:nIter

        #sample non-marker part
        ycorr = ycorr + mme.X*sol
        rhs = mme.X'ycorr

        MMEModule.Gibbs(mme.mmeLhs,sol,rhs,vRes)

        ycorr = ycorr - mme.X*sol
        solMean += (sol - solMean)/iter

        #sample marker effects
        if mme.M!=0
            MMEModule.sample_effects_ycorr!(M,mArray,mpm,ycorr,α,meanAlpha,vRes,vEff,iter)
        end

        #output samples for different effects
        outputSamples(mme,sol,iter)

        if iter%outFreq==0
            println("at sample: ",iter)
        end
    end


    ###OUTPUT
    output = Dict()
    output["Posterior Mean of Location Parameters"] = [MMEModule.getNames(mme) solMean]

    if mme.M != 0
        output["Posterior Mean of Marker Effects"] = meanAlpha
    end

    for i in  mme.outputSamplesVec
        trmi   = i.term
        trmStr = trmi.trmStr
        output["MCMCSamples: "*trmStr] = i.sampleArray
    end

    return output
end
