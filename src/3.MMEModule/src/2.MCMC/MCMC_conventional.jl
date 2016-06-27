function MCMC_conventional(nIter,mme,df;
                           sol=false,
                           outFreq=100)
    if size(mme.mmeRhs)==()
        MMEModule.getMME(mme,df)
    end

    if sol == false
        sol=zeros(size(mme.mmeLhs,1))
    end #starting value for sol can be provided

    #######################################################
    # PRIORS
    #######################################################
    #prior for residual variance
    vRes     = mme.RNew
    nuRes    = 4
    scaleRes = vRes*(nuRes-2)/nuRes
    #priors for genetic variance (polygenic effects;A) e.g Animal+ Maternal
    if mme.ped != 0
        ν = 4
        pedTrmVec = mme.pedTrmVec
        k = size(pedTrmVec,1)       #2
        νG0 = ν + k
        G0 = inv(mme.GiNew)
        P = G0*(νG0 - k - 1)
        S = zeros(Float64,k,k)
        G0Mean = zeros(Float64,k,k)
    end

    #####################################################
    # WORKING VECTORS (ycor, saving values)
    #####################################################
    #vectors to save solutions for conventional MME part
    p           = size(mme.mmeRhs,1)
    solMean     = zeros(p)
    #initiate vectors to save samples of MCMC
    initSampleArrays(mme,nIter)
    #variables to save variance for marker effects or residual
    meanVare    = 0.0
    meanVara    = 0.0
    #adjust y for strating values
    ycorr       = vec(full(mme.ySparse)-mme.X*sol)   #starting values for location parameters(no marker) are sol


    #######################################################
    # MCMC
    #######################################################
    for iter=1:nIter

        #sample non-marker part
        ycorr = ycorr + mme.X*sol
        rhs = mme.X'ycorr

        MMEModule.Gibbs(mme.mmeLhs,sol,rhs,vRes)

        ycorr = ycorr - mme.X*sol
        solMean += (sol - solMean)/iter

        #sample variances for polygenic effects
        if mme.ped != 0
            for (i,trmi) = enumerate(pedTrmVec)
                pedTrmi   = mme.modelTermDict[trmi]
                startPosi = pedTrmi.startPos
                endPosi   = startPosi + pedTrmi.nLevels - 1
                for (j,trmj) = enumerate(pedTrmVec)
                    pedTrmj   = mme.modelTermDict[trmj]
                    startPosj = pedTrmj.startPos
                    endPosj   = startPosj + pedTrmj.nLevels - 1
                    S[i,j]    = (sol[startPosi:endPosi]'*mme.Ai*sol[startPosj:endPosj])[1,1]
                end
            end

            pedTrm1 = mme.modelTermDict[pedTrmVec[1]]
            q  = pedTrm1.nLevels
            G0 = rand(InverseWishart(νG0 + q, P + S)) #ν+q?

            mme.GiOld = copy(mme.GiNew)
            mme.GiNew = inv(G0)

            G0Mean  += (G0  - G0Mean )/iter
            mme.genVarSampleArray[iter,:] = vec(G0)

            MMEModule.addA(mme) #add Ainverse*lambda
        end

        #sample varainces for (iid) random effects
        #if not required, if empty jump out
        sampleVCs(mme,sol,iter)
        addLambdas(mme)

        #sample residual variance
        mme.ROld = mme.RNew
        vRes  = sample_variance(ycorr, length(ycorr), nuRes, scaleRes)
        mme.RNew = vRes
        meanVare += (vRes - meanVare)/iter
        mme.resVarSampleArray[iter] = vRes

        #output samples for different effects
        outputSamples(mme,sol,iter)

        if iter%outFreq==0
            println("at sample: ",iter, " with meanVare: ",meanVare)
        end
    end


    ###OUTPUT
    output = Dict()
    output["Posterior Mean of Location Parameters"] = [MMEModule.getNames(mme) solMean]
    output["MCMC samples for residual variance"] = mme.resVarSampleArray

    if mme.ped != 0
        output["MCMC samples for genetic var-cov parameters"] = mme.genVarSampleArray
    end

    for i in  mme.outputSamplesVec
        trmi   = i.term
        trmStr = trmi.trmStr
        output["MCMC samples for: "*trmStr] = i.sampleArray
    end

    for i in  mme.rndTrmVec
        trmi   = i.term
        trmStr = trmi.trmStr
        output["MCMC samples for: variance of "*trmStr] = i.sampleArray
    end

    return output
end
