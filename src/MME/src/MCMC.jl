function sampleMCMC(nIter,mme,df;sol=zeros(size(mme.mmeRhs,1)),outFreq=100)
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
    output["posteriorMeanLocationParms"] = [MMEModule.getNames(mme) solMean]
   
    if mme.M != 0
        output["posteriorMeanMarkerEffects"] = meanAlpha
    end
    
    for i in  mme.outputSamplesVec
        trmi   = i.term
        trmStr = trmi.trmStr
        output["MCMCSamples: "*trmStr] = i.sampleArray
    end
        
    return output
end


function sampleMCMCv(nIter,mme,df;sol=zeros(size(mme.mmeRhs,1)),outFreq=100)
    if size(mme.mmeRhs)==() 
        MMEModule.getMME(mme,df)
    end
    p = size(mme.mmeRhs,1)
    solMean = zeros(p)
    
    initSampleArrays(mme,nIter)
    
    #prior for residual variance
    vRes     = mme.RNew
    nuRes    = 4
    scaleRes = vRes*(nuRes-2)/nuRes 

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
    
    #priors for genetic variance (polygenic effects;A)
    #e.g Animal"&"Animal*Age"  
    G0Mean   = 0.0 #only for println
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
        vRes  = sampleVariance(ycorr, length(ycorr), nuRes, scaleRes)
        mme.RNew = vRes
        
        meanVare += (vRes - meanVare)/iter
        mme.resVarSampleArray[iter] = vRes       
        
        #sample marker effect variance
        if mme.M != 0
            vEff  = sampleVariance(α, nLoci, dfEffectVar, scaleVar)
            meanVara += (vEff - meanVara)/iter
        end

        #output samples for different effects
        outputSamples(mme,sol,iter)

        if iter%outFreq==0
            println("at sample: ",iter, " with meanVare: ",meanVare," meanVara ",meanVara, " and G0Mean: ",G0Mean)
        end
    end
    
    
    ###OUTPUT
    output = Dict()
    output["posteriorMeanLocationParms"] = [MMEModule.getNames(mme) solMean]
    output["MCMCSamples for residual variance"] = mme.resVarSampleArray
    
    if mme.ped != 0
        output["MCMCSamples for genetic var-cov parameters"] = mme.genVarSampleArray
    end
    
    if mme.M != 0
        output["posteriorMeanMarkerEffects"] = meanAlpha
    end
    
    for i in  mme.outputSamplesVec
        trmi   = i.term
        trmStr = trmi.trmStr
        output["MCMCSamples: "*trmStr] = i.sampleArray
    end
    
    for i in  mme.rndTrmVec
        trmi   = i.term
        trmStr = trmi.trmStr
        output["MCMCSamples for variance of :"*trmStr] = i.sampleArray
    end
        
    return output
end
