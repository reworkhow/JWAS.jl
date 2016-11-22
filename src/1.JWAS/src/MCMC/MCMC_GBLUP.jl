#y = μ + a + e with var(a)=G
#G = LDL'
#y = μ + Lα +e with var(α)=D*σ² : L orthogonal

function MCMC_GBLUP(nIter,mme,df;
                     sol       =false,
                     outFreq   =1000,
                     output_samples_frequency=0)

   #######################################################
   # Pre-Check
   #######################################################
    if size(mme.mmeRhs)==()
        getMME(mme,df)
    end

    if sol == false #no starting values
       sol = zeros(size(mme.mmeLhs,1))
    else #besure type is Float64
       sol = map(Float64,sol)
    end
    solMean = fill(0.0,size(mme.mmeLhs,1))

    #######################################################
    # PRIORS
    #######################################################
    #prior for residual variance
    vRes        = mme.RNew
    nuRes       = mme.df.residual
    scaleRes    = vRes*(nuRes-2)/nuRes
    meanVare    = 0.0
    #priors for genetic variance (polygenic effects;A) e.g Animal+ Maternal
    if mme.ped != 0
       ν         = mme.df.polygenic
       pedTrmVec = mme.pedTrmVec
       k         = size(pedTrmVec,1)  #2
       νG0       = ν + k
       G0        = inv(mme.GiNew)
       P         = G0*(νG0 - k - 1)
       S         = zeros(Float64,k,k)
       G0Mean    = zeros(Float64,k,k)
    end

    #priors for ~~marker effect~~ genetic variance
    mGibbs      = GibbsMats(mme.M.genotypes)
    nObs        = mGibbs.nrows
    nMarkers    = mGibbs.ncols
    M           = mGibbs.X ./ sqrt(2*mme.M.alleleFreq.*(1-mme.M.alleleFreq))
    #G           = M*M'/nMarkers
    G           = (M*M'+eye(size(M,1))*0.00001)/nMarkers

    eigenG      = eig(G)
    L           = eigenG[2]#eigenvectros
    D           = eigenG[1]#eigenvalues

    dfEffectVar = mme.df.marker #actually for genetic effect here
    vEff        = mme.M.G
    scaleVar    = vEff*(dfEffectVar-2)/dfEffectVar   #scale factor for locus effects
    meanVara    = 0.0
    meanVarg    = 0.0

    α           = zeros(nObs)#starting values for marker effeccts are zeros
    meanAlpha   = zeros(nObs)#vectors to save solutions for marker effects

    #####################################################
    #  WORKING VECTORS (ycor)
    #####################################################
    ##starting values for location parameters(no marker) are sol
    #and adjust y for starting values
    ycorr       = vec(full(mme.ySparse)-mme.X*sol)

    #######################################################
    #  SET UP OUTPUT
    #######################################################
    if output_samples_frequency != 0
      #initialize arrays to save MCMC samples
      num_samples = Int(floor(nIter/output_samples_frequency))
      init_sample_arrays(mme,num_samples)

      out_i = 1
    end

    #######################################################
    # MCMC
    #######################################################
    @showprogress "running MCMC for GBLUP..." for iter=1:nIter

        #####################################
        # 1.1. Non-Marker Location Parameters
        #####################################
        ycorr = ycorr + mme.X*sol
        rhs = mme.X'ycorr

        Gibbs(mme.mmeLhs,sol,rhs,vRes)

        ycorr = ycorr - mme.X*sol
        solMean += (sol - solMean)/iter

        #####################################
        # 1.2 genetic Effects G
        #####################################
        ycorr = ycorr + L*α
        lhs   = 1 + vRes./(vEff*D)
        mean  = L'ycorr./lhs
        α     = mean + randn(nObs).*sqrt(vRes./lhs)
        meanAlpha += (α - meanAlpha)/iter
        ycorr = ycorr - L*α

        ##########################################################################
        # 2.1 Genetic Covariance Matrix (Polygenic Effects)
        ##########################################################################
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
            G0 = rand(InverseWishart(νG0 + q, round(P + S,7))) #ν+q?

            mme.GiOld = copy(mme.GiNew)
            mme.GiNew = inv(G0)

            G0Mean  += (G0  - G0Mean )/iter

            addA(mme) #add Ainverse*lambda
        end
        ##########################################################################
        # 2.2 varainces for (iid) random effects;not required(empty)=>jump out
        ##########################################################################
        sampleVCs(mme,sol)
        addLambdas(mme)

        ###############################################
        # 2.3 Residual Variance
        ###############################################
        mme.ROld = mme.RNew
        vRes     = sample_variance(ycorr, length(ycorr), nuRes, scaleRes)
        mme.RNew = vRes
        meanVare += (vRes - meanVare)/iter

        ###############################################
        # 2.4 Marker Effects Variamce (modified new alpha)
        ###############################################
        vEff  = sample_variance(α./sqrt(D), nObs, dfEffectVar, scaleVar)
        meanVara += (vEff - meanVara)/iter

        varg     = var(L*α)
        meanVarg += (varg - meanVarg)/iter

        ###############################################
        # OUTPUT
        ###############################################
        if output_samples_frequency != 0
          if iter%output_samples_frequency==0
            outputSamples(mme,sol,out_i)
            mme.samples4R[:,out_i]=vRes
            if mme.ped != 0
              mme.samples4G[:,out_i]=vec(G0)
            end
            out_i +=1
          end
        end

        if iter%outFreq==0
            println("\nPosterior means at iteration: ",iter)
            println("Residual variance: ",round(meanVare,6))
            if mme.ped !=0
              println("Polygenic effects covariance matrix \n",round(G0Mean,3))
            end
            if mme.M != 0
              println("Genetic variance (G matrix): ",round(meanVara,6))
              println("Genetic variance (GenSel): ",round(meanVarg,6))
            end
        end
    end

    #######################################################
    # After MCMC
    #######################################################
    output = Dict()
    output["Posterior mean of location parameters"] = [getNames(mme) solMean]
    if output_samples_frequency != 0
        output["MCMC samples for residual variance"]    = mme.samples4R
        if mme.ped != 0
            output["MCMC samples for polygenic effects var-cov parameters"] = mme.samples4G
        end
        for i in  mme.outputSamplesVec
            trmi   = i.term
            trmStr = trmi.trmStr
            output["MCMC samples for: "*trmStr] = [getNames(trmi) i.sampleArray]
        end
        for i in  mme.rndTrmVec
            trmi   = i.term
            trmStr = trmi.trmStr
            output["MCMC samples for: variance of "*trmStr] = i.sampleArray
        end
    end

    output["Posterior mean of genetic variance (G matrix)"]  = meanVara
    output["Posterior mean of genetic variance (GenSel)"]  = meanVarg
    output["Posterior mean of residual variance"] = meanVare

    return output
end
