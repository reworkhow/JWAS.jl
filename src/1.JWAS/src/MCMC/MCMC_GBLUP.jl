#Efficient GBLUP
#y = μ + a + e with var(a)=G**σ²
#G = LDL'
#y = μ + Lα +e with var(α)=D*σ² : L orthogonal

function MCMC_GBLUP(nIter,mme,df;
                     sol       =false,
                     outFreq   =1000,
                     output_samples_frequency=0)

    #############################################################################
    # Pre-Check
    #############################################################################
    #starting values for location parameters(no marker) are sol
    sol,solMean = pre_check(mme,df,sol)

    #######################################################
    # PRIORS
    #######################################################
    #prior for residual variance
    vRes        = mme.RNew
    nuRes       = mme.df.residual
    scaleRes    = vRes*(nuRes-2)/nuRes
    meanVare    = 0.0

    #priors for ~~marker effect~~ genetic variance
    mGibbs      = GibbsMats(mme.M.genotypes)
    nObs,nMarkers,mArray,mpm,M = mGibbs.nrows,mGibbs.ncols,mGibbs.xArray,mGibbs.xpx,mGibbs.X
    M           = mGibbs.X ./ sqrt(2*mme.M.alleleFreq.*(1-mme.M.alleleFreq))
    #G           = M*M'/nMarkers
    G           = (M*M'+eye(size(M,1))*0.00001)/nMarkers

    eigenG      = eig(G)
    L           = eigenG[2]#eigenvectros
    D           = eigenG[1]#eigenvalues

    dfEffectVar = mme.df.marker #actually for genetic effect here
    vEff        = mme.M.G       #genetic variance
    scaleVar    = vEff*(dfEffectVar-2)/dfEffectVar   #scale factor for locus effects
    meanVara    = 0.0
    meanVarg    = 0.0

    α           = zeros(nObs) #starting values for marker effeccts are zeros
    meanAlpha   = zeros(nObs) #vectors to save solutions for marker effects

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
    ############################################################################
    #  WORKING VECTORS (ycor)
    ############################################################################
    #adjust y for starting values
    ycorr       = vec(full(mme.ySparse)-mme.X*sol)

    ############################################################################
    #  SET UP OUTPUT MCMC samples
    ############################################################################
    if output_samples_frequency != 0
      out_i =output_MCMC_samples_setup(mme,nIter,output_samples_frequency,false)
    end

    ############################################################################
    # MCMC
    ############################################################################
    @showprogress "running MCMC for GBLUP..." for iter=1:nIter

        ########################################################################
        # 1.1. Non-Marker Location Parameters
        ########################################################################
        ycorr = ycorr + mme.X*sol
        rhs = mme.X'ycorr

        Gibbs(mme.mmeLhs,sol,rhs,vRes)

        ycorr = ycorr - mme.X*sol
        solMean += (sol - solMean)/iter

        ########################################################################
        # 1.2 Marker Effects (modified new alpha)
        ########################################################################
        ycorr = ycorr + L*α
        lhs   = 1 + vRes./(vEff*D)
        mean1 = L'ycorr./lhs
        α     = mean1 + randn(nObs).*sqrt(vRes./lhs)
        meanAlpha += (α - meanAlpha)/iter
        ycorr = ycorr - L*α

        ########################################################################
        # 2.1 Genetic Covariance Matrix (Polygenic Effects) (variance.jl)
        ########################################################################
        if mme.ped != 0
          G0=sample_variance_pedigree(mme,pedTrmVec,sol,P,S,νG0)
          G0Mean  += (G0  - G0Mean )/iter
        end
        ########################################################################
        # 2.2 varainces for (iid) random effects;not required(empty)=>jump out
        ########################################################################
        sampleVCs(mme,sol)
        addLambdas(mme)
        ########################################################################
        # 2.3 Residual Variance
        ########################################################################
        mme.ROld = mme.RNew
        vRes     = sample_variance(ycorr, length(ycorr), nuRes, scaleRes)
        mme.RNew = vRes
        meanVare += (vRes - meanVare)/iter
        ###############################################
        # 2.4 Marker Effects Variance (modified new alpha)
        ###############################################
        vEff  = sample_variance(α./sqrt(D), nObs, dfEffectVar, scaleVar)
        meanVara += (vEff - meanVara)/iter

        varg     = var(L*α,corrected=false)
        meanVarg += (varg - meanVarg)/iter

        ########################################################################
        # 3.1 Save MCMC samples
        ########################################################################
        if output_samples_frequency != 0 && iter%output_samples_frequency==0
            out_i=output_MCMC_samples(mme,out_i,sol,vRes,(mme.ped!=0?G0:false))
        end
        ########################################################################
        # 3.2 Printout
        ########################################################################
        if iter%outFreq==0
            println("\nPosterior means at iteration: ",iter)
            println("Residual variance: ",round(meanVare,6))
            if mme.ped !=0
              println("Polygenic effects covariance matrix \n",round(G0Mean,3))
            end
            println("Genetic variance (G matrix): ",round(meanVara,6))
            println("Genetic variance (GenSel): ",round(meanVarg,6))
        end
    end

    ############################################################################
    # After MCMC
    ############################################################################
    output=output_result(mme,solMean,output_samples_frequency)
    output["Posterior mean of estimated breeding values"]  = L*meanAlpha
    output["Posterior mean of genetic variance (G matrix)"]= meanVara
    output["Posterior mean of genetic variance (GenSel)"]  = meanVarg
    return output
end
