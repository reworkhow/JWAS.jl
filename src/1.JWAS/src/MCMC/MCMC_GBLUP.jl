#Efficient GBLUP
#y = μ + a + e with mean(a)=0,var(a)=G**σ²
#G = LDL'
#<=>
#y = μ + Lα +e with mean(α)=0,var(α)=D*σ² : L orthogonal

function MCMC_GBLUP(nIter,mme,df;
                    burnin                     = 0,
                    sol                        = false,
                    outFreq                    = 1000,
                    output_samples_frequency   = 0,
                    output_file                = "MCMC_samples")


    #TURN OFF OUTPUT MCMC SAMPLES
    output_samples_frequency=0

    #############################################################################
    # Pre-Check
    #############################################################################
    #starting values for location parameters(no marker) are sol
    solMean     = zero(sol)

    #######################################################
    # PRIORS
    #######################################################
    #prior for residual variance
    vRes        = mme.RNew
    nuRes       = mme.df.residual
    scaleRes    = vRes*(nuRes-2)/nuRes
    meanVare    = 0.0

    #priors for genetic variance (polygenic effects;A) e.g Animal+ Maternal
    if mme.pedTrmVec != 0
       ν         = mme.df.polygenic
       pedTrmVec = mme.pedTrmVec
       k         = size(pedTrmVec,1)  #2
       νG0       = ν + k
       G0        = inv(mme.GiNew)
       P         = G0*(νG0 - k - 1)
       S         = zeros(Float64,k,k)
       G0Mean    = zeros(Float64,k,k)
    end

    #priors for breeding values, genetic variance
    mGibbs      = GibbsMats(mme.M.genotypes)
    nObs,nMarkers,mArray,mpm,M = mGibbs.nrows,mGibbs.ncols,mGibbs.xArray,mGibbs.xpx,mGibbs.X
    M           = mGibbs.X ./ sqrt.(2*mme.M.alleleFreq.*(1 .- mme.M.alleleFreq))
    #G           = M*M'/nMarkers
    G           = (M*M'+Matrix{Float64}(I, size(M,1), size(M,1))*0.00001)/nMarkers

    eigenG      = eigen(G)
    L           = eigenG.vectors#eigenvectros
    D           = eigenG.values#eigenvalues

    dfEffectVar = mme.df.marker                #actually for genetic effect here
    vEff        = mme.M.genetic_variance             #genetic variance
    scaleVar    = vEff*(dfEffectVar-2)/dfEffectVar   #scale factor for locus effects
    meanVara    = 0.0
    meanVarg    = 0.0
    meanh2      = 0.0

    α           = zeros(nObs) #starting values for pseudo breeding values are zeros
    meanAlpha   = zeros(nObs) #vectors to save solutions for pseudo breeding values

    ############################################################################
    #  WORKING VECTORS (ycor)
    ############################################################################
    #adjust y for starting values
    ycorr       = vec(Matrix(mme.ySparse)-mme.X*sol)

    ############################################################################
    #  SET UP OUTPUT MCMC samples
    ############################################################################
    if output_samples_frequency != 0
        outfile=output_MCMC_samples_setup(mme,nIter-burnin,output_samples_frequency,output_file)
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
        if iter > burnin
            solMean += (sol - solMean)/(iter-burnin)
        end
        ########################################################################
        # 1.2 pseudo breeding values (modified new alpha)
        ########################################################################
        ycorr = ycorr + L*α
        lhs   = 1 .+ vRes./(vEff*D)
        mean1 = L'ycorr./lhs
        α     = mean1 + randn(nObs).*sqrt.(vRes./lhs)
        ycorr = ycorr - L*α
        if iter > burnin
            meanAlpha += (α - meanAlpha)/(iter-burnin)
        end
        ########################################################################
        # 2.1 Genetic Covariance Matrix (Polygenic Effects) (variance.jl)
        ########################################################################
        if mme.pedTrmVec != 0
          G0=sample_variance_pedigree(mme,pedTrmVec,sol,P,S,νG0)
          addA(mme)
          if iter > burnin
              G0Mean  += (G0  - G0Mean )/(iter-burnin)
          end
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
        if iter > burnin
            meanVare += (vRes - meanVare)/(iter-burnin)
        end
        ###############################################
        # 2.4 Genomic Genetic Variance (modified new alpha)
        ###############################################
        vEff  = sample_variance(α./sqrt.(D), nObs, dfEffectVar, scaleVar)
        if iter > burnin
            meanVara += (vEff - meanVara)/(iter-burnin)
        end

        if iter > burnin
            varg      = var(L*α,corrected=false)
            meanVarg += (varg - meanVarg)/(iter-burnin)

            h2       = varg/(varg+vRes)
            meanh2  += (h2 - meanh2)/(iter-burnin)
        end

        ########################################################################
        # 3.1 Save MCMC samples
        ########################################################################
        if output_samples_frequency != 0 && iter%output_samples_frequency==0 && iter>burnin
            output_MCMC_samples(mme,sol,vRes,(mme.pedTrmVec!=0 ? G0 : false),0)
        end
        ########################################################################
        # 3.2 Printout
        ########################################################################
        if iter%outFreq==0
            println("\nPosterior means at iteration: ",iter)
            println("Residual variance: ",round(meanVare,digits=6))
            if mme.pedTrmVec !=0
              println("Polygenic effects covariance matrix \n",round.(G0Mean,digits=3))
            end
            println("Genetic variance (G matrix): ",round(meanVara,digits=6))
            println("Genetic variance (GenSel): ",round(meanVarg,digits=6))
            println("Heritability: ",round(meanh2,digits=6))
        end
    end

    ############################################################################
    # After MCMC
    ############################################################################

    #output=output_result(mme,solMean,meanVare,(mme.pedTrmVec!=0 ? G0Mean : false),output_samples_frequency,
    #                    false,false,false,output_file)
    output = Dict()
    output["Posterior mean of estimated breeding values"]  = L*meanAlpha
    output["Posterior mean of genetic variance (G matrix)"]= meanVara
    output["Posterior mean of genetic variance (GenSel)"]  = meanVarg
    output["Posterior mean of heritability"]  = meanh2
    return output
end
