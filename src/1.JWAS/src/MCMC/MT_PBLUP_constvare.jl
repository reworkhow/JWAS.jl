
#This function is for Multi-trait Pedigree-Based BLUP when Residual Covariance
#Matrix is *CONSTANT*. In this situation, Ri (big), the inverse of residual covariance
#matrix is modify based on phenotype missing patterns such that no imputation for
#missing phenotypes is required.
function MT_MCMC_PBLUP_constvare(nIter,mme,df;
                        burnin                     = 0,
                        estimate_variance          = true,
                        sol                        = false,
                        outFreq                    = 1000,
                        output_samples_frequency   = 0,
                        missing_phenotypes         = false,
                        constraint                 = false,
                        update_priors_frequency    = 0,
                        output_file                = "MCMC_samples")

    ############################################################################
    # Pre-Check
    ############################################################################
    #starting values for location parameters(no marker) are sol
    solMean     = zero(sol)

    ############################################################################
    # PRIORS
    ############################################################################
    #Priors for residual covariance matrix
    ν       = mme.df.residual
    nObs    = size(df,1)
    nTraits = size(mme.lhsVec,1)
    R0      = mme.R
    R0Mean  = zeros(Float64,nTraits,nTraits)
    Ri      = mkRi(mme,df)#create Ri of size (ntrait*nobs)^2, no residual imputation

    #Priors for polygenic effect covariance matrix
    if mme.pedTrmVec != 0
      ν         = mme.df.polygenic
      pedTrmVec = mme.pedTrmVec
      k         = size(pedTrmVec,1)
      νG0       = ν + k
      G0        = inv(mme.Gi)
      P         = G0*(νG0 - k - 1)
      S         = zeros(Float64,k,k)
      G0Mean    = zeros(Float64,k,k)
    end

    ############################################################################
    # SET UP OUTPUT MCMC samples
    ############################################################################
    if output_samples_frequency != 0
          outfile=output_MCMC_samples_setup(mme,nIter-burnin,output_samples_frequency,output_file)
    end

    ############################################################################
    #MCMC
    ############################################################################
    @showprogress "running MCMC for MT-PBLUP" for iter=1:nIter
        ########################################################################
        # 1.1. Non-Marker Location Parameters
        ########################################################################
        Gibbs(mme.mmeLhs,sol,mme.mmeRhs)

        if iter > burnin
            solMean += (sol - solMean)/(iter-burnin)
        end

        ########################################################################
        # 2.1 Residual Covariance Matrix
        ########################################################################
        if iter > burnin
            R0Mean  += (R0  - R0Mean )/(iter-burnin)
        end
        ########################################################################
        # -- LHS and RHS for conventional MME (No Markers)
        # -- Position: between new Ri and new Ai
        ########################################################################
        X          = mme.X
        mme.mmeLhs = X'Ri*X
        mme.mmeRhs = X'Ri*mme.ySparse

        ########################################################################
        # 2.2 Genetic Covariance Matrix (Polygenic Effects)
        ########################################################################
        if mme.pedTrmVec != 0 && estimate_variance == true #will optimize taking advantage of symmetry
            for (i,trmi) = enumerate(pedTrmVec)
                pedTrmi   = mme.modelTermDict[trmi]
                startPosi = pedTrmi.startPos
                endPosi   = startPosi + pedTrmi.nLevels - 1
                for (j,trmj) = enumerate(pedTrmVec)
                    pedTrmj    = mme.modelTermDict[trmj]
                    startPosj  = pedTrmj.startPos
                    endPosj    = startPosj + pedTrmj.nLevels - 1
                    S[i,j]     = (sol[startPosi:endPosi]'*mme.Ai*sol[startPosj:endPosj])[1,1]
                end
            end
            pedTrm1 = mme.modelTermDict[pedTrmVec[1]]
            q = pedTrm1.nLevels

            G0 = rand(InverseWishart(νG0 + q, convert(Array,Symmetric(P + S))))
            mme.Gi = inv(G0)
            addA(mme)
            if iter > burnin
                G0Mean  += (G0  - G0Mean )/(iter-burnin)
            end
        end

        ########################################################################
        # 2.2 varainces for (iid) random effects;not required(empty)=>jump out
        ########################################################################
        if estimate_variance == true
            sampleVCs(mme,sol)
            addLambdas(mme)
        end

        ########################################################################
        # 2.4 Update priors using posteriors (empirical)
        ########################################################################
        if update_priors_frequency !=0 && iter%update_priors_frequency==0
            if mme.pedTrmVec != 0
                P  = G0Mean*(νG0 - k - 1)
            end
            println("\n Update priors from posteriors.")
        end

        ########################################################################
        # OUTPUT
        ########################################################################
        if output_samples_frequency != 0 && (iter-burnin)%output_samples_frequency==0 && iter>burnin
            output_MCMC_samples(mme,sol,R0,(mme.pedTrmVec!=0 ? G0 : false),false,false,false,outfile)
        end

        if iter%outFreq==0 && iter>burnin
            println("\nPosterior means at iteration: ",iter)
            println("Residual covariance matrix: \n",round.(R0Mean,digits=6))

            if mme.pedTrmVec !=0
              println("Polygenic effects covariance matrix \n",round.(G0Mean,digits=6))
            end

            println()
        end
    end

    ############################################################################
    # After MCMC
    ############################################################################
    if output_samples_frequency != 0
      for (key,value) in outfile
        close(value)
      end
    end
    output=output_result(mme,solMean,R0Mean,(mme.pedTrmVec!=0 ? G0Mean : false),output_samples_frequency,
                             false,false,false,false,false,false,output_file)
    return output
end
