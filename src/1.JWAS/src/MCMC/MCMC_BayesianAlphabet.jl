################################################################################
#MCMC for RR-BLUP, GBLUP, BayesABC, and conventional (no markers) methods
################################################################################
function MCMC_BayesianAlphabet(mme,df)
    ############################################################################
    chain_length             = mme.MCMCinfo.chain_length
    burnin                   = mme.MCMCinfo.burnin
    output_samples_frequency = mme.MCMCinfo.output_samples_frequency
    output_folder            = mme.MCMCinfo.output_folder
    invweights               = mme.invweights
    update_priors_frequency  = mme.MCMCinfo.update_priors_frequency
    has_categorical_trait    = "categorical"         ∈ mme.traits_type
    has_binary_trait         = "categorical(binary)" ∈ mme.traits_type
    has_censored_trait       = "censored"            ∈ mme.traits_type
    missing_phenotypes       = mme.MCMCinfo.missing_phenotypes
    causal_structure         = mme.causal_structure
    is_multi_trait           = mme.nModels != 1
    is_nnbayes_partial       = mme.nonlinear_function != false && mme.is_fully_connected==false
    is_activation_fcn        = mme.is_activation_fcn
    nonlinear_function       = mme.nonlinear_function
    fast_blocks              = mme.MCMCinfo.fast_blocks
    ############################################################################
    # Categorical Traits (starting values for maker effects defaulting to 0s)
    ############################################################################
    if has_categorical_trait || has_censored_trait || has_binary_trait
        lower_bound,upper_bound,category_obs = categorical_censored_traits_setup!(mme,df) #initialize: mme.thresholds, lower_bound, upper_bound, liability(=mme.ySparse)
    end
    ############################################################################
    # Working Variables
    # 1) samples at current iteration (starting values default to zeros)
    # 2) posterior mean and variance at current iteration (zeros at the beginning)
    # 3) ycorr: phenotypes corrected for all effects
    ############################################################################
    #location parameters
    #mme.sol (its starting values were set in runMCMC)
    mme.solMean, mme.solMean2  = zero(mme.sol),zero(mme.sol)
    #residual variance
    mme.meanVare  = zero(mme.R.val)
    mme.meanVare2 = zero(mme.R.val)

    #polygenic effects (pedigree), e.g, Animal+ Maternal
    if mme.pedTrmVec != 0
       mme.G0Mean,mme.G0Mean2  = zero(mme.Gi),zero(mme.Gi)
    end
    #marker effects
    if mme.M != 0
        for Mi in mme.M
            #Mi.α  (starting values were set in get_genotypes)
            mGibbs    = GibbsMats(Mi.genotypes,invweights,fast_blocks=mme.MCMCinfo.fast_blocks)
            Mi.mArray,Mi.mRinvArray,Mi.mpRinvm = mGibbs.xArray,mGibbs.xRinvArray,mGibbs.xpRinvx
            if fast_blocks != false
                Mi.MArray  = mGibbs.XArray
                Mi.MRinvArray = mGibbs.XRinvArray
                Mi.MpRinvM = mGibbs.XpRinvX
            end

            if Mi.method=="BayesB" #α=β.*δ
                Mi.G.val        = fill(Mi.G.val,Mi.nMarkers) #a scalar in BayesC but a vector in BayeB
            end
            if Mi.method=="BayesL"         #in the MTBayesLasso paper
                if mme.nModels == 1
                    Mi.G.val   /= 8           #mme.M.G.val is the scale Matrix, Sigma
                    Mi.G.scale /= 8
                    gammaDist  = Gamma(1, 8) #8 is the scale parameter of the Gamma distribution (1/8 is the rate parameter)
                else
                    Mi.G.val         /= 4*(Mi.ntraits+1)
                    Mi.G.scale     /= 4*(Mi.ntraits+1)
                    gammaDist     = Gamma((Mi.ntraits+1)/2, 8) #8 (1/8): the scale (rate) parameter
                end
                Mi.gammaArray = rand(gammaDist,Mi.nMarkers)
            end
            if Mi.method=="GBLUP"
                GBLUP_setup(Mi)
            end
            Mi.β                  = [copy(Mi.α[traiti]) for traiti = 1:Mi.ntraits] #partial marker effeccts used in BayesB
            Mi.δ                  = [ones(typeof(Mi.α[traiti][1]),Mi.nMarkers) for traiti = 1:Mi.ntraits] #inclusion indicator for marker effects
            Mi.meanAlpha          = [zero(Mi.α[traiti]) for traiti = 1:Mi.ntraits] #marker effects
            Mi.meanAlpha2         = [zero(Mi.α[traiti]) for traiti = 1:Mi.ntraits] #marker effects
            Mi.meanDelta          = [zero(Mi.δ[traiti]) for traiti = 1:Mi.ntraits] #inclusion indicator for marker effects
            Mi.meanVara           = zero(mme.R.val)  #posterir mean of variance for marker effect
            Mi.meanVara2          = zero(mme.R.val)  #variable to save variance for marker effect
            Mi.meanScaleVara      = zero(mme.R.val) #variable to save Scale parameter for prior of marker effect variance
            Mi.meanScaleVara2     = zero(mme.R.val)  #variable to save Scale parameter for prior of marker effect variance
            if is_multi_trait
                # if is_mega_trait
                if Mi.G.constraint==true
                    Mi.π        = zeros(Mi.ntraits)
                    Mi.mean_pi  = zeros(Mi.ntraits)
                    Mi.mean_pi2 = zeros(Mi.ntraits)
                else
                    Mi.π,Mi.mean_pi,Mi.mean_pi2 = copy(Mi.π),copy(Mi.π),copy(Mi.π)
                    if Mi.estimatePi == true
                      for key in keys(Mi.mean_pi)
                        Mi.mean_pi[key]=0.0
                        Mi.mean_pi2[key]=0.0
                      end
                    end
                    #if methods == "BayesCC"  labels,BigPi,BigPiMean=setPi(Pi)
                end
            else
                Mi.mean_pi,Mi.mean_pi2 = 0.0,0.0      #inclusion probability
            end
        end
    end
    if is_nnbayes_partial
        nnbayes_partial_para_modify3(mme)
    end

    #phenotypes corrected for all effects
    ycorr = vec(Matrix(mme.ySparse)-mme.X*mme.sol)
    if mme.M != 0
        for Mi in mme.M
            for traiti in 1:Mi.ntraits
                if Mi.α[traiti] != zero(Mi.α[traiti])
                    ycorr[(traiti-1)*Mi.nObs+1 : traiti*Mi.nObs] = ycorr[(traiti-1)*Mi.nObs+1 : traiti*Mi.nObs]
                                                                 - Mi.genotypes*Mi.α[traiti]
                end
            end
        end
    end
    ############################################################################
    #More on Multi-Trait
    ############################################################################
    if is_multi_trait
        wArray = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,mme.nModels)
        for traiti = 1:mme.nModels
            startPosi             = (traiti-1)*length(mme.obsID)  + 1
            ptr                   = pointer(ycorr,startPosi)
            wArray[traiti]        = unsafe_wrap(Array,ptr,length(mme.obsID))
        end

        #Starting value for Ri is made based on missing value pattern
        #(imputed phenotypes will not used to compute first mmeRhs)
        Ri         = mkRi(mme,df,invweights)
        dropzeros!(Ri)
    end
    ############################################################################
    # Starting values for SEM
    ############################################################################
    if causal_structure != false
        Y,Λy,causal_structure_outfile = SEM_setup(wArray,causal_structure,mme)
    end
    ############################################################################
    #  SET UP OUTPUT MCMC samples
    ############################################################################
    if output_samples_frequency != 0
        outfile=output_MCMC_samples_setup(mme,chain_length-burnin,
                                          output_samples_frequency,
                                          output_folder*"/MCMC_samples")
    end
    ############################################################################
    # MCMC (starting values for sol (zeros);  mme.RNew; G0 are used)
    ############################################################################
    # # Initialize mme for hmc before Gibbs
    if nonlinear_function != false
        mme.weights_NN    = vcat(mean(mme.ySparse),zeros(mme.nModels))
    end
    if mme.pedTrmVec!=0
        polygenic_pos = findfirst(i -> i.randomType=="A", mme.rndTrmVec)
    end
    @showprogress "running MCMC ..." for iter=1:chain_length
        ########################################################################
        # 0. Categorical and censored traits
        ########################################################################
        if has_categorical_trait || has_censored_trait || has_binary_trait
            sample_liabilities!(mme,ycorr,lower_bound,upper_bound) #update liability(=mme.ySparse) and ycorr
            categorical_trait_sample_threshold!(mme,lower_bound,upper_bound,category_obs) #update mme.thresholds,lower_bound,upper_bound
        end
        ########################################################################
        # 1. Non-Marker Location Parameters
        ########################################################################
        # 1.1 Update Left-hand-side of MME
        if is_multi_trait
            mme.mmeLhs = mme.X'Ri*mme.X #normal equation, Ri is changed
            dropzeros!(mme.mmeLhs)
        end
        addVinv(mme)
        # 1.2 Update Right-hand-side of MME
        if is_multi_trait
            if mme.MCMCinfo.missing_phenotypes==true
              ycorr[:]=sampleMissingResiduals(mme,ycorr)
            end
        end
        ycorr[:] = ycorr + mme.X*mme.sol
        if is_multi_trait
            mme.mmeRhs =  mme.X'Ri*ycorr
        else
            mme.mmeRhs = (invweights == false) ? mme.X'ycorr : mme.X'Diagonal(invweights)*ycorr
        end
        # 1.3 Gibbs sampler
        if is_multi_trait
            Gibbs(mme.mmeLhs,mme.sol,mme.mmeRhs)
        else
            Gibbs(mme.mmeLhs,mme.sol,mme.mmeRhs,mme.R.val)
        end

        ycorr[:] = ycorr - mme.X*mme.sol
        ########################################################################
        # 2. Marker Effects
        ########################################################################
        if mme.M !=0
            for i in 1:length(mme.M)
                Mi=mme.M[i]
                ########################################################################
                # Marker Effects
                ########################################################################
                if Mi.method in ["BayesC","BayesB","BayesA"]
                    locus_effect_variances = (Mi.method == "BayesC" ? fill(Mi.G.val,Mi.nMarkers) : Mi.G.val)
                    if is_multi_trait && !is_nnbayes_partial
                        if Mi.G.constraint==true
                            megaBayesABC!(Mi,wArray,mme.R.val,locus_effect_variances)
                        else
                            if fast_blocks == false
                                MTBayesABC!(Mi,wArray,mme.R.val,locus_effect_variances,mme.nModels)
                            else
                                MTBayesABC_block!(Mi,wArray,mme.R.val,locus_effect_variances)
                            end
                        end
                    elseif is_nnbayes_partial
                        BayesABC!(Mi,wArray[i],mme.R.val[i,i],locus_effect_variances) #this can be parallelized (conflict with others)
                    else
                        if fast_blocks == false
                            BayesABC!(Mi,ycorr,mme.R.val,locus_effect_variances)
                        else
                            BayesABC_block!(Mi,ycorr,mme.R.val,locus_effect_variances)
                        end
                    end
                elseif Mi.method =="RR-BLUP"
                    if is_multi_trait && !is_nnbayes_partial
                        if Mi.G.constraint==true
                            megaBayesC0!(Mi,wArray,mme.R.val)
                        else
                            MTBayesC0!(Mi,wArray,mme.R.val)
                        end
                    elseif is_nnbayes_partial
                        BayesC0!(Mi,wArray[i],mme.R.val[i,i])
                    else
                        BayesC0!(Mi,ycorr,mme.R.val)
                    end
                elseif Mi.method == "BayesL"
                    if is_multi_trait && !is_nnbayes_partial
                        #problem with sampleGammaArray
                        if Mi.G.constraint==true
                            megaBayesL!(Mi,wArray,mme.R.val)
                        else
                            MTBayesL!(Mi,wArray,mme.R.val)
                        end
                    elseif is_nnbayes_partial
                        BayesC0!(Mi,wArray[i],mme.R.val[i,i])
                    else
                        BayesL!(Mi,ycorr,mme.R.val)
                    end
                elseif Mi.method == "GBLUP"
                    if is_multi_trait && !is_nnbayes_partial
                        if Mi.G.constraint==true
                            megaGBLUP!(Mi,wArray,mme.R.val,invweights)
                        else
                            MTGBLUP!(Mi,wArray,ycorr,mme.R.val,invweights)
                        end
                    elseif is_nnbayes_partial
                        GBLUP!(Mi,wArray[i],mme.R.val[i,i],invweights)
                    else
                        GBLUP!(Mi,ycorr,mme.R.val,invweights)
                    end
                end
                ########################################################################
                # Marker Inclusion Probability
                ########################################################################
                if Mi.estimatePi == true
                    if is_multi_trait && !is_nnbayes_partial
                        if Mi.G.constraint==true
                            Mi.π = [samplePi(sum(Mi.δ[i]), Mi.nMarkers) for i in 1:mme.nModels]
                        else
                            samplePi(Mi.δ,Mi.π) #samplePi(deltaArray,Mi.π,labels)
                        end
                    else
                        Mi.π = samplePi(sum(Mi.δ[1]), Mi.nMarkers)
                    end
                end
                ########################################################################
                # Variance of Marker Effects
                ########################################################################
                if Mi.G.estimate_variance == true #methd specific estimate_variance
                    sample_marker_effect_variance(Mi)
                    if mme.MCMCinfo.double_precision == false && Mi.method != "BayesB"
                        Mi.G.val = Float32.(Mi.G.val)
                    end
                end
                ########################################################################
                # Scale Parameter in Priors for Marker Effect Variances
                ########################################################################
                if Mi.G.estimate_scale == true
                    if !is_multi_trait
                        a = size(Mi.G.val,1)*Mi.G.df/2   + 1
                        b = sum(Mi.G.df ./ (2*Mi.G.val)) + 1
                        Mi.G.scale = rand(Gamma(a,1/b))
                    end
                end
            end
        end
        ########################################################################
        # 3. Non-marker Variance Components
        ########################################################################

        ########################################################################
        # 3.1 Variance of Non-marker Random Effects
        # e.g, i.i.d; polygenic effects (pedigree)
        ########################################################################
        if length(mme.rndTrmVec)>0
            if mme.rndTrmVec[1].Gi.estimate_variance == true
                sampleVCs(mme,mme.sol)
            end
        end
        ########################################################################
        # 3.2 Residual Variance
        ########################################################################
        if mme.R.estimate_variance == true
            if is_multi_trait
                mme.R.val = sample_variance(wArray, length(mme.obsID),
                                        mme.R.df, mme.R.scale,
                                        invweights,mme.R.constraint;
                                        binary_trait_index=has_binary_trait ? findall(x->x=="categorical(binary)", mme.traits_type) : false)
                Ri    = kron(inv(mme.R.val),spdiagm(0=>invweights))
            else #single trait
                if !has_categorical_trait && !has_binary_trait # fixed =1 for single categorical/binary trait
                    mme.ROld  = mme.R.val
                    mme.R.val = sample_variance(ycorr,length(ycorr), mme.R.df, mme.R.scale, invweights)
                end
            end
            if mme.MCMCinfo.double_precision == false
                mme.R.val = Float32.(mme.R.val)
            end
        end
        ########################################################################
        # 4. Causal Relationships among Phenotypes (Structure Equation Model)
        ########################################################################
        if is_multi_trait && causal_structure != false
            sample4λ,sample4λ_vec = get_Λ(Y,mme.R.val,ycorr,Λy,mme.ySparse,causal_structure) #no missing phenotypes
        end
        ########################################################################
        # 5. Latent Traits (NNBayes)
        ########################################################################
        if nonlinear_function != false #to update ycorr!
            sample_latent_traits(mme.yobs,mme,ycorr,nonlinear_function)
        end
        ########################################################################
        # 5. Update priors using posteriors (empirical) LATER
        ########################################################################
        if update_priors_frequency !=0 && iter%update_priors_frequency==0
            if mme.M!=0 && methods != "BayesB"
                if is_multi_trait
                    mme.M[1].scale = mme.M[1].meanVara*(mme.M[1].G.df-mme.M[1].ntraits-1)
                else
                    mme.M[1].scale   = mme.M[1].meanVara*(mme.M[1].G.df-2)/mme.M[1].G.df
                end
            end
            if mme.pedTrmVec != 0
                polygenic_pos = findfirst(i -> i.randomType=="A", mme.rndTrmVec)
                mme.scalePed  = mme.G0Mean*(mme.rndTrmVec[polygenic_pos].Gi.df - size(mme.pedTrmVec,1) - 1)
            end
            mme.R.scale  =  mme.meanVare*(mme.R.df-2)/mme.R.df
            println("\n Update priors from posteriors.")
        end
        ########################################################################
        # 3.1 Save MCMC samples
        ########################################################################
        if iter>burnin && (iter-burnin)%output_samples_frequency == 0
            #MCMC samples from posterior distributions
            nsamples       = (iter-burnin)/output_samples_frequency
            output_posterior_mean_variance(mme,nsamples)
            #mean and variance of posterior distribution
            if mme.pedTrmVec!=0
                polygenic_effects_variance = inv(mme.rndTrmVec[polygenic_pos].Gi.val) 
            else
                polygenic_effects_variance=false 
            end
            output_MCMC_samples(mme,mme.R.val,polygenic_effects_variance,outfile)
             if causal_structure != false
                 writedlm(causal_structure_outfile,sample4λ_vec',',')
             end
        end
        ########################################################################
        # 3.2 Printout
        ########################################################################
        if iter%mme.MCMCinfo.printout_frequency==0 && iter>burnin
            println("\nPosterior means at iteration: ",iter)
            println("Residual variance: ",round.(mme.meanVare,digits=6))
        end
    end

    ############################################################################
    # After MCMC
    ############################################################################
    if output_samples_frequency != 0
      for (key,value) in outfile
        close(value)
      end
      if causal_structure != false
        close(causal_structure_outfile)
      end
    end
    if methods == "GBLUP"
        for Mi in mme.M
            mv(output_folder*"/MCMC_samples_marker_effects_variances"*"_"*Mi.name*".txt",
               output_folder*"/MCMC_samples_genetic_variance(REML)"*"_"*Mi.name*".txt")
        end
    end

    output=output_result(mme,output_folder,
                         mme.solMean,mme.meanVare,
                         mme.pedTrmVec!=0 ? mme.G0Mean : false,
                         mme.solMean2,mme.meanVare2,
                         mme.pedTrmVec!=0 ? mme.G0Mean2 : false)
    return output
end
