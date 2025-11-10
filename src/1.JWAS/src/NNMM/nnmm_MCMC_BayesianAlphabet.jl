################################################################################
#MCMC for RR-BLUP, GBLUP, BayesABC, and conventional (no markers) methods
################################################################################
function nnmm_MCMC_BayesianAlphabet(mme1,df1,mme2,df2)
    
    #1->2:
    ############################################################################
    chain_length             = mme1.MCMCinfo.chain_length
    burnin                   = mme1.MCMCinfo.burnin
    output_samples_frequency = mme1.MCMCinfo.output_samples_frequency
    output_folder            = mme1.MCMCinfo.output_folder
    
    invweights1               = mme1.invweights
    update_priors_frequency1  = mme1.MCMCinfo.update_priors_frequency
    missing_phenotypes1       = mme1.MCMCinfo.missing_phenotypes
    is_multi_trait1           = mme1.nModels != 1
    is_nnbayes_partial       = mme1.nonlinear_function != false && mme1.is_fully_connected==false
    is_activation_fcn        = mme1.is_activation_fcn
    nonlinear_function       = mme1.nonlinear_function

    ############################################################################
    # Working Variables
    # 1) samples at current iteration (starting values default to zeros)
    # 2) posterior mean and variance at current iteration (zeros at the beginning)
    # 3) ycorr: phenotypes corrected for all effects
    ############################################################################
    #location parameters
    #mme.sol (its starting values were set in runMCMC)
    mme1.solMean, mme1.solMean2  = zero(mme1.sol),zero(mme1.sol)
    #residual variance
    mme1.meanVare  = zero(mme1.R.val)
    mme1.meanVare2 = zero(mme1.R.val)

    #polygenic effects (pedigree), e.g, Animal+ Maternal
    if mme1.pedTrmVec != 0
       mme1.G0Mean,mme1.G0Mean2  = zero(mme1.Gi),zero(mme1.Gi)
    end
    #marker effects
    if mme1.M != 0
        for Mi in mme1.M
            #Mi.α  (starting values were set in get_genotypes)
            mGibbs    = GibbsMats(Mi.genotypes,invweights1)
            Mi.mArray,Mi.mRinvArray,Mi.mpRinvm  = mGibbs.xArray,mGibbs.xRinvArray,mGibbs.xpRinvx

            if Mi.method=="BayesB" #α=β.*δ
                Mi.G.val        = fill(Mi.G.val,Mi.nMarkers) #a scalar in BayesC but a vector in BayeB
            end
            if Mi.method=="BayesL"         #in the MTBayesLasso paper
                if mme1.nModels == 1
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
            Mi.meanVara           = zero(mme1.R.val)  #posterir mean of variance for marker effect
            Mi.meanVara2          = zero(mme1.R.val)  #variable to save variance for marker effect
            Mi.meanScaleVara      = zero(mme1.R.val) #variable to save Scale parameter for prior of marker effect variance
            Mi.meanScaleVara2     = zero(mme1.R.val)  #variable to save Scale parameter for prior of marker effect variance
            if is_multi_trait1
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
        nnbayes_partial_para_modify3(mme1)
    end

    #phenotypes corrected for all effects
    ycorr1 = vec(Matrix(mme1.ySparse)-mme1.X*mme1.sol)
    if mme1.M != 0
        for Mi in mme1.M
            for traiti in 1:Mi.ntraits
                if Mi.α[traiti] != zero(Mi.α[traiti])
                    ycorr1[(traiti-1)*Mi.nObs+1 : traiti*Mi.nObs] = ycorr1[(traiti-1)*Mi.nObs+1 : traiti*Mi.nObs]
                                                                 - Mi.genotypes*Mi.α[traiti]
                end
            end
        end
    end
    ############################################################################
    #More on Multi-Trait
    ############################################################################
    if is_multi_trait1
        wArray1 = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,mme1.nModels)
        for traiti = 1:mme1.nModels
            startPosi             = (traiti-1)*length(mme1.obsID)  + 1
            ptr                   = pointer(ycorr1,startPosi)
            wArray1[traiti]        = unsafe_wrap(Array,ptr,length(mme1.obsID))
        end

        #Starting value for Ri is made based on missing value pattern
        #(imputed phenotypes will not used to compute first mmeRhs)
        Ri         = mkRi(mme1,df1,invweights1)
        dropzeros!(Ri)
    end
    ############################################################################
    #  SET UP OUTPUT MCMC samples
    ############################################################################
    if output_samples_frequency != 0
        outfile=output_MCMC_samples_setup(mme1,chain_length-burnin,
                                          output_samples_frequency,
                                          output_folder*"/MCMC_samples")
    end
    ############################################################################
    # MCMC (starting values for sol (zeros);  mme.RNew; G0 are used)
    ############################################################################
    # # Initialize mme for hmc before Gibbs
    if nonlinear_function != false
        # mme1.weights_NN    = vcat(mean(mme1.ySparse),zeros(mme1.nModels))
        mme1.weights_NN    = zeros(mme1.nModels) # intercept is not included
    end
    if mme1.pedTrmVec!=0
        polygenic_pos = findfirst(i -> i.randomType=="A", mme1.rndTrmVec)
    end




    #2->3:
    ############################################################################
    invweights2               = mme2.invweights
    update_priors_frequency2  = mme2.MCMCinfo.update_priors_frequency
    missing_phenotypes2       = mme2.MCMCinfo.missing_phenotypes
    is_multi_trait2           = mme2.nModels != 1
    # has_categorical_trait    = "categorical"         ∈ mme.traits_type
    # has_binary_trait         = "categorical(binary)" ∈ mme.traits_type
    # has_censored_trait       = "censored"            ∈ mme.traits_type
    has_categorical_trait = false
    has_binary_trait = false
    has_censored_trait = false

    ############################################################################
    # Working Variables
    # 1) samples at current iteration (starting values default to zeros)
    # 2) posterior mean and variance at current iteration (zeros at the beginning)
    # 3) ycorr: phenotypes corrected for all effects
    ############################################################################
    #location parameters
    #mme.sol (its starting values were set in runMCMC)
    mme2.solMean, mme2.solMean2  = zero(mme2.sol),zero(mme2.sol)
    #residual variance
    mme2.meanVare  = zero(mme2.R.val)
    mme2.meanVare2 = zero(mme2.R.val)

    #polygenic effects (pedigree), e.g, Animal+ Maternal
    if mme2.pedTrmVec != 0
       mme2.G0Mean,mme2.G0Mean2  = zero(mme2.Gi),zero(mme2.Gi)
    end
    #marker effects
    if mme2.M != 0
        for Mi in mme2.M
            #Mi.α  (starting values were set in get_genotypes)
            Mi_genotypes = convert(Matrix{Float32}, Mi.aligned_omics_w_phenotype)
            mGibbs    = GibbsMats(Mi_genotypes,invweights2)
            Mi.mArray,Mi.mRinvArray,Mi.mpRinvm  = mGibbs.xArray,mGibbs.xRinvArray,mGibbs.xpRinvx

            if Mi.method=="BayesB" #α=β.*δ
                Mi.G.val        = fill(Mi.G.val,Mi.nFeatures) #a scalar in BayesC but a vector in BayeB
            end
            if Mi.method=="BayesL"         #in the MTBayesLasso paper
                if mme2.nModels == 1
                    Mi.G.val   /= 8           #mme.M.G.val is the scale Matrix, Sigma
                    Mi.G.scale /= 8
                    gammaDist  = Gamma(1, 8) #8 is the scale parameter of the Gamma distribution (1/8 is the rate parameter)
                else
                    Mi.G.val         /= 4*(Mi.ntraits+1)
                    Mi.G.scale     /= 4*(Mi.ntraits+1)
                    gammaDist     = Gamma((Mi.ntraits+1)/2, 8) #8 (1/8): the scale (rate) parameter
                end
                Mi.gammaArray = rand(gammaDist,Mi.nFeatures)
            end
            if Mi.method=="GBLUP"
                GBLUP_setup(Mi)
            end
            Mi.β                  = [copy(Mi.α[traiti]) for traiti = 1:Mi.ntraits] #partial marker effeccts used in BayesB
            Mi.δ                  = [ones(typeof(Mi.α[traiti][1]),Mi.nFeatures) for traiti = 1:Mi.ntraits] #inclusion indicator for marker effects
            Mi.meanAlpha          = [zero(Mi.α[traiti]) for traiti = 1:Mi.ntraits] #marker effects
            Mi.meanAlpha2         = [zero(Mi.α[traiti]) for traiti = 1:Mi.ntraits] #marker effects
            Mi.meanDelta          = [zero(Mi.δ[traiti]) for traiti = 1:Mi.ntraits] #inclusion indicator for marker effects
            Mi.meanVara           = zero(mme2.R.val)  #posterir mean of variance for marker effect
            Mi.meanVara2          = zero(mme2.R.val)  #variable to save variance for marker effect
            Mi.meanScaleVara      = zero(mme2.R.val) #variable to save Scale parameter for prior of marker effect variance
            Mi.meanScaleVara2     = zero(mme2.R.val)  #variable to save Scale parameter for prior of marker effect variance
            if is_multi_trait2
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

    #phenotypes corrected for all effects
    ycorr2 = vec(Matrix(mme2.ySparse)-mme2.X*mme2.sol) #length of ycorr2 is #individuals with non-missing yobs
    if mme2.M != 0
        for Mi in mme2.M
            Mi_genotypes = mme2.M[1].aligned_omics_w_phenotype
            for traiti in 1:Mi.ntraits
                if Mi.α[traiti] != zero(Mi.α[traiti])  
                    ycorr2[(traiti-1)*Mi.aligned_nObs_w_phenotype+1 : traiti*Mi.aligned_nObs_w_phenotype] = ycorr2[(traiti-1)*Mi.aligned_nObs_w_phenotype+1 : traiti*Mi.aligned_nObs_w_phenotype]
                                                                 - Mi_genotypes*Mi.α[traiti]
                end
            end
        end
    end

    ############################################################################
    #More on Multi-Trait
    ############################################################################
    if is_multi_trait2
        error("MT 2->3 is not supported for now")
        wArray2 = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,mme2.nModels)
        for traiti = 1:mme2.nModels
            startPosi             = (traiti-1)*length(mme2.obsID)  + 1
            ptr                   = pointer(ycorr2,startPosi)
            wArray2[traiti]        = unsafe_wrap(Array,ptr,length(mme2.obsID))
        end

        #Starting value for Ri is made based on missing value pattern
        #(imputed phenotypes will not used to compute first mmeRhs)
        Ri         = mkRi(mme2,df2,invweights2)
        dropzeros!(Ri)
    end
    ############################################################################
    #  SET UP OUTPUT MCMC samples
    ############################################################################
    # if output_samples_frequency != 0
    #     outfile=output_MCMC_samples_setup(mme,chain_length-burnin,
    #                                       output_samples_frequency,
    #                                       output_folder*"/MCMC_samples")
    # end
    ############################################################################
    # MCMC (starting values for sol (zeros);  mme.RNew; G0 are used)
    ############################################################################

    if mme2.pedTrmVec!=0
        polygenic_pos = findfirst(i -> i.randomType=="A", mme2.rndTrmVec)
    end


    @showprogress "running MCMC ..." for iter=1:chain_length
        
        #1->2
        ########################################################################
        # 1. Non-Marker Location Parameters
        ########################################################################
        # 1.1 Update Left-hand-side of MME
        if is_multi_trait1
            mme1.mmeLhs = mme1.X'Ri*mme1.X #normal equation, Ri is changed
            dropzeros!(mme1.mmeLhs)
        end
        addVinv(mme1)
        # 1.2 Update Right-hand-side of MME
        if is_multi_trait1
            if mme1.MCMCinfo.missing_phenotypes==true
              ycorr1[:]=sampleMissingResiduals(mme1,ycorr1)
            end
        end
        ycorr1[:] = ycorr1 + mme1.X*mme1.sol
        if is_multi_trait1
            mme1.mmeRhs =  mme1.X'Ri*ycorr1
        else
            mme1.mmeRhs = (invweights1 == false) ? mme1.X'ycorr1 : mme1.X'Diagonal(invweights1)*ycorr1
        end
        # 1.3 Gibbs sampler
        if is_multi_trait1
            Gibbs(mme1.mmeLhs,mme1.sol,mme1.mmeRhs)
        else
            Gibbs(mme1.mmeLhs,mme1.sol,mme1.mmeRhs,mme1.R.val)
        end

        ycorr1[:] = ycorr1 - mme1.X*mme1.sol
        ########################################################################
        # 2. Marker Effects
        ########################################################################
        if mme1.M !=0
            for i in 1:length(mme1.M)
                Mi=mme1.M[i]
                ########################################################################
                # Marker Effects
                ########################################################################
                if Mi.method in ["BayesC","BayesB","BayesA"]
                    locus_effect_variances = (Mi.method == "BayesC" ? fill(Mi.G.val,Mi.nMarkers) : Mi.G.val)
                    if is_multi_trait1 && !is_nnbayes_partial
                        if Mi.G.constraint==true
                            megaBayesABC!(Mi,wArray1,mme1.R.val,locus_effect_variances)
                        else
                            MTBayesABC!(Mi,wArray1,mme1.R.val,locus_effect_variances,mme1.nModels)
                        end
                    elseif is_nnbayes_partial
                        BayesABC!(Mi,wArray1[i],mme1.R.val[i,i],locus_effect_variances) #this can be parallelized (conflict with others)
                    else
                        BayesABC!(Mi,ycorr1,mme1.R.val,locus_effect_variances)
                    end
                elseif Mi.method =="RR-BLUP"
                    if is_multi_trait1 && !is_nnbayes_partial
                        if Mi.G.constraint==true
                            megaBayesC0!(Mi,wArray1,mme1.R.val)
                        else
                            MTBayesC0!(Mi,wArray1,mme1.R.val)
                        end
                    elseif is_nnbayes_partial
                        BayesC0!(Mi,wArray1[i],mme1.R.val[i,i])
                    else
                        BayesC0!(Mi,ycorr1,mme1.R.val)
                    end
                elseif Mi.method == "BayesL"
                    if is_multi_trait1 && !is_nnbayes_partial
                        #problem with sampleGammaArray
                        if Mi.G.constraint==true
                            megaBayesL!(Mi,wArray1,mme1.R.val)
                        else
                            MTBayesL!(Mi,wArray1,mme1.R.val)
                        end
                    elseif is_nnbayes_partial
                        BayesC0!(Mi,wArray1[i],mme1.R.val[i,i])
                    else
                        BayesL!(Mi,ycorr1,mme1.R.val)
                    end
                elseif Mi.method == "GBLUP"
                    if is_multi_trait1 && !is_nnbayes_partial
                        if Mi.G.constraint==true
                            megaGBLUP!(Mi,wArray1,mme1.R.val,invweights1)
                        else
                            MTGBLUP!(Mi,wArray1,ycorr1,mme1.R.val,invweights1)
                        end
                    elseif is_nnbayes_partial
                        GBLUP!(Mi,wArray1[i],mme.R.val[i,i],invweights)
                    else
                        GBLUP!(Mi,ycorr1,mme1.R.val,invweights1)
                    end
                end
                ########################################################################
                # Marker Inclusion Probability
                ########################################################################
                if Mi.estimatePi == true
                    if is_multi_trait1 && !is_nnbayes_partial
                        if Mi.G.constraint==true
                            Mi.π = [samplePi(sum(Mi.δ[i]), Mi.nMarkers) for i in 1:mme1.nModels]
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
                    if mme1.MCMCinfo.double_precision == false && Mi.method != "BayesB"
                        Mi.G.val = Float32.(Mi.G.val)
                    end
                end
                ########################################################################
                # Scale Parameter in Priors for Marker Effect Variances
                ########################################################################
                if Mi.G.estimate_scale == true
                    if !is_multi_trait1
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
        if length(mme1.rndTrmVec)>0
            if mme1.rndTrmVec[1].Gi.estimate_variance == true
                sampleVCs(mme1,mme1.sol)
            end
        end
        ########################################################################
        # 3.2 Residual Variance
        ########################################################################
        if mme1.R.estimate_variance == true
            has_binary_trait=false #NNMM does not support binary traits now
            if is_multi_trait1
                mme1.R.val = sample_variance(wArray1, length(mme1.obsID),
                                        mme1.R.df, mme1.R.scale,
                                        invweights1,mme1.R.constraint;
                                        binary_trait_index=has_binary_trait ? findall(x->x=="categorical(binary)", mme1.traits_type) : false)
                Ri    = kron(inv(mme1.R.val),spdiagm(0=>invweights1))
            else #single trait
                if !has_categorical_trait && !has_binary_trait # fixed =1 for single categorical/binary trait
                    mme1.ROld  = mme1.R.val
                    mme1.R.val = sample_variance(ycorr1,length(ycorr1), mme1.R.df, mme1.R.scale, invweights1)
                end
            end
            if mme1.MCMCinfo.double_precision == false
                mme1.R.val = Float32.(mme1.R.val)
            end
        end

        ########################################################################
        # 5. Latent Traits (NNBayes) 
        # to update ycorr1!
        # from sample_latent_traits(yobs,mme,ycorr,nonlinear_function)
        ########################################################################
        yobs = mme1.yobs
        ylats_old = mme1.ySparse          # current values of each latent trait; [trait_1_obs;trait_2_obs;...]
        # μ_ylats   = mme1.ySparse - ycorr1 # mean of each latent trait, [trait_1_obs-residuals;trait_2_obs-residuals;...]
        #                                 # = vcat(getEBV(mme,1).+mme.sol[1],getEBV(mme,2).+mme.sol[2]))
        μ_ylats   = vcat((getEBV(mme1, i) for i in 1:mme1.nModels)...)  #cannot use "mme1.ySparse - ycorr1" becuase there may be other ramdon/fixed model terms                       
        σ2_yobs      = mme1.σ2_yobs      # residual variance of yobs (scalar)
        σ2_weightsNN = mme1.σ2_weightsNN # variance of nn weights between middle and output layers
    
        #reshape the vector to nind X ntraits
        nobs, ntraits = length(mme1.obsID), mme1.nModels
        ylats_old     = reshape(ylats_old,nobs,ntraits) #Tianjing's mme.Z
        ylats_old2    = copy(ylats_old) #save original omics data before updating
        μ_ylats       = reshape(μ_ylats,nobs,ntraits)
        ycorr_reshape        = reshape(ycorr1,nobs,ntraits)
    
        #sample latent trait
        incomplete_omics = mme1.incomplete_omics #indicator for ind with incomplete omics
        if mme1.is_activation_fcn == true #Neural Network with activation function
            if sum(incomplete_omics) != 0   #at least 1 ind with incomplete omics
                #step 1. sample latent trait (only for individuals with incomplete omics data
                Z  = map(Float32,mkmat_incidence_factor(mme2.obsID,mme2.M[1].obsID))
                ylats_new = hmc_one_iteration(10,0.1,ylats_old[incomplete_omics,:],yobs[incomplete_omics],mme1.weights_NN,mme1.R.val,σ2_yobs,ycorr_reshape[incomplete_omics,:],nonlinear_function,ycorr2[BitVector(Z*incomplete_omics),:])
            end
        else  #user-defined function, MH
            candidates       = μ_ylats+randn(size(μ_ylats))  #candidate samples
            if nonlinear_function == "Neural Network (MH)"
                error("not supported for now")
                μ_yobs_candidate = [ones(nobs) nonlinear_function.(candidates)]*weights
                μ_yobs_current   = X*weights
            else #user-defined non-linear function
                μ_yobs_candidate = nonlinear_function.(Tuple([view(candidates,:,i) for i in 1:ntraits])...)
                μ_yobs_current   = nonlinear_function.(Tuple([view(ylats_old,:,i) for i in 1:ntraits])...)
            end
            llh_current      = -0.5*(yobs - μ_yobs_current ).^2/σ2_yobs
            llh_candidate    = -0.5*(yobs - μ_yobs_candidate).^2/σ2_yobs
            mhRatio          = exp.(llh_candidate - llh_current)
            updateus         = rand(nobs) .< mhRatio
            ylats_new        = candidates.*updateus + ylats_old.*(.!updateus)
        end
        
        #step 2. update ylats with sampled latent traits
        ylats_old[incomplete_omics,:] = ylats_new
        
        #step 3. for individuals with partial omics data, put back the partial real omics.
        ylats_old[mme1.missingPattern] .= ylats_old2[mme1.missingPattern]
    
        #update ylats (i.e., omics)
        mme1.ySparse = vec(ylats_old) #omics data, note that ylats_old is the updated omics data, the old omics data is in ylats_old2

        #update ycorr1 (for 1->2)
        # ycorr1[:]    = mme1.ySparse - vec(μ_ylats) # =(ylats_new - ylats_old) + ycorr: update residuls (ycorr)
        ycorr1[:] += mme1.ySparse - vec(ylats_old2)

        #update data for 2->3
        #update omics data
        mme2.M[1].data[!,mme2.M[1].featureID] = ylats_old
        #update aligned transformed omics data (g(z))
        align_transformed_omics_with_phenotypes(mme2,nonlinear_function)
        #update ycorr2 
        ycorr2[:] = vec(Matrix(mme2.ySparse)-mme2.X*mme2.sol)
        if mme2.M != 0
            for Mi in mme2.M
                Mi_genotypes = mme2.M[1].aligned_omics_w_phenotype
                for traiti in 1:Mi.ntraits
                    if Mi.α[traiti] != zero(Mi.α[traiti])  
                        ycorr2[(traiti-1)*Mi.aligned_nObs_w_phenotype+1 : traiti*Mi.aligned_nObs_w_phenotype] = ycorr2[(traiti-1)*Mi.aligned_nObs_w_phenotype+1 : traiti*Mi.aligned_nObs_w_phenotype]
                                                                    - Mi_genotypes*Mi.α[traiti]
                    end
                end
            end
        end
        #update Mi.mArray, Mi.mRinvArray, Mi.mpRinvx for 2->3
        Mi_genotypes = convert(Matrix{Float32}, mme2.M[1].aligned_omics_w_phenotype)
        mGibbs    = GibbsMats(Mi_genotypes,invweights2)
        mme2.M[1].mArray, mme2.M[1].mRinvArray, mme2.M[1].mpRinvm  = mGibbs.xArray, mGibbs.xRinvArray, mGibbs.xpRinvx

        #2->3
        ########################################################################
        # 1. Non-Marker Location Parameters
        ########################################################################
        # 1.1 Update Left-hand-side of MME
        if is_multi_trait2
            error("MT 2->3 is not supported for now")
            mme2.mmeLhs = mme2.X'Ri*mme2.X #normal equation, Ri is changed
            dropzeros!(mme2.mmeLhs)
        end
        addVinv(mme2)
        # 1.2 Update Right-hand-side of MME
        if is_multi_trait2
            if mme2.MCMCinfo.missing_phenotypes==true
              ycorr2[:]=sampleMissingResiduals(mme2,ycorr2)
            end
        end
        ycorr2[:] = ycorr2 + mme2.X*mme2.sol
        if is_multi_trait2
            mme2.mmeRhs =  mme2.X'Ri*ycorr2
        else
            mme2.mmeRhs = (invweights2 == false) ? mme2.X'ycorr2 : mme2.X'Diagonal(invweights2)*ycorr2
        end
        # 1.3 Gibbs sampler
        if is_multi_trait2
            Gibbs(mme2.mmeLhs,mme2.sol,mme2.mmeRhs)
        else
            Gibbs(mme2.mmeLhs,mme2.sol,mme2.mmeRhs,mme2.R.val)
        end

        ycorr2[:] = ycorr2 - mme2.X*mme2.sol
        ########################################################################
        # 2. Marker Effects
        ########################################################################
        if mme2.M !=0
            for i in 1:length(mme2.M)
                Mi=mme2.M[i]
                ########################################################################
                # Marker Effects
                ########################################################################
                if Mi.method in ["BayesC","BayesB","BayesA"]
                    locus_effect_variances = (Mi.method == "BayesC" ? fill(Mi.G.val,Mi.nFeatures) : Mi.G.val)
                    if is_multi_trait2 && !is_nnbayes_partial
                        if Mi.G.constraint==true
                            megaBayesABC!(Mi,wArray2,mme2.R.val,locus_effect_variances)
                        else
                            MTBayesABC!(Mi,wArray2,mme2.R.val,locus_effect_variances,mme2.nModels)
                        end
                    elseif is_nnbayes_partial
                        BayesABC!(Mi,wArray2[i],mme2.R.val[i,i],locus_effect_variances) #this can be parallelized (conflict with others)
                    else
                        BayesABC!(Mi,ycorr2,mme2.R.val,locus_effect_variances)
                    end
                elseif Mi.method =="RR-BLUP"
                    if is_multi_trait2 && !is_nnbayes_partial
                        if Mi.G.constraint==true
                            megaBayesC0!(Mi,wArray2,mme2.R.val)
                        else
                            MTBayesC0!(Mi,wArray2,mme2.R.val)
                        end
                    elseif is_nnbayes_partial
                        BayesC0!(Mi,wArray2[i],mme2.R.val[i,i])
                    else
                        BayesC0!(Mi,ycorr2,mme2.R.val)
                    end
                elseif Mi.method == "BayesL"
                    if is_multi_trait2 && !is_nnbayes_partial
                        #problem with sampleGammaArray
                        if Mi.G.constraint==true
                            megaBayesL!(Mi,wArray2,mme2.R.val)
                        else
                            MTBayesL!(Mi,wArray2,mme2.R.val)
                        end
                    elseif is_nnbayes_partial
                        BayesC0!(Mi,wArray2[i],mme2.R.val[i,i])
                    else
                        BayesL!(Mi,ycorr2,mme2.R.val)
                    end
                elseif Mi.method == "GBLUP"
                    if is_multi_trait2 && !is_nnbayes_partial
                        if Mi.G.constraint==true
                            megaGBLUP!(Mi,wArray2,mme2.R.val,invweights2)
                        else
                            MTGBLUP!(Mi,wArray2,ycorr2,mme2.R.val,invweights2)
                        end
                    elseif is_nnbayes_partial
                        GBLUP!(Mi,wArray2[i],mme2.R.val[i,i],invweights2)
                    else
                        GBLUP!(Mi,ycorr2,mme2.R.val,invweights2)
                    end
                end
                ########################################################################
                # Marker Inclusion Probability
                ########################################################################
                if Mi.estimatePi == true
                    if is_multi_trait2 && !is_nnbayes_partial
                        if Mi.G.constraint==true
                            Mi.π = [samplePi(sum(Mi.δ[i]), Mi.nFeatures) for i in 1:mme2.nModels]
                        else
                            samplePi(Mi.δ,Mi.π) #samplePi(deltaArray,Mi.π,labels)
                        end
                    else
                        Mi.π = samplePi(sum(Mi.δ[1]), Mi.nFeatures)
                    end
                end
                ########################################################################
                # Variance of Marker Effects
                ########################################################################
                if Mi.G.estimate_variance == true #methd specific estimate_variance
                    sample_marker_effect_variance(Mi)
                    if mme2.MCMCinfo.double_precision == false && Mi.method != "BayesB"
                        Mi.G.val = Float32.(Mi.G.val)
                    end
                end
                ########################################################################
                # Scale Parameter in Priors for Marker Effect Variances
                ########################################################################
                if Mi.G.estimate_scale == true
                    if !is_multi_trait1
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
        if length(mme2.rndTrmVec)>0
            if mme2.rndTrmVec[1].Gi.estimate_variance == true
                sampleVCs(mme2,mme2.sol)
            end
        end
        ########################################################################
        # 3.2 Residual Variance
        ########################################################################
        if mme2.R.estimate_variance == true
            has_binary_trait=false #NNMM does not support binary traits now
            if is_multi_trait2
                mme2.R.val = sample_variance(wArray2, length(mme2.obsID),
                                        mme2.R.df, mme2.R.scale,
                                        invweights2,mme2.R.constraint;
                                        binary_trait_index=has_binary_trait ? findall(x->x=="categorical(binary)", mme2.traits_type) : false)
                Ri    = kron(inv(mme2.R.val),spdiagm(0=>invweights2))
            else #single trait
                if !has_categorical_trait && !has_binary_trait # fixed =1 for single categorical/binary trait
                    mme2.ROld  = mme2.R.val
                    mme2.R.val = sample_variance(ycorr2,length(ycorr2), mme2.R.df, mme2.R.scale, invweights2)
                end
            end
            if mme2.MCMCinfo.double_precision == false
                mme2.R.val = Float32.(mme2.R.val)
            end
        end


        #update σ2_yobs, σ2_weightsNN, and weights_NN
        mme1.σ2_yobs = mme2.R.val
        mme1.σ2_weightsNN = mme2.M[1].G.val
        mme1.weights_NN = mme2.M[1].α[1]

        ########################################################################
        # 3.1 Save MCMC samples
        ########################################################################
        if iter>burnin && (iter-burnin)%output_samples_frequency == 0
            #MCMC samples from posterior distributions
            nsamples       = (iter-burnin)/output_samples_frequency
            output_posterior_mean_variance(mme1,nsamples)
            #mean and variance of posterior distribution
            if mme1.pedTrmVec!=0
                polygenic_effects_variance = inv(mme1.rndTrmVec[polygenic_pos].Gi.val) 
            else
                polygenic_effects_variance=false 
            end
            output_MCMC_samples(mme1,mme1.R.val,polygenic_effects_variance,outfile)
            #  if causal_structure != false
            #      writedlm(causal_structure_outfile,sample4λ_vec',',')
            #  end
        end
        ########################################################################
        # 3.2 Printout
        ########################################################################
        if iter%mme1.MCMCinfo.printout_frequency==0 && iter>burnin
            println("\nPosterior means at iteration: ",iter)
            println("Residual variance: ",round.(mme1.meanVare,digits=6))
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
        for Mi in mme1.M
            mv(output_folder*"/MCMC_samples_marker_effects_variances"*"_"*Mi.name*".txt",
               output_folder*"/MCMC_samples_genetic_variance(REML)"*"_"*Mi.name*".txt")
        end
    end

    output=output_result(mme1,output_folder,
                         mme1.solMean,mme1.meanVare,
                         mme1.pedTrmVec!=0 ? mme1.G0Mean : false,
                         mme1.solMean2,mme1.meanVare2,
                         mme1.pedTrmVec!=0 ? mme1.G0Mean2 : false)
    return output
end
