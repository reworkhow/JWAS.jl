#                                |---- ylat1 ----- Mα
#yobs ---f(ylat1,ylat2,ylat3)----|---- ylat2 ----- Mα
#                                |---- ylat3 ----- Mα
#
#nonlinear function: e.g.,
# (1) pig_growth(x1,x2) = sqrt(x1^2 / (x1^2 + x2^2))
# (2) neural network: a1*tan(x1)+a2*tan(x2)

#nonlinear_function: #user-provide function, "tanh"

#ycorr: residual for omics data
function sample_latent_traits(yobs,mme,ycorr,nonlinear_function)
    ylats_old = mme.ySparse         # current values of each latent trait; [trait_1_obs;trait_2_obs;...]
    μ_ylats   = mme.ySparse - ycorr # mean of each latent trait, [trait_1_obs-residuals;trait_2_obs-residuals;...]
                                    # = vcat(getEBV(mme,1).+mme.sol[1],getEBV(mme,2).+mme.sol[2]))
    σ2_yobs      = mme.NNMM.σ2_yobs      # residual variance of yobs (scalar)
    σ2_weightsNN = mme.NNMM.σ2_weightsNN # variance of nn weights between middle and output layers

    incomplete_omics = mme.NNMM.incomplete_omics #indicator for ind with incomplete omics

    #reshape the vector to nind X ntraits
    nobs, ntraits = length(mme.obsID), mme.nModels
    ylats_old     = reshape(ylats_old,nobs,ntraits) #Tianjing's mme.Z
    ylats_old2    = copy(ylats_old)
    μ_ylats       = reshape(μ_ylats,nobs,ntraits)
    ycorr2        = reshape(ycorr,nobs,ntraits)

    #sample latent trait
    if sum(incomplete_omics) != 0   #at least 1 ind with incomplete omics
        #step 1. sample latent trait for individuals with incomplete omics data
        if mme.NNMM.is_activation_fcn == true #Neural Network with activation function
            ylats_new = hmc_one_iteration(10,0.1,ylats_old[incomplete_omics,:],yobs[incomplete_omics],mme.NNMM.weights_NN,mme.R,σ2_yobs,ycorr2[incomplete_omics,:],nonlinear_function)
        else  #user-defined function, MH
            candidates       = μ_ylats[incomplete_omics,:]+randn(size(μ_ylats[incomplete_omics,:]))  #candidate samples
            if nonlinear_function == "Neural Network (MH)" # current NN-MM do not allow "Neural Network (MH)"
                μ_yobs_candidate = [ones(sum(incomplete_omics)) nonlinear_function.(candidates)]*mme.NNMM.weights_NN
                μ_yobs_current   = X*mme.NNMM.weights_NN
            else #user-defined non-linear function
                μ_yobs_candidate = nonlinear_function.(Tuple([view(candidates[incomplete_omics,:],:,i) for i in 1:ntraits])...)
                μ_yobs_current   = nonlinear_function.(Tuple([view(ylats_old[incomplete_omics,:],:,i) for i in 1:ntraits])...)
            end
            llh_current      = -0.5*(yobs - μ_yobs_current ).^2/σ2_yobs
            llh_candidate    = -0.5*(yobs - μ_yobs_candidate).^2/σ2_yobs
            mhRatio          = exp.(llh_candidate - llh_current)
            updateus         = rand(sum(incomplete_omics)) .< mhRatio
            ylats_new        = candidates.*updateus + ylats_old.*(.!updateus)
        end
        #step 2. update ylats with sampled latent traits
        ylats_old[incomplete_omics,:] = ylats_new
        #step 3. for individuals with partial omics data, put back the partial real omics.
        ylats_old[mme.missingPattern] .= ylats_old2[mme.missingPattern]
    end
    #update samples for saving results
    mme.NNMM.middle_nodes = ylats_old

    #sample neural network weights between hidden and output layer
    ############################################################################
    #remove testing individuals
    # testing individuals will be included in Bayesian Alphabet, but they should
    # be removed when we sample weights between hidden and output layer
    ############################################################################
    # the yobs for testing ind is missing, e.g., [0.1,0.2,missing]
    trainingInd      = .!ismissing.(yobs)
    yobs_train       = yobs[trainingInd]
    ylats_old_train  = ylats_old[trainingInd,:]
    nTrain           = sum(trainingInd)
    nOmics           = size(ylats_old_train,2)

    #RR-BLUP to sample neural network weights: omics~yobs (no needed for user-defined function)
    if mme.NNMM.is_activation_fcn == true #Neural Network with activation function

        #calcluate g(Z), e.g., tanh(Z)
        g_z       = nonlinear_function.(ylats_old_train)
        yobs_corr = yobs_train .- mme.NNMM.weights_NN[1] - g_z*mme.NNMM.weights_NN[2:end] #y corrected for all

        #sample mean of phenotype
        # μ1~N( sum(y-g(Z)*w1)/n , σ2_e/n )
        yobs_corr = yobs_corr .+ mme.NNMM.weights_NN[1]
        μ1 = sum(yobs_corr)/nTrain + randn()*sqrt(σ2_yobs/nTrain)
        mme.NNMM.weights_NN[1] = μ1
        yobs_corr = yobs_corr .- mme.NNMM.weights_NN[1]

        #sample omics effects
        # w1~N( lhs^-1 * rhs, lhs^-1 σ2_e  )
        for j in 1:nOmics
            x=g_z[:,j]
            rhs=dot(x,yobs_corr)+dot(x,x)*mme.NNMM.weights_NN[j+1] #the first element is μ1
            lhs=dot(x,x)+σ2_yobs/mme.NNMM.σ2_weightsNN
            invLhs=1/lhs
            meanj=invLhs*rhs
            oldAlpha = mme.NNMM.weights_NN[j+1]
            mme.NNMM.weights_NN[j+1]=meanj+randn()*sqrt(invLhs*σ2_yobs)
            yobs_corr=yobs_corr+(oldAlpha-mme.NNMM.weights_NN[j+1])*x
        end
    end

    #update ylats
    mme.ySparse = vec(ylats_old)
    ycorr[:]    = mme.ySparse - vec(μ_ylats) # =(ylats_new - ylats_old) + ycorr: update residuls (ycorr)

    #sample σ2_yobs
    if mme.NNMM.is_activation_fcn == false  # user-defined nonlinear function
        residuals = yobs_train-nonlinear_function.(Tuple([view(ylats_old_train,:,i) for i in 1:ntraits])...)
    else   # Neural Network with activation function
        residuals = yobs_corr
    end

    if mme.NNMM.fixed_σ2_NN==false
        mme.NNMM.σ2_yobs      = sample_variance(residuals, nTrain, 4, 1) #(dot(x,x) + df*scale)/rand(Chisq(n+df))
        mme.NNMM.σ2_weightsNN = sample_variance(mme.NNMM.weights_NN[2:end], nOmics, 4, 1)
    end
end
