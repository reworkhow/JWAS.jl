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
    σ2_yobs      = mme.σ2_yobs      # residual variance of yobs (scalar)
    σ2_weightsNN = mme.σ2_weightsNN # variance of nn weights between middle and output layers

    #reshape the vector to nind X ntraits
    nobs, ntraits = length(mme.obsID), mme.nModels
    ylats_old     = reshape(ylats_old,nobs,ntraits) #Tianjing's mme.Z
    ylats_old2    = copy(ylats_old)
    μ_ylats       = reshape(μ_ylats,nobs,ntraits)
    ycorr2        = reshape(ycorr,nobs,ntraits)

    #sample latent trait
    incomplete_omics = mme.incomplete_omics #indicator for ind with incomplete omics
    if mme.is_activation_fcn == true #Neural Network with activation function
        if sum(incomplete_omics) != 0   #at least 1 ind with incomplete omics
            #step 1. sample latent trait (only for individuals with incomplete omics data
            ylats_new = hmc_one_iteration(10,0.1,ylats_old[incomplete_omics,:],yobs[incomplete_omics],mme.weights_NN,mme.R.val,σ2_yobs,ycorr2[incomplete_omics,:],nonlinear_function)
            #step 2. update ylats with sampled latent traits
            ylats_old[incomplete_omics,:] = ylats_new
            #step 3. for individuals with partial omics data, put back the partial real omics.
            ylats_old[mme.missingPattern] .= ylats_old2[mme.missingPattern]
        end
    else  #user-defined function, MH
        candidates       = μ_ylats+randn(size(μ_ylats))  #candidate samples
        if nonlinear_function == "Neural Network (MH)"
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

        #step 2. update ylats with sampled latent traits
        ylats_old[incomplete_omics,:] = ylats_new
        #step 3. for individuals with partial omics data, put back the partial real omics.
        ylats_old[mme.missingPattern] .= ylats_old2[mme.missingPattern]
    end



    #update ylats
    mme.ySparse = vec(ylats_old)
    ycorr[:]    = mme.ySparse - vec(μ_ylats) # =(ylats_new - ylats_old) + ycorr: update residuls (ycorr)

    #sample σ2_yobs
    if mme.is_activation_fcn == false  # user-defined nonlinear function
        residuals = yobs_train-nonlinear_function.(Tuple([view(ylats_old_train,:,i) for i in 1:ntraits])...)
    else   # Neural Network with activation function
        residuals = yobs_corr
    end

    if mme.fixed_σ2_NN==false
        mme.σ2_yobs      = sample_variance(residuals, nTrain, 4, 1) #(dot(x,x) + df*scale)/rand(Chisq(n+df))
        mme.σ2_weightsNN = sample_variance(mme.weights_NN[2:end], nOmics, 4, 1)
    end
end
