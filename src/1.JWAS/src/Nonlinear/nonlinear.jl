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
    σ2_yobs   = mme.σ2_yobs         # residual variance of yobs (scalar)

    #reshape the vector to nind X ntraits
    nobs, ntraits = length(mme.obsID), mme.nModels
    ylats_old     = reshape(ylats_old,nobs,ntraits) #Tianjing's mme.Z
    ylats_old2    = copy(ylats_old)
    μ_ylats       = reshape(μ_ylats,nobs,ntraits)
    ycorr2        = reshape(ycorr,nobs,ntraits)

    #sample latent trait
    ############################################################################
    #remove individuals with full omics
    # in hmc, individuals are assumed to be independent, so individuals with full
    # omics data can be ignored. If a individual has at least one missing omics,
    # we will sample all omics, then only update the missing omics.
    ############################################################################
    # step1. find individuals with full omics data
    n_real_omis = sum(mme.missingPattern,dims=2) #number of true omics for each ind
    n_omics     = length(mme.lhsVec)             #number of omics
    full_omics = n_real_omis .== n_omics         #indicator for ind with full omics
    no_full_omics = vec(.!full_omics)            #indicator for ind with no/partial omics

    # step2. remove individuals with full omics data, only keep ind with no/partial omics
    ylats_old_tmp = ylats_old[no_full_omics,:]
    yobs_tmp     =yobs[no_full_omics]
    ycorr2_tmp    = ycorr2[no_full_omics,:]

    if mme.is_activation_fcn == true #Neural Network with activation function
        if sum(no_full_omics) != 0   #have at least 1 ind with no/partial omics
            #step 3. sample latent trait for individuals with no/partial omics data
            ylats_new = hmc_one_iteration(10,0.1,ylats_old_tmp,yobs_tmp,mme.weights_NN,mme.R,σ2_yobs,ycorr2_tmp,nonlinear_function)
            #step 4. put the sampled latent trait back to latent trait for all individuals (i.e., update ylats_old)
            ylats_old[no_full_omics,:] = ylats_new
            #step 5. for individuals with partial omics data, put back the partial real omics.
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
    end

    #sample neural network weights between hidden and output layer
    ############################################################################
    #remove testing individuals
    # testing individuals will be included in Bayesian Alphabet, but they should
    # be removed when we sample weights between hidden and output layer
    ############################################################################
    # the yobs for testing ind is missing, e.g., [0.1,0.2,missing]
    trainingInd = .!ismissing.(yobs)
    yobs_train       = yobs[trainingInd]
    ylats_old_train  = ylats_old[trainingInd,:]

    if mme.is_activation_fcn == true #Neural Network with activation function
        if mme.nnweight_lambda == false #fixed effect: 0.000001; known random:e.g.,2.5;unknown random: false
            nnweight_lambda = mme.σ2_yobs/mme.σ2_weightsNN
        else
            nnweight_lambda = mme.nnweight_lambda
        end

        X       = Matrix([ones(length(yobs_train)) nonlinear_function.(ylats_old_train)])
        lhs     = X'X + I*0.00001 + I*nnweight_lambda
        lhs[1,1]= lhs[1,1] - nnweight_lambda

        try
           Ch      = cholesky(lhs)
        catch
           println("--------PosDefException: matrix is not positive definite; Cholesky factorization failed.")
           writedlm("X.txt",X)
           writedlm("nnweight_lambda.txt",nnweight_lambda)
           writedlm("lhs.txt",lhs)

           error("----------PosDefException")
        end
        Ch      = cholesky(lhs)
        L       = Ch.L
        iL      = inv(L)
        rhs     = X'yobs_train
        weights = Ch\rhs + iL'randn(size(X,2))*sqrt(σ2_yobs)
        mme.weights_NN = weights
    end

    #update ylats
    mme.ySparse = vec(ylats_old)
    ycorr[:]    = mme.ySparse - vec(μ_ylats) # =(ylats_new - ylats_old) + ycorr: update residuls (ycorr)

    #sample σ2_yobs
    if mme.is_activation_fcn == false  # user-defined nonlinear function
        residuals = yobs_train-nonlinear_function.(Tuple([view(ylats_old_train,:,i) for i in 1:ntraits])...)
    else   # Neural Network with activation function
        residuals = yobs_train-[ones(length(yobs_train)) nonlinear_function.(ylats_old_train)]*mme.weights_NN
    end
    mme.σ2_yobs= dot(residuals,residuals)/rand(Chisq(length(yobs_train))) #(dot(x,x) + df*scale)/rand(Chisq(n+df))
    mme.σ2_weightsNN = dot(mme.weights_NN,mme.weights_NN)/rand(Chisq(length(mme.weights_NN)))
end
