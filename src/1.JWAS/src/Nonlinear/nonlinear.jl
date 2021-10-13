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
    μ_ylats       = reshape(μ_ylats,nobs,ntraits)
    ycorr2        = reshape(ycorr,nobs,ntraits)

    #testing pattern
    # for NN-Bayes Omics: the yobs for testing ind is missing, e.g., [0.1,0.2,missing]
    # need to remove testing ind from yobs, ylats_old, μ_ylats
    trainingInd = .!ismissing.(yobs)
    yobs2       = yobs[trainingInd]  #(90,)
    ylats_old2  = ylats_old[trainingInd,:]  #(90,5)
    ycorr2      = ycorr2[trainingInd,:] #(90,5)
    missingPattern = mme.missingPattern[trainingInd,:] #(90,5)

    if mme.full_omics==false #skip if we have full omics data
        if mme.is_activation_fcn == true #Neural Network with activation function
            ylats_new = hmc_one_iteration(10,0.1,ylats_old2,yobs2,mme.weights_NN,mme.R,σ2_yobs,ycorr2,nonlinear_function) #(90,5)
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

        if mme.latent_traits !== false  #gene1, gene2, ..
            ylats_new[missingPattern] .= ylats_old2[missingPattern]
        end
    else #mme.full_omics==true
        ylats_new = ylats_old2
    end

    #sample weights
    if mme.is_activation_fcn == true #Neural Network with activation function
        if mme.nnweight_lambda == false #flat prior
            X       = Matrix([ones(length(yobs2)) nonlinear_function.(ylats_new)])
            lhs     = X'X + I*0.00001
            Ch      = cholesky(lhs)
            L       = Ch.L
            iL      = inv(L)
            rhs     = X'yobs2
            weights = Ch\rhs + iL'randn(size(X,2))*sqrt(σ2_yobs)
            mme.weights_NN = weights
        else  #normal prior with fixed variance, e.g., mme.nnweight_lambda == 3002
            X       = nonlinear_function.(ylats_new)  #the weight for last col is 1
            lhs     = [length(yobs2)            sum(X,dims=1)
                       sum(X,dims=1)'  X'X + I*mme.nnweight_lambda]
            rhs     = [sum(yobs2); X'yobs2]
            weights = inv(lhs)*rhs
            mme.weights_NN = weights
        end
    end

    #update ylats
    ylats_old[trainingInd,:]=ylats_new
    mme.ySparse = vec(ylats_old)
    ycorr[:]    = mme.ySparse - vec(μ_ylats) # =(ylats_new - ylats_old) + ycorr: update residuls (ycorr)

    #sample σ2_yobs
    if mme.is_activation_fcn == false  # user-defined nonlinear function
        residuals = yobs2-nonlinear_function.(Tuple([view(ylats_new,:,i) for i in 1:ntraits])...)
    else   # Neural Network with activation function
        residuals = yobs2-[ones(length(yobs2)) nonlinear_function.(ylats_new)]*weights
    end
    mme.σ2_yobs= dot(residuals,residuals)/rand(Chisq(length(yobs2))) #(dot(x,x) + df*scale)/rand(Chisq(n+df))
end
