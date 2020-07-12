#                                |---- ylat1 ----- Mα
#yobs ---f(ylat1,ylat2,ylat3)----|---- ylat2 ----- Mα
#                                |---- ylat3 ----- Mα
#
#nonlinear function: e.g.,
# (1) pig_growth(x1,x2) = sqrt(x1^2 / (x1^2 + x2^2))
# (2) neural network: a1*tan(x1)+a2*tan(x2)
function sample_latent_traits(yobs,mme,ycorr,var,nonlinear_function)
    ylats_old = mme.ySparse         # current values of each latent trait
    μ_ylats   = mme.ySparse - ycorr # mean of each latent trait
    σ2_yobs   = mme.σ2_yobs         # residual variance of yobs (scalar)

    #reshape the vector to nind X ntraits
    nobs, ntraits = length(mme.obsID), mme.nModels
    ylats_old     = reshape(ylats_old,nobs,ntraits)
    μ_ylats       = reshape(μ_ylats,nobs,ntraits)
    ylats_new     = copy(ylats_old)


    if nonlinear_function == "Neural Network" #need an intercept??
        #sample weights
        X       = [ones(nobs) tanh.(ylats_old)]
        lhs     = X'X + I*0.00001
        Ch      = cholesky(lhs)
        L       = Ch.L
        iL      = inv(L)
        rhs     = X'yobs
        weights = Ch\rhs + iL'randn(size(X,2))*sqrt(σ2_yobs)
        #weights = rand(MvNormal(lhs\rhs,inv(lhs)*σ2_yobs))

        μ_yobs_current_all    = X*weights
        candidates_all        = μ_ylats+randn(size(μ_ylats))
        μ_yobs_candidate_all  = [ones(nobs) tanh.(candidates_all)]*weights
    end

    for i = 1:nobs
        if nonlinear_function != "Neural Network"
            candidates       = randn(mme.nModels)+ μ_ylats[i,:]  #candidate samples
            μ_yobs_current   = nonlinear_function(Tuple(ylats_old[i,:])...)
            μ_yobs_candidate = nonlinear_function(Tuple(candidates)...)
        else
            μ_yobs_current   = μ_yobs_current_all[i]
            μ_yobs_candidate = μ_yobs_candidate_all[i]
        end
        llh_current      = -0.5*(yobs[i] - μ_yobs_current )^2/σ2_yobs
        llh_candidate    = -0.5*(yobs[i] - μ_yobs_candidate)^2/σ2_yobs
        mhRatio = exp(llh_candidate - llh_current)
        if rand() < mhRatio
            ylats_new[i,:] = candidates_all[i,:]
        end
    end

    ycorr[:]    = ycorr + vec(ylats_new - ylats_old)
    mme.ySparse = vec(ylats_new)

    #sample σ2_yobs
    if nonlinear_function != "Neural Network"
        sse=0
        for i = 1:nobs
            sse += (yobs[i]-nonlinear_function(Tuple(ylats_new[i,:])...))^2
        end
        mme.σ2_yobs= sse/rand(Chisq(nobs)) #(dot(x,x) + df*scale)/rand(Chisq(n+df))
    else
        yobs_corr = yobs-[ones(nobs) ylats_new]*weights
        mme.σ2_yobs= dot(yobs_corr,yobs_corr)/rand(Chisq(nobs))
    end
end
