#                                |---- ylat1 ----- Mα
#yobs ---f(ylat1,ylat2,ylat3)----|---- ylat2 ----- Mα
#                                |---- ylat3 ----- Mα
#
#nonlinear function: e.g.,
# (1) pig_growth(x1,x2) = sqrt(x1^2 / (x1^2 + x2^2))
# (2) neural network: a1*tan(x1)+a2*tan(x2)
function sample_latent_traits(yobs,mme,ycorr,nonlinear_function)
    ylats_old = mme.ySparse         # current values of each latent trait
    μ_ylats   = mme.ySparse - ycorr # mean of each latent trait
                                    # = vcat(getEBV(mme,1).+mme.sol[1],getEBV(mme,2).+mme.sol[2]))
    σ2_yobs   = mme.σ2_yobs         # residual variance of yobs (scalar)

    #reshape the vector to nind X ntraits
    nobs, ntraits = length(mme.obsID), mme.nModels
    ylats_old     = reshape(ylats_old,nobs,ntraits)
    μ_ylats       = reshape(μ_ylats,nobs,ntraits)

    if nonlinear_function == "Neural Network" #sample weights
        X       = [ones(nobs) tanh.(ylats_old)]
        lhs     = X'X + I*0.00001
        Ch      = cholesky(lhs)
        L       = Ch.L
        iL      = inv(L)
        rhs     = X'yobs
        weights = Ch\rhs + iL'randn(size(X,2))*sqrt(σ2_yobs)
        mme.weights_NN = weights
    end

    candidates       = μ_ylats+randn(size(μ_ylats))  #candidate samples
    if nonlinear_function == "Neural Network"
        μ_yobs_candidate = [ones(nobs) tanh.(candidates)]*weights
        μ_yobs_current   = X*weights
    else
        μ_yobs_candidate = nonlinear_function.(Tuple([view(candidates,:,i) for i in 1:ntraits])...)
        μ_yobs_current   = nonlinear_function.(Tuple([view(ylats_old,:,i) for i in 1:ntraits])...)
    end
    llh_current      = -0.5*(yobs - μ_yobs_current ).^2/σ2_yobs
    llh_candidate    = -0.5*(yobs - μ_yobs_candidate).^2/σ2_yobs
    mhRatio          = exp.(llh_candidate - llh_current)
    updateus         = rand(nobs) .< mhRatio
    ylats_new        = candidates.*updateus + ylats_old.*(.!updateus)

    mme.ySparse = vec(ylats_new)
    ycorr[:]    = mme.ySparse - vec(μ_ylats)

    #sample σ2_yobs
    if nonlinear_function != "Neural Network"
        residuals = yobs-nonlinear_function.(Tuple([view(ylats_new,:,i) for i in 1:ntraits])...)
    else
        residuals = yobs-[ones(nobs) tanh.(ylats_new)]*weights
    end
    mme.σ2_yobs= dot(residuals,residuals)/rand(Chisq(nobs)) #(dot(x,x) + df*scale)/rand(Chisq(n+df))
end
