#                                |---- ylat1 ----- Mα
#yobs ---f(ylat1,ylat2,ylat3)----|---- ylat2 ----- Mα
#                                |---- ylat3 ----- Mα
#
#nonlinear function: e.g.,
# (1) pig_growth(x1,x2) = sqrt(x1^2 / (x1^2 + x2^2))
# (2) neural network: a1*tan(x1)+a2*tan(x2)
function sample_latent_traits(yobs,mme,ycorr,var,nonlinear_function)
    ylats_old = mme.ySparse  # current values of each latent trait
    μ_ylats   = mme.ySparse - ycorr # mean of each latent traits
    #reshape the vector to nind X ntraits
    nobs, ntraits = length(mme.obsID), mme.nModels
    ylats_old     = reshape(ylats_old,nobs,ntraits)
    μ_ylats       = reshape(μ_ylats,nobs,ntraits)
    ylats_new     = copy(ylats_old)

    for i = 1:nobs
        candidates    = randn(mme.nModels)+ μ_ylats[i,:]  #candidate samples
        llh_current   = -0.5*(yobs[i] - nonlinear_function(Tuple(ylats_old[i,:])...))^2/var
        llh_candidate = -0.5*(yobs[i] - nonlinear_function(Tuple(candidates)...))^2/var
        mhRatio = exp(llh_candidate - llh_current)
        if rand() < mhRatio
            ylats_new[i,:] = candidates
        end
    end
    ycorr[:]    = ycorr + vec(ylats_new - ylats_old)
    mme.ySparse = vec(ylats_new)
end
