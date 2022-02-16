function censored_trait_setup!(mme)
    lower_bound  = Float64.(vec(mme.ySparse))                 #[lower_bound_trait1; lower_bound_trait2; ...]
    upper_bound  = vec(Float64.(mme.MCMCinfo.censored_trait)) #[upper_bound_trait1; upper_bound_trait2; ...]
    #liabilities are sampled and stored in mme.ySparse
    starting_value = mme.sol
    cmean          = mme.X*starting_value[1:size(mme.mmeLhs,1)] #maker effects defaulting to all zeros
    for i in 1:length(mme.ySparse) #1,2,2,3,1...
        #mme.ySparse[i] = rand(truncated(Normal(cmean[i], sqrt(mme.R)), lower_bound[i], upper_bound[i]))
        if lower_bound[i] == -Inf
            mme.ySparse[i] = upper_bound[i]
        elseif upper_bound[i] == Inf
            mme.ySparse[i] = lower_bound[i]
        else
            mme.ySparse[i] = rand(lower_bound[i]:upper_bound[i])
        end
    end
    return lower_bound,upper_bound
end


function MT_censored_trait_sample_liabilities(mme,ycorr,lower_bound,upper_bound)
    ########################################################################
    # MT: multiple censored traits and continuous traits
    # continuous traits: upper_bound=lower_bound=y
    ########################################################################
    cmean       = mme.ySparse - ycorr     #liabilities - residuals
    nInd        = length(mme.obsID)
    nTrait      = mme.nModels
    cmean       = reshape(cmean,       nInd,nTrait)
    lower_bound = reshape(lower_bound, nInd,nTrait)
    upper_bound = reshape(upper_bound, nInd,nTrait)
    ySparse     = reshape(mme.ySparse, nInd,nTrait)
    R           = mme.R #residual variance

    #### sample liabilities
    for i in 1:nInd
        for iter=1:5 #Gibbs sampler to sample liabilities for one individual
            for t in 1:nTrait
                if lower_bound[i,t] == upper_bound[i,t] # t-th trait is continuous
                    ySparse[i,t] = lower_bound[i,t]
                else # t-th trait is censored: "1" denotes a censored trait, "2" denotes all other traits.
                    index1 = t                                 #index for the censored trait
                    index2 = deleteat!(collect(1:nTrait),t)    #index for all other traits
                    d      = ySparse[i,index2]-cmean[i,index2] #current residuals for all other traits (d)
                    # sample residual for censored trait "1"
                    μ_1    = R[index1,index2]'inv(R[index2,index2])*d
                    σ2_1   = R[index1,index1]-R[index1,index2]'inv(R[index2,index2])*R[index2,index1] #will be pre-calculated for speedup
                    ϵ1_lower_bound = lower_bound[i,t] - cmean[i,t]
                    ϵ1_upper_bound = upper_bound[i,t] - cmean[i,t]
                    ϵ1             = rand(truncated(Normal(μ_1,sqrt(σ2_1)), ϵ1_lower_bound, ϵ1_upper_bound))
                    ySparse[i,t]   = cmean[i,t] + ϵ1                     #update ySparse
                end
            end
        end
    end
    mme.ySparse = vec(ySparse)
    ycorr = mme.ySparse - vec(cmean)
    return ycorr  #vector of ycorr, ordered by trait
end


function censored_trait_sample_liabilities(mme,ycorr,lower_bound,upper_bound)
    ########################################################################
    # liabilities
    ########################################################################
    cmean = mme.ySparse - ycorr #liabilities - residuals
    for i in 1:length(mme.ySparse) #1,2,2,3,1...
        if lower_bound[i] != upper_bound[i]
            mme.ySparse[i] = rand(truncated(Normal(cmean[i], sqrt(mme.R)), lower_bound[i], upper_bound[i]))
        else
            mme.ySparse[i] = lower_bound[i]
        end
    end
    ycorr = mme.ySparse - cmean
    return ycorr
end
