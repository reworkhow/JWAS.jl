function censored_trait_setup!(mme)
    lower_bound  = Float64.(vec(mme.ySparse))
    upper_bound  = Float64.(mme.MCMCinfo.censored_trait)
    #liabilities are sampled and stored in mme.ySparse
    starting_value = mme.sol
    cmean          = mme.X*starting_value[1:size(mme.mmeLhs,1)] #maker effects defaulting to all zeros
    for i in 1:length(mme.ySparse) #1,2,2,3,1...
        #mme.ySparse[i] = rand(truncated(Normal(cmean[i], sqrt(mme.R)), lower_bound[i], upper_bound[i]))
        mme.ySparse[i] = rand((lower_bound[i]:upper_bound[i]))
    end
    return lower_bound,upper_bound
end

function censored_trait_sample_liabilities(mme,ycorr,lower_bound,upper_bound)
    ########################################################################
    # liabilities
    ########################################################################
    cmean = mme.ySparse - ycorr #liabilities - residuals
    for i in 1:length(mme.ySparse) #1,2,2,3,1...
        mme.ySparse[i] = rand(truncated(Normal(cmean[i], sqrt(mme.R)), lower_bound[i], upper_bound[i]))
    end
    ycorr = mme.ySparse - cmean
    return ycorr
end
