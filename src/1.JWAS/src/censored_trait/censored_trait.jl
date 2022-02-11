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

#### sample liabilities
function MT_censored_trait_sample_liabilities(mme,ycorr,lower_bound,upper_bound)
    ########################################################################
    # MT: ONLY for one censored trait + one continuous trait
    # continouse traits: UB=LB=y
    ########################################################################
    cmean = mme.ySparse - ycorr #liabilities - residuals
    nInd=length(mme.obsID)
    nTrait=mme.nModels
    cmean = reshape(cmean,nInd,nTrait)
    lower_bound = reshape(lower_bound,nInd,nTrait)
    upper_bound = reshape(upper_bound,nInd,nTrait)
    ySparse = reshape(mme.ySparse,nInd,nTrait)
    R=mme.R

    for i in 1:nInd
        if lower_bound[i,:] == upper_bound[i,:] #two traits are continuous; trait1_lb=trait1_ub && trait2_lb=trait2_ub
            ySparse[i,:] = lower_bound[i,:]
        else #one trait is censored: sample censored traits | continuouse trait
            for t in 1:nTrait
                if lower_bound[i,t] == upper_bound[i,t] # t-th trait is continuous
                    ySparse[i,t] = lower_bound[i,t]
                else # t-th trait is censored, use "1" denote censored trait, "2" denote continuous trait
                     # residual for continuous trait = d
                    index1=[t]  #index for censored trait
                    index2=deleteat!(collect(1:2),t) #index for continuous trait
                    d=ySparse[i,index2]-cmean[i,index2]
                    #sample residual from its conditional distribution
                    μ_1 = R[index1,index2]*inv(R[index2,index2])*d
                    σ2_1= R[index1,index1]-R[index1,index2]*inv(R[index2,index2])*R[index2,index1]
                    ϵ1_lower_bound=lower_bound[i,t]-cmean[i,t]
                    ϵ1_upper_bound=upper_bound[i,t]-cmean[i,t]
                    μ_1=μ_1[1]   #a vector of length 1
                    σ2_1=σ2_1[1] #a vector of length 1
                    ϵ1 = rand(truncated(Normal(μ_1,sqrt(σ2_1)), ϵ1_lower_bound, ϵ1_upper_bound))
                    #update ySparse
                    ySparse[i,t] = cmean[i,t] + ϵ1
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
