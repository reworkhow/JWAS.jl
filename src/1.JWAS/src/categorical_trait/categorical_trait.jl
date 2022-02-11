################################################################################
#(single) Threshold trait:
#1)Sorensen and Gianola, Likelihood, Bayesian, and MCMC Methods in Quantitative
#Genetics
#2)Wang et al.(2013). Bayesian methods for estimating GEBVs of threshold traits.
#Heredity, 110(3), 213–219.
################################################################################
function categorical_trait_setup!(mme)
    #starting values for thresholds  -∞ < t1=0 < t2 < ... < t_{#category-1} < +∞
    # where t1=0 (must be fixed to center the distribution) and t_{#category-1}<1.
    #Then liabilities are sampled and stored in mme.ySparse
    category_obs  = map(Int,mme.ySparse) # categories (1,2,2,3,1...)
    ncategories = length(unique(category_obs))
    #-Inf,t1,t2,...,t_{#c-1},Inf, where t1=0 and t_{#c-1}<1
    thresholds = [-Inf;range(0, length=ncategories,stop=1)[1:(end-1)];Inf]

    starting_value = mme.sol
    cmean      = mme.X*starting_value[1:size(mme.mmeLhs,1)] #maker effects defaulting to all zeros
    for i in 1:length(category_obs) #1,2,2,3,1...
        whichcategory = category_obs[i]
        mme.ySparse[i] = rand(truncated(Normal(cmean[i], 1), thresholds[whichcategory], thresholds[whichcategory+1]))
    end
    return category_obs,thresholds
end



function categorical_trait_sample_liabilities(mme,ycorr,category_obs,thresholds)
    ########################################################################
    # liabilities
    ########################################################################
    cmean = vec(mme.ySparse) - ycorr #liabilities - residuals
    for i in 1:length(category_obs) #1,2,2,3,1...
        whichcategory = category_obs[i]
        mme.ySparse[i] = rand(truncated(Normal(cmean[i], 1), thresholds[whichcategory], thresholds[whichcategory+1]))
    end
    ########################################################################
    # Thresholds
    ########################################################################
    #thresholds -∞,t1,t2,...,t_{categories-1},+∞
    #the threshold t1=0 must be fixed to center the distribution
    for i in 3:(length(thresholds)-1) #e.g., t2 between categories 2 and 3
        lowerboundry  = maximum(mme.ySparse[category_obs .== (i-1)])
        upperboundry  = minimum(mme.ySparse[category_obs .== i])
        thresholds[i] = rand(Uniform(lowerboundry,upperboundry))
    end
    ycorr = vec(mme.ySparse) - cmean
    return ycorr
end
