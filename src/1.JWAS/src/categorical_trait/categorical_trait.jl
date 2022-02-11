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

################################################################################
#Multi-trait: one continuous + one threshold:
#1)Sorensen and Gianola, Likelihood, Bayesian, and MCMC Methods in Quantitative
#Genetics
#2)Wang, Chonglong, et al. Bayesian methods for jointly estimating genomic
#breeding values of one continuous and one threshold trait. Plos one 12.4 (2017)
################################################################################
function MT_categorical_trait_setup!(mme)
    #trait1 is continuouse, trait2 is ordinal
    nInd=length(mme.obsID)
    category_obs  = map(Int,mme.ySparse[(nInd+1):end]) #trait2 is ordinal
    ncategories = length(unique(category_obs))
    ntrait = mme.nModels
    starting_value = mme.sol
    cmean      = mme.X*starting_value[1:size(mme.mmeLhs,1)] #maker effects defaulting to all zeros
    μ = mean(cmean[(nInd+1):end]) #trait2 is ordinal
    σ = mme.R[2,2]
    tmin = μ-10σ
    tmax = μ+10σ

    #starting values for thresholds
    if ncategories > 2
        #t_min=μ-10σ < t1=0 < t2=1 < t3 < ... < t_{#category-1} < t_max=μ+10σ
        thresholds = [[tmin, 0]; collect(1:(tmax-1)/(ncategories-2):tmax)]
    else
        error("ncategories shoule > 2.")
    end
    #starting values for liability
    Mu     = cmean[(nInd+1):end] + (mme.R[1,2]/mme.R[1,1])*(mme.ySparse[1:nInd] - cmean[1:nInd]) # vector of full conditional mean for all inds
    var    = mme.R[2,2]-(mme.R[1,2])^2/mme.R[1,1] #full conditional variance, same to all inds
    for i in 1:length(category_obs) #[ind1,ind2,...]
        whichcategory = category_obs[i]
        mme.ySparse[nInd+i] = rand(truncated(Normal(Mu[i], var), thresholds[whichcategory], thresholds[whichcategory+1]))
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

function MT_categorical_trait_sample_liabilities(mme,ycorr,category_obs,thresholds)
    ########################################################################
    # MT: ONLY for one ordinal trait + one continuous trait
    # no missing data for both ordinal trait and continuous trait
    ########################################################################
    nInd=length(mme.obsID)
    # sample liabilities
    cmean = vec(mme.ySparse) - ycorr #liabilities - residuals
    Mu     = cmean[(nInd+1):end] + (mme.R[1,2]/mme.R[1,1])*(mme.ySparse[1:nInd] - cmean[1:nInd]) # vector of full conditional mean for all inds
    var    = mme.R[2,2]-(mme.R[1,2])^2/mme.R[1,1] #full conditional variance, same to all inds
    for i in 1:length(category_obs) #ind1,ind2,...indn
        whichcategory = category_obs[i]
        mme.ySparse[nInd+i] = rand(truncated(Normal(Mu[i], var), thresholds[whichcategory], thresholds[whichcategory+1]))
    end

    # sample thresholds (same as single trait, except t2=1)
    #thresholds: t_min=μ-10σ < t1=0 < t2=1 < t3 < ... < t_{#category-1} < t_max=μ+10σ
    for i in 4:(length(thresholds)-1) #e.g., t2 between categories 2 and 3
        lowerboundry  = maximum(mme.ySparse[category_obs .== (i-1)])
        upperboundry  = minimum(mme.ySparse[category_obs .== i])
        thresholds[i] = rand(Uniform(lowerboundry,upperboundry))
    end
    ycorr = vec(mme.ySparse) - cmean
    return ycorr
end
