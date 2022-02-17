################################################################################
# Ordinal trait
#1)Sorensen and Gianola, Likelihood, Bayesian, and MCMC Methods in Quantitative
#Genetics
#2)Wang, Chonglong, et al. Bayesian methods for jointly estimating genomic
#breeding values of one continuous and one threshold trait. Plos one 12.4 (2017)
#3)Wang et al.(2013). Bayesian methods for estimating GEBVs of threshold traits.
#Heredity, 110(3), 213–219.
################################################################################
function categorical_trait_setup!(mme,lower_bound,upper_bound)
    nInd                    = length(mme.obsID)
    nTrait                  = mme.nModels
    categorical_trait_index = mme.categorical_trait_index #[2,4] means the 2nd and 4th traits are ordinal
    n_categorical_trait     = length(categorical_trait_index)
    is_multi_trait          = nTrait>1
    R                       = mme.R

    ySparse      = reshape(mme.ySparse,nInd,nTrait)

    category_obs = map(Int,ySparse[:,categorical_trait_index])
    ncategories  = [length(unique(col)) for col in eachcol(category_obs)] #[4,5] means 4 and 5 categories for two ordinal traits, respectively

    starting_value = mme.sol
    cmean = mme.X*starting_value[1:size(mme.mmeLhs,1)] #maker effects defaulting to all zeros
    cmean = reshape(cmean,nInd,nTrait)
    ############################################################################
    # thresholds
    ############################################################################
    thresholds_all = Array{Array{Float64},1}(undef,n_categorical_trait) #each element is thresholds for one trait
    if is_multi_trait
        if all(ncategories .> 2) #all traits have 3 categories
            μ = vec(mean(cmean[:,categorical_trait_index],dims=1)) # #categoricalTrait-by-1
            σ = diag(R)[categorical_trait_index]                   # #categoricalTrait-by-1
            tmin = μ-10σ                                           # #categoricalTrait-by-1
            tmax = μ+10σ                                           # #categoricalTrait-by-1
            for t in 1:n_categorical_trait #t_min=μ-10σ < t1=0 < t2=1 < t3 < ... < t_{#category-1} < t_max=μ+10σ
                thresholds_all[t] = [[tmin[t], 0]; collect(1:(tmax[t]-1)/(ncategories[t]-2):tmax[t])]
                # thresholds_all[t] = [-Inf, 0, 1, Inf]
            end
        else
            error("ncategories shoule > 2.")
        end
    else #single ordinal trait: -Inf,t1,t2,...,t_{#c-1},Inf, where t1=0 and t_{#c-1}<1
        thresholds_all[1] = [-Inf;range(0, length=ncategories[1],stop=1)[1:(end-1)];Inf]
    end
    # update lower bound and upper bound from thresholds
    lower_bound,upper_bound = update_bounds_from_threshold(lower_bound,upper_bound,category_obs,thresholds_all,mme.categorical_trait_index)

    ############################################################################
    # liability
    ############################################################################
    for t in 1:nTrait
        if mme.lhsTag[t] == "categorical"
            for i in 1:nInd
                ySparse[i,t] = rand(truncated(Normal(cmean[i,t], is_multi_trait ? R[t,t] : 1), lower_bound[i,t],upper_bound[i,t]))
            end
        end
    end

    return category_obs,thresholds_all,lower_bound,upper_bound #category_obs: nInt-by-#categoricaTrait;
end                                    #thresholds_all: vector of length #categoricalTrait, each element is its threshold


function categorical_trait_sample_threshold(mme, thresholds_all, category_obs)
    ########################################################################
    # Thresholds
    ########################################################################
    # single categorical trait:
    #   thresholds: t_min=-∞ < t1=0 < t2 < ... < t_{categories-1} < t_max=+∞
    # multiple traits:
    #   thresholds: t_min=μ-10σ < t1=0 < t2=1 < t3 < ... < t_{#category-1} < t_max=μ+10σ
    nInd                    = length(mme.obsID)
    nTrait                  = mme.nModels
    categorical_trait_index = mme.categorical_trait_index
    n_categorical_trait     = length(categorical_trait_index)
    is_multi_trait          = nTrait>1
    start_index             = is_multi_trait ? 4 : 3  #single trait: sample t_2,...,t_{#category-1}; multi-trait: sample t_3,...,t_{#category-1}
    ySparse                 = reshape(mme.ySparse, nInd, nTrait)

    for t in 1:n_categorical_trait
        thresholds = thresholds_all[t]
        categorical_trait = categorical_trait_index[t]
        for i in start_index:(length(thresholds)-1) #e.g., t2 between categories 2 and 3
            lowerboundry  = maximum(ySparse[:,categorical_trait][category_obs[:,t].==(i-1)])
            upperboundry  = minimum(ySparse[:,categorical_trait][category_obs[:,t].== i])
            thresholds[i] = rand(Uniform(lowerboundry,upperboundry))
        end
        thresholds_all[t]=thresholds
    end
    return thresholds_all
end


function sample_liabilities!(mme,ycorr,lower_bound,upper_bound;nGibbs=50)
    cmean = mme.ySparse - ycorr #liabilities - residuals
    # lower_bound:  # nInd-by-nTrait
    # upper_bound:  # nInd-by-nTrait
    nInd           = length(mme.obsID)
    nTrait         = mme.nModels
    is_multi_trait = nTrait>1
    R              = mme.R
    cmean          = reshape(cmean              ,nInd,nTrait)
    ySparse        = reshape(Matrix(mme.ySparse),nInd,nTrait)
    ycorr          = reshape(ycorr              ,nInd,nTrait)

    for iter in 1:nGibbs
        for t in 1:nTrait
            if mme.lhsTag[t] ∈ ["categorical","censored"]
                index1 = t  # "1" denotes the trait for sampling liability, "2" denotes all other traits.
                index2 = deleteat!(collect(1:nTrait),index1)
                d      = ySparse[:,index2]-cmean[:,index2] #current residuals for all other traits (d)
                #sample residual trait "1"
                μ_1    = is_multi_trait ? vec(R[index1,index2]'inv(R[index2,index2])*d') : zeros(nInd)
                σ2_1   = is_multi_trait ? R[index1,index1]-R[index1,index2]'inv(R[index2,index2])*R[index2,index1] : R #R=1 for single categorical trait
                ϵ1_lower_bound = lower_bound[:,t] - cmean[:,index1] #thresholds[t][whichcategory] - cmean[:,index1]
                ϵ1_upper_bound = upper_bound[:,t] - cmean[:,index1] #thresholds[t][whichcategory+1] - cmean[:,index1]
                for i in 1:nInd
                    if ϵ1_lower_bound[i] != ϵ1_upper_bound[i]
                        ϵ1                = rand(truncated(Normal(μ_1[i], sqrt(σ2_1)), ϵ1_lower_bound[i], ϵ1_upper_bound[i]))
                        ySparse[i,index1] = cmean[i,index1] + ϵ1
                        ycorr[i,index1]   = ϵ1
                    end
                end
            end
        end
    end
    mme.ySparse = reshape(ySparse,nInd*nTrait,1)
    ycorr=vec(ycorr)
    # ycorr       = mme.ySparse - vec(cmean)
end


#Below function is to update upper and lower bounds due to the update of thresholds
function update_bounds_from_threshold(lower_bound,upper_bound,category_obs,thresholds,categorical_trait_index)
    # category_obs: nInd-by-#categorical_trait
    # thresholds: #categorical_trait-by-1, each element is a threshold vector
    # lower_bound, upper_bound: nInd-by-nTrait
    for t in 1:length(categorical_trait_index) # number of categorical_trait
        index_t       = categorical_trait_index[t]
        whichcategory = category_obs[:,t]
        lower_bound[:,index_t] = thresholds[t][whichcategory]
        upper_bound[:,index_t] = thresholds[t][whichcategory.+1]
    end

    return lower_bound,upper_bound
end
