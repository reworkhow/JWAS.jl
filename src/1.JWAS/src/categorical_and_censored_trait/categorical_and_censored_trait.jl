################################################################################
# Categorial and censored traits
#1)Sorensen and Gianola, Likelihood, Bayesian, and MCMC Methods in Quantitative
#Genetics
#2)Wang, Chonglong, et al. Bayesian methods for jointly estimating genomic
#breeding values of one continuous and one threshold trait. Plos one 12.4 (2017)
#3)Wang et al.(2013). Bayesian methods for estimating GEBVs of threshold traits.
#Heredity, 110(3), 213–219.
################################################################################
function categorical_censored_traits_setup!(mme,df)
    ############################################################################
    # initialize: mme.thresholds, lower_bound, upper_bound, liability(=mme.ySparse)
    ############################################################################
    nInd           = length(mme.obsID)
    nTrait         = mme.nModels
    R              = mme.R

    starting_value = mme.sol
    cmean          = mme.X*starting_value[1:size(mme.mmeLhs,1)] #maker effects defaulting to all zeros
    cmean          = reshape(cmean,nInd,nTrait)

    ySparse        = reshape(mme.ySparse,nInd,nTrait) #mme.ySparse will also change since reshape is a reference, not copy

    category_obs   = Dict() #e.g., [1] => [3,2,2,3,1,..] saves the observations for the 1st trait
    upper_bound    = Dict() #e.g., [1] => [Inf,Inf,9,..] saves the upper bound for the 1st trait
    lower_bound    = Dict() #e.g., [1] => [-Inf,-Inf,1,..] saves the lower bound for the 1st trait

    for t in 1:nTrait
        if mme.traits_type[t] ∈ ["categorical","categorical(binary)","censored"]
            if mme.traits_type[t] == "categorical" #categorial traits with >2 categories
                category_obs[t] = map(Int,ySparse[:,t])
                ncategories     = length(unique(category_obs[t]))
                # step1. initialize thresholds
                if nTrait==1 #single categorical trait: [t_min=-Inf < t1=0 < t2 <...< t_{#category-1} < t_max=+Inf], where t_{#category-1}<1
                    mme.thresholds[t] = [-Inf;range(0, length=ncategories,stop=1)[1:(end-1)];Inf]
                else #multiple traits: [t_min=-Inf < t1=0 < t2=1 < t3 <...< t_{#category-1} < t_max=+Inf], where t_{#category-1}<μ+10σ
                    μ = mean(cmean[:,t])
                    σ = R[t,t]
                    mme.thresholds[t] = [-Inf; 0; range(1,length=ncategories-1,stop=μ+10σ)[1:(end-1)];Inf]
                end
                # step2. update lower_bound and upper_bound
                whichcategory  = category_obs[t]
                lower_bound[t] = mme.thresholds[t][whichcategory]
                upper_bound[t] = mme.thresholds[t][whichcategory.+1]
            elseif mme.traits_type[t] == "categorical(binary)"
                category_obs[t] = map(Int,ySparse[:,t])
                # step1. initialize thresholds
                mme.thresholds[t] = [-Inf, 0, Inf]
                # step2. update lower_bound and upper_bound
                whichcategory  = category_obs[t]
                lower_bound[t] = mme.thresholds[t][whichcategory]
                upper_bound[t] = mme.thresholds[t][whichcategory.+1]
            elseif mme.traits_type[t] == "censored"
                # (step1. initialize thresholds is not required for censored trait)
                # step2. update lower_bound and upper_bound
                lower_bound[t]  = df[!,Symbol.(mme.lhsVec[t],"_l")]
                upper_bound[t]  = df[!,Symbol.(mme.lhsVec[t],"_u")]
            end
            #step3. set up liability (=mme.ySparse)
            for i in 1:nInd
                if lower_bound[t][i] != upper_bound[t][i]
                    ySparse[i,t] = rand(truncated(Normal(cmean[i,t], sqrt(R[t,t])), lower_bound[t][i], upper_bound[t][i])) #mme.R has been fixed to 1.0 for category single-trait analysis
                else
                    ySparse[i,t] = lower_bound[t][i]
                end
            end
        end
    end

    ySparse = reshape(ySparse,nInd*nTrait,1)

    return lower_bound,upper_bound,category_obs
end



function categorical_trait_sample_threshold!(mme,lower_bound,upper_bound,category_obs)
    ############################################################################
    # update mme.thresholds for categorical traits with >2 categories
    # lower_bound and upper_bound should also be updated due to the update of mme.thresholds
    # note that threashold for binary trait is fixed as [-Inf,0,Inf]
    ############################################################################
    # single-trait: (only t1 is fixed)
    #   t_min=-Inf < t1=0 < t2 < ... < t_{categories-1} < t_max=Inf
    # multi-trait: (both t1 and t2 are fixed)
    #   t_min=-Inf < t1=0 < t2=1 < t3 < ... < t_{#category-1} < t_max=Inf
    nInd                    = length(mme.obsID)
    nTrait                  = mme.nModels
    start_index             = nTrait==1 ? 3 : 4  #some thresholds were fixed
    ySparse                 = reshape(mme.ySparse, nInd, nTrait)

    for t in 1:nTrait
        if mme.traits_type[t] == "categorical" #categorial traits with >2 categories
            for i in start_index:(length(mme.thresholds[t])-1) #e.g., t2 between categories 2 and 3; will be skipped for binary trait
                lowerboundry_threshold  = maximum(ySparse[:,t][category_obs[t] .== (i-1)]) #lower bound for current threshold
                upperboundry_threshold  = minimum(ySparse[:,t][category_obs[t] .== i])     #upper bound for current threshold
                mme.thresholds[t][i] = rand(Uniform(lowerboundry_threshold,upperboundry_threshold))
            end
            ####################################################################
            # update lower_bound, upper_bound from thresholds
            ####################################################################
            whichcategory  = category_obs[t]
            lower_bound[t] = mme.thresholds[t][whichcategory]
            upper_bound[t] = mme.thresholds[t][whichcategory.+1]
        end
    end
end


function sample_liabilities!(mme,ycorr,lower_bound,upper_bound)
    ############################################################################
    # update mme.ySparse (i.e., liability) for censored and categorical traits
    # ycorr should be updated due to the update of mme.ySparse
    ############################################################################
    cmean          = mme.ySparse - ycorr
    nInd           = length(mme.obsID)
    nTrait         = mme.nModels
    is_multi_trait = nTrait>1
    R              = mme.R
    cmean          = reshape(cmean,      nInd,nTrait)
    ySparse        = reshape(mme.ySparse,nInd,nTrait) #mme.ySparse will also be updated since reshape is a reference, not copy
    ycorr_ref      = reshape(ycorr,      nInd,nTrait) #ycorr will be updated since reshape is a reference, not copy

    if length(findall(t -> t ∈ ["categorical","categorical(binary)","censored"], mme.traits_type)) >1 #at least two traits with liability
        nGibbs = 5     #need Gibbs sampling because liability is not from the stationary distribution
    else  #only one trait with liability, Gibbs sampling is not required
        nGibbs = 1
    end
    for iter in 1:nGibbs
        for t in 1:nTrait
            if mme.traits_type[t] ∈ ["categorical","categorical(binary)","censored"]
                index1 = t  # "1" denotes the trait for sampling liability, "2" denotes all other traits.
                index2 = deleteat!(collect(1:nTrait),index1)
                d      = ycorr_ref[:,index2] #current residuals for all other traits (d)
                #sample residual trait "1"
                μ_1    = is_multi_trait ? vec(R[index1,index2]'inv(R[index2,index2])*d') : zeros(nInd)
                σ2_1   = is_multi_trait ? R[index1,index1]-R[index1,index2]'inv(R[index2,index2])*R[index2,index1] : R #R=1 for single categorical trait
                ϵ1_lower_bound = lower_bound[t] - cmean[:,index1] #thresholds[t][whichcategory] - cmean[:,index1]
                ϵ1_upper_bound = upper_bound[t] - cmean[:,index1] #thresholds[t][whichcategory+1] - cmean[:,index1]
                for i in 1:nInd
                    if ϵ1_lower_bound[i]   != ϵ1_upper_bound[i]
                        ϵ1                  = rand(truncated(Normal(μ_1[i], sqrt(σ2_1)), ϵ1_lower_bound[i], ϵ1_upper_bound[i]))
                        ySparse[i,index1]   = cmean[i,index1] + ϵ1 #mme.ySparse will also be updated since reshape is a reference, not copy
                        ycorr_ref[i,index1] = ϵ1                   #ycorr will also be updated since reshape is a reference, not copy
                    end
                end
            end
        end
    end
end



# print an example for deprecated JWAS function runMCMC(categorical_trait,censored_trait)
function print_single_categorical_censored_trait_example()
    @error("The arguments 'categorical_trait' and  'censored_trait' has been moved to build_model(). Please check our latest example.")
    printstyled("1. Example to build model for single categorical trait:\n"; color=:red)
    printstyled("      model_equation  = \"y = intercept + x1 + x2 + x2*x3 + ID + dam + genotypes\"
      model = build_model(model_equation,categorical_trait=[\"y\"])\n"; color=:red)
    printstyled("2. Example to build model for single censored trait:\n"; color=:red)
    printstyled("      model_equation  = \"y = intercept + x1 + x2 + x2*x3 + ID + dam + genotypes\"
      model = build_model(model_equation,censored_trait=[\"y\"])\n"; color=:red)
end



# # conditional inverse Wishart distribution, where binary traits are independent with unit variance
# # in practice, this may cause convergence problem for residual covariance matrix
# function sample_from_conditional_inverse_Wishart(df,scale,binary_trait_index)
#     #df: degree of freedom
#     #scale: scale
#     #binary_trait_index: index of binary trait, e.g,[1,3] means the 1st and 3rd traits are binary
#     ntraits= size(scale,2)
#     index2 = binary_trait_index         #index for binary traits "2"
#     index1 = setdiff(1:ntraits,index2)  #index for non-binary traits "1"
#     n1     = length(index1)
#     n2     = length(index2)
#
#     V11    = scale[index1,index1]
#     V12    = scale[index1,index2]
#     V22    = scale[index2,index2]
#     V22_1  = V22 - V12'*inv(V11)*V12    #n2 - by - n2
#     X1     = rand(Wishart(df,V11))  # X1~W(V11,n) for non-binary trait, n1-by-n1
#     μ      = vec(inv(V11)*V12)          # n1*n2 - by - 1
#     Σ      = Symmetric(kron(V22_1,inv(X1))) # n1*n2 - by - n1*n2, the order of kron in paper is incorrect, check wikipedia for correct order
#     X2     = rand(MvNormal(μ, Σ))       # n1*n2 - by - 1
#     X2     = reshape(X2,n1,n2)          # n1 - by - n2
#     T11    = inv(X1) + X2*X2'           # n1 - by - n1
#     T12    = -X2                        # n1 - by - n2
#     R      = [T11  T12
#               T12' I(n2)]               #traits are ordered as [index1;index2]
#     re_order = sortperm([index1;index2])
#     R        = R[re_order,re_order]     #traits are ordered as [1,2,3,...]
#     return R
# end
