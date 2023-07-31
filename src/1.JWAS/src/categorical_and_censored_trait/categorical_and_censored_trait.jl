################################################################################
# Categorial and censored traits
#1)Sorensen and Gianola, Likelihood, Bayesian, and MCMC Methods in Quantitative
#Genetics
#2)Wang et al.(2013). Bayesian methods for estimating GEBVs of threshold traits.
#Heredity, 110(3), 213–219.
#3)Korsgaard et al,(2003). Multivariate bayesian analysis of gaussian, right censored gaussian,
#ordered categorical and binary traits using gibbs sampling.
#4)Sorensen et al,(1995). Bayesian inference in threshold models using gibbs sampling.
#5)Sorensen et al,(1998). Bayesian mixed-effects model analysis of a censored normal
# distribution with animal breeding applications.
#6)Korsgaard et al,(1999). A useful reparameterisation to obtain samples from conditional
#inverse wishart distributions.
################################################################################

# |---------------|-------------|-------------------------------------------------------|
# | Analysis Type | Trait Type  |  Parameterization                                     |
# |---------------|-------------|-------------------------------------------------------|
# | single-trait  | censored    | NONE                                                  |
# |               | categorical | vare=1; -Inf < t1=0 < t2 <...< t_{#category-1} < +Inf |
# |               | binary      | vare=1; -Inf < t1=0 < +Inf                            |
# |---------------|-------------|-------------------------------------------------------|
# | multi-trait   | censored    | NONE                                                  |
# |               | categorical | -Inf < t1=0 < t2=1 < t3 <...< t_{#category-1} < +Inf  |
# |               | binary      | vare of binary traits=I; -Inf < t1=0 < +Inf           |
# |---------------|-------------|-------------------------------------------------------|


function categorical_censored_traits_setup!(mme,df)
    ############################################################################
    # Goal:
    #   initialize mme.thresholds, lower_bound, upper_bound, liability(=mme.ySparse)
    ############################################################################
    nInd           = length(mme.obsID)
    nTrait         = mme.nModels
    R              = mme.R.val

    starting_value = mme.sol
    cmean          = mme.X*starting_value[1:size(mme.mmeLhs,1)] #maker effects defaulting to all zeros
    cmean          = reshape(cmean,nInd,nTrait)

    ySparse        = reshape(mme.ySparse,nInd,nTrait) #mme.ySparse will also change since reshape is a reference, not copy

    category_obs   = Dict() #e.g., [1] => [3,1,2,...]         all observations for 1st trait
    upper_bound    = Dict() #e.g., [1] => [Inf,Inf,9.8,...]   upper bounds for all observations for 1st trait
    lower_bound    = Dict() #e.g., [1] => [-Inf,-Inf,1.5,...] lower bounds for all observations for 1st trait

    for t in 1:nTrait
        if mme.traits_type[t] ∈ ["categorical","categorical(binary)","censored"]
            if mme.traits_type[t] ∈ ["categorical","categorical(binary)"]
                ################################################################
                # step1. initialize thresholds for categorical and binary traits
                ################################################################
                category_obs[t] = collect(map(Int,ySparse[:,t]))  #may have 0 for missing categorical trait
                if mme.traits_type[t] == "categorical" #categorial traits with >2 categories
                    ncategories     = length(filter!(x->x!=0,unique(category_obs[t]))) #revome 0 for missing categorical trait
                    if nTrait==1 #single categorical trait: [t_min=-Inf < t1=0 < t2 <...< t_{#category-1} < t_max=+Inf], where t_{#category-1}<1
                        mme.thresholds[t] = [-Inf;range(0, length=ncategories,stop=1)[1:(end-1)];Inf]
                    else #multiple traits: [t_min=-Inf < t1=0 < t2=1 < t3 <...< t_{#category-1} < t_max=+Inf], where t_{#category-1}<μ+10σ
                        μ = mean(cmean[:,t])
                        σ = R[t,t]
                        mme.thresholds[t] = [-Inf; 0; range(1,length=ncategories-1,stop=μ+10σ)[1:(end-1)];Inf]
                    end
                else # binary trait
                    mme.thresholds[t] = [-Inf, 0, Inf]
                end
                #################################################################################
                # step2. initialize lower_bound and upper_bound for categorical and binary traits
                #################################################################################
                lower_bound[t],upper_bound[t]=update_lower_upper_bound_with_threshold(mme.thresholds[t],category_obs[t])
            end
            ###################################################################
            # step3. initialize lower_bound and upper_bound for censored traits
            ###################################################################
            if mme.traits_type[t] == "censored"
                lower_bound[t]  = df[!,Symbol.(mme.lhsVec[t],"_l")]
                upper_bound[t]  = df[!,Symbol.(mme.lhsVec[t],"_u")]
            end
            ##################################################################################
            #step4. initialize liability (=mme.ySparse) for categorical,binary,censored traits
            ##################################################################################
            for i in 1:nInd
                if lower_bound[t][i] != upper_bound[t][i]
                    ySparse[i,t] = rand(truncated(Normal(cmean[i,t], sqrt(R[t,t])), lower_bound[t][i], upper_bound[t][i])) #mme.R.val has been fixed to 1.0 for category single-trait analysis
                else
                    ySparse[i,t] = lower_bound[t][i] #mme.ySparse will also change since reshape is a reference, not copy
                end
            end
        end
    end

    ySparse = reshape(ySparse,nInd*nTrait,1)

    return lower_bound,upper_bound,category_obs
end


function update_lower_upper_bound_with_threshold(thresholds, whichcategory)
    ############################################################################
    # Goal:
    #  for a categorical trait (i.e., whichcategory), give the thresholds, return
    #  lower and upper bound of all individuals.
    # Arguments:
    #  - thresholds: array, the threasholds for current trait. e.g, [-Inf,0,Inf]
    #  - whichcategory: array, the observed categories for current trait. e.g., [1,0,1,0,...]
    # Output:
    #  - lower_bound, upper_bound: vector, length = number of individuals
    # Notes:
    #  - category "0" represents the missing categorical trait
    ############################################################################
    nInd        = length(whichcategory)
    lower_bound = zeros(nInd) # initialize lower bound
    upper_bound = zeros(nInd) # initialize upper bound
    for j in 1:nInd # for each individual
        if whichcategory[j]==0  # "0" represents a missing categorical trait
            lower_bound[j]=-Inf # liability should be sampled from untruncated distribution
            upper_bound[j]=Inf  # thus lower and upper bounds are infinity
        else
            lower_bound[j]=thresholds[whichcategory[j]]
            upper_bound[j]=thresholds[whichcategory[j]+1]
        end
    end
    return lower_bound,upper_bound
end


function categorical_trait_sample_threshold!(mme,lower_bound,upper_bound,category_obs)
    ############################################################################
    # Goal:
    #   update mme.thresholds for categorical traits with >2 categories. Then
    #   lower and upper bounds should also be updated due to the update of mme.thresholds.
    # Arguments:
    #   - lower_bound, upper_bound, category_obs: dictionary, key is the index of the trait
    # Notes:
    #   - threashold for binary trait is fixed as [-Inf,0,Inf]
    #   - "0" in category_obs indicate missing categorical trait
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
            ####################################################################
            # step1. update thresholds
            ####################################################################
            for i in start_index:(length(mme.thresholds[t])-1) #e.g., t2 between categories 2 and 3; will be skipped for binary trait
                lowerboundry_threshold  = maximum(ySparse[:,t][category_obs[t] .== (i-1)]) #lower bound for current threshold
                upperboundry_threshold  = minimum(ySparse[:,t][category_obs[t] .== i])     #upper bound for current threshold
                mme.thresholds[t][i]    = rand(Uniform(lowerboundry_threshold,upperboundry_threshold))
            end
            ####################################################################
            # step2. update lower_bound, upper_bound from thresholds
            ####################################################################
            lower_bound[t],upper_bound[t]=update_lower_upper_bound_with_threshold(mme.thresholds[t],category_obs[t])
        end
    end
end


function sample_liabilities!(mme,ycorr,lower_bound,upper_bound)
    ############################################################################
    # Goal:
    #   update liability (i.e., mme.ySparse) for censored and categorical traits.
    #   Then ycorr should be updated due to the update of mme.ySparse
    # Arguments:
    #  - ycorr: residual vector of length nInd-by-nTraits
    #  - lower_bound,upper_bound: dictionary, key is the index of the trait
    ############################################################################
    cmean          = mme.ySparse - ycorr #mean = liability - residual
    nInd           = length(mme.obsID)
    nTrait         = mme.nModels
    is_multi_trait = nTrait>1
    R              = mme.R.val
    cmean          = reshape(cmean,      nInd,nTrait)
    ySparse        = reshape(mme.ySparse,nInd,nTrait) #mme.ySparse will also be updated since reshape is a reference, not copy
    ycorr_ref      = reshape(ycorr,      nInd,nTrait) #ycorr will be updated since reshape is a reference, not copy

    if length(findall(t -> t ∈ ["categorical","categorical(binary)","censored"], mme.traits_type)) >1 #at least two traits with liability
        nGibbs = 5     #need Gibbs sampling because liability is not from the stationary distribution
    else  #only one trait with liability, Gibbs sampling is not required
        nGibbs = 1
    end
    for iter in 1:nGibbs
        for t in 1:nTrait #for each trait
            if mme.traits_type[t] ∈ ["categorical","categorical(binary)","censored"]
                index1 = t  # "1" denotes the trait for sampling liability, "2" denotes all other traits.
                index2 = deleteat!(collect(1:nTrait),index1)
                d      = ycorr_ref[:,index2] #current residuals for all other traits (d)
                #sample residual trait "1"
                μ_1    = is_multi_trait ? vec(R[index1,index2]'inv(R[index2,index2])*d') : zeros(nInd)
                σ2_1   = is_multi_trait ? R[index1,index1]-R[index1,index2]'inv(R[index2,index2])*R[index2,index1] : R #R=1 for single categorical trait
                ϵ1_lower_bound = lower_bound[t] - cmean[:,index1] #thresholds[t][whichcategory] - cmean[:,index1]
                ϵ1_upper_bound = upper_bound[t] - cmean[:,index1] #thresholds[t][whichcategory+1] - cmean[:,index1]
                for i in 1:nInd #for each individual
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


function print_single_categorical_censored_trait_example()
    ############################################################################
    # Goal:
    #   print an example for deprecated JWAS function runMCMC(categorical_trait,censored_trait)
    ############################################################################
    @error("The arguments 'categorical_trait' and  'censored_trait' has been moved to build_model(). Please check our latest example.")
    printstyled("1. Example to build model for single categorical trait:\n"; color=:red)
    printstyled("      model_equation  = \"y = intercept + x1 + x2 + x2*x3 + ID + dam + genotypes\"
      model = build_model(model_equation,categorical_trait=[\"y\"])\n"; color=:red)
    printstyled("2. Example to build model for single censored trait:\n"; color=:red)
    printstyled("      model_equation  = \"y = intercept + x1 + x2 + x2*x3 + ID + dam + genotypes\"
      model = build_model(model_equation,censored_trait=[\"y\"])\n"; color=:red)
end


function sample_from_conditional_inverse_Wishart(df,scale,binary_trait_index)
    ############################################################################
    # Goal:
    #   sample from conditional inverse Wishart distribution, where binary traits
    #   are independent with unit variance.
    # Arguments:
    #  - df: degree of freedom of the conditional inverse Wishart
    #  - scale: scale of the conditional inverse Wishart
    #  - binary_trait_index: index of binary trait, e.g,[1,3] means the 1st and 3rd traits are binary
    ############################################################################
    ntraits= size(scale,2)              #number of total traits, e.g., 5
    index2 = binary_trait_index         #index for binary traits, e.g., [1, 3 , 4]
    index1 = setdiff(1:ntraits,index2)  #index for non-binary traits, e.g.,[2, 5]
    n1     = length(index1)             #number of non-binary traits
    n2     = length(index2)             #number of binary traits

    V11    = scale[index1,index1]
    V12    = scale[index1,index2]
    V22    = scale[index2,index2]
    V22_1  = V22 - V12'*inv(V11)*V12    #n2 - by - n2
    X1     = rand(Wishart(df,V11))      # X1~W(V11,n) for non-binary trait, n1-by-n1
    μ      = vec(inv(V11)*V12)          # n1*n2 - by - 1
    Σ      = Symmetric(kron(V22_1,inv(X1))) # n1*n2 - by - n1*n2, the order of kron in paper is incorrect, check wikipedia for correct order
    X2     = rand(MvNormal(μ, Σ))       # n1*n2 - by - 1
    X2     = reshape(X2,n1,n2)          # n1 - by - n2
    T11    = inv(X1) + X2*X2'           # n1 - by - n1
    T12    = -X2                        # n1 - by - n2
    R      = [T11  T12
              T12' I(n2)]               #traits are ordered as [index1;index2]
    re_order = sortperm([index1;index2])
    R        = R[re_order,re_order]     #traits are ordered as [1,2,3,...]
    return R
end


function add_censored_trait_column!(mme,df)
    ############################################################################
    # Goal:
    #   add the column named "traitname" using trait's lower bound and upper bound
    # Arguments:
    #   - df: dataframe of phenotypes
    # Notes:
    #   - for missing censored trait, the lower bound=-Inf and higher bound=Inf,
    #     or both are missing.
    ############################################################################
    for t in 1:mme.nModels #for each trait
        if mme.traits_type[t] == "censored"
            trait_name       = mme.lhsVec[t]           #e.g., :y
            lower_bound_name = Symbol(trait_name,"_l") #e.g., :y_l
            upper_bound_name = Symbol(trait_name,"_u") #e.g., :y_u
            df[!,trait_name] = Array{Union{Missing, Float64}}(undef,nrow(df))  #initialize
            for i in 1:nrow(df) #for each individual
                if ismissing(df[i,lower_bound_name]) && ismissing(df[i,upper_bound_name]) #[missing,missing]
                    df[i,trait_name]       = missing #missing censored trait
                    df[i,lower_bound_name] = -Inf    #non-truncated
                    df[i,upper_bound_name] = Inf     #non-truncated
                elseif df[i,lower_bound_name]==-Inf && df[i,upper_bound_name]==Inf #[-Inf,Inf]
                    df[i,trait_name] = missing #missing censored trait
                elseif df[i,lower_bound_name]==-Inf #e.g.,[-Inf,3]
                    df[i,trait_name] = df[i,upper_bound_name]
                elseif df[i,upper_bound_name]==Inf #e.g.,[1,Inf]
                    df[i,trait_name] = df[i,lower_bound_name]
                else #e.g., [1,2]
                    df[i,trait_name] = rand(df[i,lower_bound_name]:df[i,upper_bound_name])
                end
            end
        end
    end
end
