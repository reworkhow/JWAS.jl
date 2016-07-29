function get_BV_samples(M::Array{Float64,2},marker_file;header=true)
    if header==true
        alpha=readdlm(marker_file,header=true)[1]
    else
        alpha=readdlm(marker_file)
    end
    M*alpha'
end

function get_breeding_values(M::Array{Float64,2},files...)
    nTraits = length(files)
    BV  = Array(Array{Float64,2},nTraits)
    traiti=1

    EBV = Array(Any,nTraits)
    for i in files
        BV[traiti] = get_BV_samples(M,i)
        EBV[traiti]= DataFrame(EBV=vec(mean(BV[traiti],2)),PEV=vec(var(BV[traiti],2)))
        traiti +=1
    end
    return EBV
end

function get_additive_genetic_variances(M::Array{Float64,2},files...)
    nTraits = length(files)
    BV  = Array(Array{Float64,2},nTraits)
    traiti=1
    for i in files
        BV[traiti] = get_BV_samples(M,i)
        traiti +=1
    end

    num_samples = size(BV[1],2)
    G = Array(Array{Float64,2},num_samples)
    for i = 1:num_samples
        BVi = vec(BV[1][:,i])
        for j = 2:nTraits
            BVi = [BVi BV[j][:,i]]
        end
        G[i] = cov(BVi)
    end
    Gmean = mean(G)
    Gvar  = mean(G.^2)-Gmean.^2
    df_G  = Dict("mean_of_additive_genetic_covariance_matrix"=>Gmean,
            "variance_of_additive_genetic_covariance_matrix"=>Gvar)

    return G,df_G
end


function get_BV_samples(M::Array{Float64,2},marker_file;header=true)
    if header==true
        alpha=readdlm(marker_file,header=true)[1]
    else
        alpha=readdlm(marker_file)
    end
    M*alpha'
end

function get_breeding_values(M::Array{Float64,2},files...)
    nTraits = length(files)
    BV  = Array(Array{Float64,2},nTraits)
    traiti=1

    EBV = Array(Any,nTraits)
    for i in files
        BV[traiti] = get_BV_samples(M,i)
        EBV[traiti]= DataFrame(EBV=vec(mean(BV[traiti],2)),PEV=vec(var(BV[traiti],2)))
        traiti +=1
    end
    return EBV
end

"""
    get_additive_genetic_variances(model::MME,files...)

* calculate distributions of additive genetic variances with output_marker_effects **files**
* return samples, mean, variance for additive genetic variances
"""
function get_additive_genetic_variances(model,files...)
    return get_additive_genetic_variances(model.M.genotypes,files...)
end


"""
    get_breeding_values(model::MME,files...)

* calculate distributions of breeding values with output_marker_effects **files**
* return esitimated breeding values and prediction error variances
"""
function get_breeding_values(model,files...)
    EBV=get_breeding_values(model.M.genotypes,files...)
    if model.M.obsID[1] != "NA"
        EBVwithID = Any(EBV)
        traiti=1
        for i in EBV
            ID        = DataFrame(ID=model.M.obsID)
            EBVwithID[traiti] = [ID i]
            traiti +=1
        end
        return EBVwithID
    end
    return EBV
end

function get_heritability(G::Array{Array{Float64,2},1},R::Array{Array{Float64,2},1})
  if size(G,1)!=size(R,1)
    error("Number of MCMC samples for genetic variances and residual vairances are not equal!")
  end
  h2 = Array(Array{Float64,2},size(G,1))
  for i in 1:size(G,1)
    h2[i]=G[i] ./ (G[i]+R[i])
  end

  h2mean = mean(h2)
  h2var  = mean(h2.^2)-h2mean.^2
  df_h2  = Dict("mean_of_heritability"=>h2mean,
          "variance_of_heritability"=>h2var)
  h2, df_h2
end

export get_additive_genetic_variances
export get_breeding_values
export get_heritability
