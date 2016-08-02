function get_BV_samples(M::Array{Float64,2},marker_file;header=true)
    if header==true
        alpha=readdlm(marker_file,header=true)[1]
    else
        alpha=readdlm(marker_file)
    end
    M*alpha'
end

function get_breeding_values(M::Array{Float64,2},files...;header=true)
    nTraits = length(files)
    BV  = Array(Array{Float64,2},nTraits)
    traiti=1

    EBV = Array(Any,nTraits)
    for i in files
        BV[traiti] = get_BV_samples(M,i,header=header)
        EBV[traiti]= DataFrame(EBV=vec(mean(BV[traiti],2)),PEV=vec(var(BV[traiti],2)))
        traiti +=1
    end
    return EBV
end

function get_additive_genetic_variances(M::Array{Float64,2},files...;header=true)
    nTraits = length(files)
    BV  = Array(Array{Float64,2},nTraits)
    traiti=1
    for i in files
        BV[traiti] = get_BV_samples(M,i,header=header)
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
    return G
end

"""
    get_breeding_values(model::MME,files...;header=true)

* Get esitimated breeding values and prediction error variances using samples of marker effects stored in **files**.
"""
function get_breeding_values(model,files...;header=true)
    EBV=get_breeding_values(model.M.genotypes,files...,header=header)
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

"""
    get_additive_genetic_variances(model::MME,files...;header=true)

* Get MCMC samples for additive genetic variances using samples of marker effects stored in **files**.
"""
function get_additive_genetic_variances(model,files...;header=true)
    return get_additive_genetic_variances(model.M.genotypes,files...,header=header)
end

"""
    get_heritability(samples_for_genetic_variances::Array{Array{Float64,2},1},samples_for_residual_vairances::Array{Array{Float64,2},1}))

* Get MCMC samples for heritabilities.
"""
function get_heritability(G::Array{Array{Float64,2},1},R::Array{Array{Float64,2},1})
  if size(G,1)!=size(R,1)
    error("Number of MCMC samples for genetic variances and residual vairances are not equal!")
  end
  h2 = Array(Array{Float64,1},size(G,1))
  for i in 1:size(G,1)
    h2[i]=diag(G[i] ./ (G[i]+R[i]))
  end
  h2
end

"""
    get_genetic_correlations(samples_for_genetic_variances::Array{Array{Float64,2},1})

* Get MCMC samples for genetic_correlations
"""
function get_genetic_correlations(G::Array{Array{Float64,2},1})
  gentic_correlation = Array(Array{Float64,2},size(G,1))
  for i in 1:size(G,1)
    varg=diag(G[i])
    gentic_correlation[i]=G[i]./sqrt(varg*varg')
  end
  gentic_correlation
end


export get_additive_genetic_variances
export get_breeding_values
export get_heritability
export get_genetic_correlations
