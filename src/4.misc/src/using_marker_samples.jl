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
    BV  = Array{Array{Float64,2}}(nTraits)
    traiti=1

    EBV = Array{Any}(nTraits)
    for i in files
        BV[traiti] = get_BV_samples(M,i,header=header)
        EBV[traiti]= DataFrame(EBV=vec(mean(BV[traiti],2)),PEV=vec(var(BV[traiti],2)))
        traiti +=1
    end
    return EBV
end

function get_additive_genetic_variances(M::Array{Float64,2},files...;header=true)
    nTraits = length(files)
    BV  = Array{Array{Float64,2}}(nTraits)
    traiti=1
    for i in files
        BV[traiti] = get_BV_samples(M,i,header=header)
        traiti +=1
    end

    num_samples = size(BV[1],2)
    if nTraits == 1
        G=Array{Float64,1}()
    else
        G=Array{Array{Float64,2},1}()
    end
#    G = zeros(nTraits^2,num_samples)
    for i = 1:num_samples
        BVi = vec(BV[1][:,i])
        for j = 2:nTraits
            BVi = [BVi BV[j][:,i]]
        end
        push!(G,cov(BVi))
        #G[:,i] = vec(cov(BVi))
    end
    return G
end

"""
    getEBV(model,files...;header=true)

* Get esitimated breeding values and prediction error variances using samples of marker effects stored in **files**
    for individuals defined by `outputEBV(model,IDs::Array{String,1})`, defaulting to all phenotyped individuals.
"""
function getEBV(model,files...;header=true)
    EBV = get_breeding_values(model.output_genotypes,files...,header=header)
    EBVwithID = Any(EBV)
    traiti=1
    for i in EBV
        ID = DataFrame(ID=model.output_ID)
        EBVwithID[traiti] = [ID i]
        traiti +=1
    end
    return EBVwithID
end

"""
    get_breeding_values(model,filename;header=true)

* Get esitimated breeding values and prediction error variances using samples of marker effects stored in **files**
    for individuals defined by `outputEBV(model,IDs::Array{String,1})`, defaulting to all phenotyped individuals.
"""
function get_breeding_values(model,files...;header=true)
    for traiti in 1:model.nModels
        file_marker= output_file*"_"*"marker_effects_"*string(mme.lhsVec[traiti])*".txt"
        file_J     = output_file*"_"*"J"*string(mme.lhsVec[traiti])*".txt"
        BVsamples_markers = get_BV_samples(model.output_genotypes,file,header=header)



        res["EBV"*"_"*string(mme.lhsVec[traiti])]=misc.get_breeding_values(mme,file,header=(mme.M.markerID[1]!="NA")?true:false)[1]

        EBV = get_breeding_values(model.output_genotypes,files...,header=header)
        EBVwithID = Any(EBV)
        traiti=1
        for i in EBV
            ID = DataFrame(ID=model.output_ID)
            EBVwithID[traiti] = [ID i]
            traiti +=1
        end
        return EBVwithID
    end
end

# function getEBV(model)
#     for traiti in 1:mme.nModels
#
#         res["Posterior mean of marker effects"]
#
#        file= output_file*"_"*"marker_effects_"*string(mme.lhsVec[traiti])*".txt"
#        res["EBV"*"_"*string(mme.lhsVec[traiti])]=
#           misc.get_breeding_values(mme,file,header=(mme.M.markerID[1]!="NA")?true:false)[1]
#     end
#     output_others(mme)
# end




"""
    get_additive_genetic_variances(model::MME,files...;header=true)

* Get MCMC samples for additive genetic variances using samples for marker effects stored in **files**.
* Return a vector for single-trait analysis and an array of matrices for multi-trait analysis
"""
function get_additive_genetic_variances(model,files...;header=true)
    return get_additive_genetic_variances(model.M.genotypes,files...,header=header)
end

"""
    get_heritability(samples_for_genetic_variances::Array{Array{Float64,2},1},samples_for_residual_variances::Array{Array{Float64,2},1}))
* Get MCMC samples for heritabilities using MCMC samples for genetic variances and residual variances for multi-trait analysis
"""
function get_heritability(G::Array{Array{Float64,2},1},R::Array{Array{Float64,2},1})
  if size(G,1)!=size(R,1)
    error("Number of MCMC samples for genetic variances and residual variances are not equal!")
  end
  h2 = Array{Array{Float64,1}}(size(G,1))
  for i in 1:size(G,1)
    h2[i]=diag(G[i] ./ (G[i]+R[i]))
  end
  h2
end

"""
    get_heritability(samples_for_genetic_variances::Array{Float64,1},samples_for_residual_variances::Array{Float64,1}))
* Get MCMC samples for heritabilities using MCMC samples for genetic variances and residual variances for single-trait analysis
"""
function get_heritability(G::Array{Float64,1},R::Array{Float64,1})
  if length(G)!=length(R)
    error("Number of MCMC samples for genetic variances and residual variances are not equal!")
  end
  h2=G./(G+R)
  return h2
end

"""
    get_correlations(samples_for_genetic_variances::Array{Array{Float64,2},1})

* Get MCMC samples for correlations using MCMC samples for covariance matrces
"""
function get_correlations(G::Array{Array{Float64,2},1})
  gentic_correlation = Array{Array{Float64,2}}(size(G,1))
  for i in 1:size(G,1)
    varg=diag(G[i])
    gentic_correlation[i]=G[i]./sqrt.(varg*varg')
  end
  gentic_correlation
end

"""
    reformat(G::Array{Float64,2},ntraits=false)

* convert the format from a `Matrix` to a `Array of Matrices`, when the matrix is
read from text file storing samples for covariance matrices or heritabilities
* the size of the `Matrix` is
    * ntraits x nsamples for MCMC samples for heritbilities
    * (ntraits^2) x nsamples for MCMC samples for covariance matrices
* the `Array of Matrices` is an array (length=nsamples) of
    *  matrices of size ntraits X ntraits for covariance matrices
    *  vectors of length ntraits for heritabilities
"""
function reformat(G::Array{Float64,2},ntraits=false)
    if ntraits==size(G,2) #for heritability size(G) = ntraits x nsamples
        Gnew = Array{Array{Float64,1}}(size(G,1))
        for i in 1:size(G,1)
            Gnew[i]=vec(G[i,:])
        end
        return Gnew
    end

    if ntraits==false || ntraits==Int(sqrt(size(G,2)))
    #for covariance matrix, size(G) = (ntrait^2) x nsamples
        Gnew = Array{Array{Float64,2}}(size(G,1))
        ntraits=Int(sqrt(size(G,2)))
        for i in 1:size(G,1)
            Gnew[i]=reshape(G[i,:],ntraits,ntraits)
        end
        return Gnew
    end
end

"""
    reformat(G::Array{Array{Float64,2},1},1})

* convert the format from a `Array of Matrices` to a `Matrix`, then the `Matrix` chance
be write to text files to save samples for covariance matrices or heritabilities
* the `Array of Matrices` is an array (length=nsamples) of
    *  matrices of size ntraits X ntraits for covariance matrices (`Array{Array{Float64,2},1}`)
    *  vectors of length ntraits for heritabilities (`Array{Array{Float64,1},1}`)
* the size of the `Matrix` is
    * ntraits x nsamples for MCMC samples for heritbilities
    * (ntraits^2) x nsamples for MCMC samples for covariance matrices
"""
function reformat(G::Array{Array{Float64,2},1})
    Gnew = zeros(length(G),length(G[1]))
    for i in 1:length(G)
        Gnew[i,:]=vec(G[i])
    end
    Gnew
end

function reformat(G::Array{Array{Float64,1},1})
    Gnew = zeros(length(G),length(G[1]))
    for i in 1:length(G)
        Gnew[i,:]=vec(G[i])
    end
    Gnew
end

export get_additive_genetic_variances
export get_breeding_values
export get_heritability
export get_correlations
export reformat
