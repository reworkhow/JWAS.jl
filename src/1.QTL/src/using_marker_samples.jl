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
    G = Array(Any,num_samples)
    for i = 1:num_samples
        BVi = vec(BV[1][:,i])
        for j = 2:nTraits
            BVi = [BVi BV[j][:,i]]
        end
        G[i] = cov(BVi)
    end
    Gmean = mean(G)
    Gvar  = mean(G.^2)-Gmean.^2
    return G,Gmean,Gvar
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

export get_additive_genetic_variances
export get_breeding_values
