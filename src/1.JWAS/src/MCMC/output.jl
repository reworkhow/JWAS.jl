################################################################################
#User-interface to set training and testing individuals
################################################################################
#function set_training(model,IDs)
#    IDs = map(String,vec(IDs)) #Array{String,1}
#    model.training_ID=IDs
#end


"""
    outputEBV(model,IDs::Array;PEV=false)

Output estimated breeding values and prediction error variances (defaulting to false) for IDs.
"""
function outputEBV(model,IDs;PEV=false)
    #print("Estimated breeding values and prediction error variances will be included in the output.")
    IDs = map(string,vec(IDs)) #Array{String,1}
    model.output_ID=IDs
end

function set_testing(model,IDs;PEV=false)
    outputEBV(model,IDs,PEV=PEV)
end

"""
    get_outputX_others(model)

Make incidence matrices for effects involve in EBV inclung J, ϵ, pedTrmVec except marker covariates
"""
function get_outputX_others(model,single_step_analysis)
    #trick to avoid errors (PedModule.getIDs(ped) [nongeno ID;geno ID])
    if single_step_analysis == true
        model.output_X["ϵ"]=mkmat_incidence_factor(model.output_ID,
                            PedModule.getIDs(model.ped))[:,1:length(model.ped.setNG)]
        #model.output_X["J"] is in SSBRrun
    end
    if model.pedTrmVec != 0
        for i in model.pedTrmVec
            model.output_X[i]=mkmat_incidence_factor(model.output_ID,
                                                     model.modelTermDict[i].names)
        end
    end
end

function getEBV(model::MME)
    if model.output_ID != 0 &&  (model.pedTrmVec != 0 || model.M != 0)
        if model.nModels == 1
            EBV = model.output["EBV"*"_"*string(model.lhsVec[1])]
        else
            EBV =[] #Array{Any,1}(undef,0)
            for traiti in 1:model.nModels
                push!(EBV,model.output["EBV"*"_"*string(model.lhsVec[traiti])])
            end
        end
    end
    return EBV
end

#only for genotyped individuals
#same order for markers
#first column is individual ID
function getEBV(model::MME,genotypefile::AbstractString;header=true,separator=',')
    myfile = open(genotypefile)
    #get number of columns
    row1   = split(readline(myfile),[separator,'\n'],keepempty=false)
    #set types for each column and get number of markers
    ncol= length(row1)
    etv = Array{DataType}(undef,ncol)
    fill!(etv,Float64)
    etv[1]=String
    close(myfile)
    df = readtable(genotypefile, eltypes=etv, separator = separator, header=header)
    obsID     = map(String,df[1]) #convert from Array{Union{String, Missings.Missing},1} to String #redundant actually
    genotypes = map(Float64,convert(Array,df[2:end]))
    genotypes = genotypes .- model.M.alleleFreq

    if model.nModels == 1
        marker_effects=map(Float64,
                        model.output["Posterior mean of marker effects"][:,end])
        EBV = [obsID genotypes*marker_effects]
    else
        EBV =[] #Array{Any,1}(undef,0)
        for traiti in 1:model.nModels
            marker_effects=map(Float64,
                model.output["Posterior mean of marker effects"][traiti][:,end])

            push!(EBV,[obsID genotypes*marker_effects])
        end
    end
    return EBV
end

function getEBV(model::MME,genotypes::Array{Float64,2})
    genotypes = genotypes .- model.M.alleleFreq

    if model.nModels == 1
        marker_effects=map(Float64,
                        model.output["Posterior mean of marker effects"][:,end])
        EBV = genotypes*marker_effects
    else
        EBV =[] #Array{Any,1}(undef,0)
        for traiti in 1:model.nModels
            marker_effects=map(Float64,
                model.output["Posterior mean of marker effects"][traiti][:,end])

            push!(EBV,genotypes*marker_effects)
        end
    end
    return EBV
end



export outputEBV
