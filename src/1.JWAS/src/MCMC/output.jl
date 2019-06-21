################################################################################
# Return Output Results (Dictionary)
################################################################################
function output_result(mme,solMean,meanVare,G0Mean,output_samples_frequency,
                       meanAlpha,meanVara,estimatePi,mean_pi,estimateScale,meanScaleVar,output_file)
  output = Dict()
  if mme.output_ID != 0 &&  (mme.pedTrmVec != 0 || mme.M != 0)
      for traiti in 1:mme.nModels
          output["EBV"*"_"*string(mme.lhsVec[traiti])]=zeros(length(mme.output_ID))
      end
  end

  location_parameters = reformat2DataFrame([getNames(mme) solMean])
  output["Posterior mean of location parameters"] = location_parameters
  output["Posterior mean of residual variance"]   = meanVare
  if mme.pedTrmVec != 0
    output["Posterior mean of polygenic effects covariance matrix"]=G0Mean
    if mme.output_ID != 0
        for pedtrm in mme.pedTrmVec
            traiti, effect = split(pedtrm,':')
            sol_pedtrm     = map(Float64,location_parameters[(location_parameters[:Effect].==effect).&(location_parameters[:Trait].==traiti),:Estimate])
            EBV_pedtrm     = mme.output_X[pedtrm]*sol_pedtrm
            output["EBV"*"_"*string(mme.lhsVec[parse(Int64,traiti)])] += EBV_pedtrm
        end
    end
  end

  if mme.M != 0
    if mme.nModels == 1
        #make single-trait estimated marker effects as type
        #*array of array* (same to multi-trait) for coding simplicity
        meanAlpha=[meanAlpha]
    end
    markerout        = []
    if mme.M.markerID[1]!="NA"
        for markerArray in meanAlpha
          push!(markerout,[mme.M.markerID markerArray])
        end
    else
        for markerArray in meanAlpha
          push!(markerout,markerArray)
        end
    end

    output["Posterior mean of marker effects"] = (mme.nModels==1) ? markerout[1] : markerout
    output["Posterior mean of marker effects variance"] = meanVara
    if estimatePi == true
        output["Posterior mean of Pi"] = mean_pi
    end
    if estimateScale == true
        output["Posterior mean of ScaleEffectVar"] = meanScaleVar
    end

    if mme.output_ID != 0
        for traiti in 1:mme.nModels
            EBV_markers  = mme.output_genotypes*meanAlpha[traiti] #fixed for mt
            output["EBV"*"_"*string(mme.lhsVec[traiti])] += EBV_markers
        end
    end
  end

  if mme.output_ID != 0 && haskey(mme.output_X,"J") #single-step analyis
      for traiti in 1:mme.nModels
          sol_J        = map(Float64,location_parameters[(location_parameters[:Effect].=="J").&(location_parameters[:Trait].==string(traiti)),:Estimate])[1]
          sol_ϵ        = map(Float64,location_parameters[(location_parameters[:Effect].=="ϵ").&(location_parameters[:Trait].==string(traiti)),:Estimate])
          EBV_J        = mme.output_X["J"]*sol_J
          EBV_ϵ        = mme.output_X["ϵ"]*sol_ϵ
          output["EBV"*"_"*string(mme.lhsVec[traiti])] += (EBV_J+EBV_ϵ)
      end
  end

  if mme.output_ID != 0 &&  (mme.pedTrmVec != 0 || mme.M != 0)
      for traiti in 1:mme.nModels
          EBV = output["EBV"*"_"*string(mme.lhsVec[traiti])]
          if EBV != zeros(length(mme.output_ID))
              myEBV = "EBV"*"_"*string(mme.lhsVec[traiti])
              output[myEBV]= DataFrame([mme.output_ID EBV],[:ID,:Estimate])
          end
          if mme.MCMCinfo.output_PEV == true
              EBVsamplesfile = output_file*"_"*myEBV*".txt"
              EBVsamples,IDs = readdlm(EBVsamplesfile,',',header=true)
              if vec(IDs) == mme.output_ID
                  output[myEBV][:PEV]= vec(var(EBVsamples,dims=1))
              else
                  error("The EBV file is wrong.")
              end
          end
      end
  end

  return output
end

################################################################################
# Reformat Output Array to DataFrame
################################################################################
function reformat2DataFrame(res::Array)
    ##SLOW, 130s
    #out_names=[strip(i) for i in split(res[1,1],':',keepempty=false)]
    #for rowi in 2:size(res,1)
    #    out_names=[out_names [strip(i) for i in split(res[rowi,1],':',keepempty=false)]]#hcat two vectors
    #end
    out_names = Array{String}(undef,size(res,1),3)
    for rowi in 1:size(res,1)
        out_names[rowi,:]=[strip(i) for i in split(res[rowi,1],':',keepempty=false)]
    end

    if size(out_names,2)==1 #convert vector to matrix
        out_names = reshape(out_names,length(out_names),1)
    end
    #out_names=permutedims(out_names,[2,1]) #rotate
    out_values=map(Float64,res[:,2])
    out=[out_names out_values]
    out = DataFrame(out, [:Trait, :Effect, :Level, :Estimate])
    return out
end

"""
    outputEBV(model,IDs::Array;PEV=false)

Output estimated breeding values and prediction error variances (defaulting to false) for IDs.
"""
function outputEBV(model,IDs;PEV=false)
    #print("Estimated breeding values and prediction error variances will be included in the output.")
    IDs = map(string,vec(IDs)) #Array{String,1}
    model.output_ID=IDs
end

"""
    get_outputX_others(model)

Make incidence matrices for effects involved in calculation of EBV including J, ϵ, pedTrmVec except marker covariates.
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

"""
    getEBV(model::MME)

Get estimated breeding values for individuals defined by `outputEBV()` (defaulting to all genotyped individuals).
"""
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

"""
    getEBV(model::MME,genotypefile::AbstractString;header=true,separator=',')

Get estimated breeding values for individuals with genotypes defined by `genotypefile`.
The genotype file format is same to that used in `add_genotypes()`. The same order for
markers (columns) in  `genotypefile` as in `add_genotypes()` is assumed for simplicity
now. The file format is as follows,

```
Animal,marker1,marker2,marker3,marker4,marker5
A1,1,0,1,1,1
C3,1,2,0,1,0
B2,0,0,2,1,1
```
"""
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
    obsID     = map(string,df[1]) #convert from Array{Union{String, Missings.Missing},1} to String #redundant actually
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

"""
    getEBV(model::MME,genotypes::Array{Float64,2})

Get estimated breeding values for individuals with genotypes. The genotypes are
provided as an Array. The same order for markers (columns)in `genotypes::Array{Float64,2}`
as in `add_genotypes()` is assumed for simplicity now.
"""
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

function getEBV(mme,sol,α,traiti)
    traiti = string(traiti)
    EBV=zeros(length(mme.output_ID))

    location_parameters = reformat2DataFrame([getNames(mme) sol])
    if mme.pedTrmVec != 0
        for pedtrm in mme.pedTrmVec
            mytrait, effect = split(pedtrm,':')
            if mytrait == traiti
                sol_pedtrm     = map(Float64,location_parameters[(location_parameters[:Effect].==effect).&(location_parameters[:Trait].==traiti),:Estimate])
                EBV_pedtrm     = mme.output_X[pedtrm]*sol_pedtrm
                EBV += EBV_pedtrm
            end
        end
    end

    if mme.M != 0
        EBV += mme.output_genotypes*α
    end

    if haskey(mme.output_X,"J") #single-step analyis
        for traiti in 1:mme.nModels
            sol_J  = map(Float64,location_parameters[(location_parameters[:Effect].=="J").&(location_parameters[:Trait].==string(traiti)),:Estimate])[1]
            sol_ϵ  = map(Float64,location_parameters[(location_parameters[:Effect].=="ϵ").&(location_parameters[:Trait].==string(traiti)),:Estimate])
            EBV_J  = mme.output_X["J"]*sol_J
            EBV_ϵ  = mme.output_X["ϵ"]*sol_ϵ
            EBV   += (EBV_J+EBV_ϵ)
        end
    end

    return EBV
end
