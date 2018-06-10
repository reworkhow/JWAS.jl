################################################################################
# Pre-Check
################################################################################
function pre_check(mme,df,sol)
#  if size(mme.mmeRhs)==()
      getMME(mme,df)
#  end
  #starting value for sol can be provided
  if sol == false #no starting values
      sol = zeros(size(mme.mmeLhs,1))
  else            #besure type is Float64
      sol = map(Float64,sol)
  end
  return sol
end

################################################################################
# Return Output Results (Dictionary)
################################################################################
function output_result(mme,solMean,meanVare,G0Mean,output_samples_frequency,
                       meanAlpha,meanVara,estimatePi,mean_pi,output_file="MCMC_samples")
  output = Dict()
  for traiti in 1:mme.nModels
      output["EBV"*"_"*string(mme.lhsVec[traiti])]=zeros(length(mme.output_ID))
  end

  location_parameters = reformat2DataFrame([getNames(mme) solMean])
  output["Posterior mean of location parameters"] = location_parameters
  output["Posterior mean of residual variance"]   = meanVare
  if mme.pedTrmVec != 0
    output["Posterior mean of polygenic effects covariance matrix"]=G0Mean

    for pedtrm in mme.pedTrmVec
        traiti, effect = split(pedtrm,':')
        sol_pedtrm     = map(Float64,location_parameters[(location_parameters[:Effect].==effect).&(location_parameters[:Trait].==traiti),:Estimate])
        EBV_pedtrm     = mme.output_X[pedtrm]*sol_pedtrm
        output["EBV"*"_"*string(mme.lhsVec[parse(Int64,traiti)])] += EBV_pedtrm
    end
  end

  #samples for non-marker effects
  if output_samples_frequency != 0
      for i in  mme.outputSamplesVec
          trmi   = i.term
          trmStr = trmi.trmStr
          writedlm(output_file*"_"*trmStr*".txt",[transubstrarr(getNames(trmi))
                                                  i.sampleArray])
      end
  end

  if mme.M != 0
    if mme.nModels == 1
        meanAlpha=[meanAlpha] # make st array of array
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

    output["Posterior mean of marker effects"] = (mme.nModels==1)?markerout[1]:markerout
    output["Posterior mean of marker effects variance"] = meanVara
    if estimatePi == true
        output["Posterior mean of Pi"] = mean_pi
    end

    for traiti in 1:mme.nModels
        EBV_markers  = mme.output_genotypes*meanAlpha[traiti] #fixed for mt
        output["EBV"*"_"*string(mme.lhsVec[traiti])] += EBV_markers
    end
  end

  if haskey(mme.output_X,"J") #single-step analyis
      for traiti in 1:mme.nModels
          sol_J        = map(Float64,location_parameters[(location_parameters[:Effect].=="J").&(location_parameters[:Trait].==string(traiti)),:Estimate])[1]
          sol_ϵ        = map(Float64,location_parameters[(location_parameters[:Effect].=="ϵ").&(location_parameters[:Trait].==string(traiti)),:Estimate])
          EBV_J        = mme.output_X["J"]*sol_J
          EBV_ϵ        = mme.output_X["ϵ"]*sol_ϵ
          output["EBV"*"_"*string(mme.lhsVec[traiti])] += (EBV_J+EBV_ϵ)
      end
  end

  for traiti in 1:mme.nModels
      EBV = output["EBV"*"_"*string(mme.lhsVec[traiti])]
      if EBV != zeros(length(mme.output_ID))
          output["EBV"*"_"*string(mme.lhsVec[traiti])]= [mme.output_ID EBV]
      end
  end
  return output
end
################################################################################
# Reformat Output Array to DataFrame
################################################################################
function reformat2DataFrame(res::Array)
    out_names=[strip(i) for i in split(res[1,1],':',keep=false)]
    for rowi in 2:size(res,1)
        out_names=[out_names [strip(i) for i in split(res[rowi,1],':',keep=false)]]
    end
    out_names=permutedims(out_names,[2,1])
    out_values=map(Float64,res[:,2])
    out=[out_names out_values]
    out = DataFrame(out, [:Trait, :Effect, :Level, :Estimate],)
    return out
end

################################################################################
#Convert Genetic variance to marker effect variance based on Pi
################################################################################
function genetic2marker(M::Genotypes,Pi::Dict)
  G=M.G #genetic variance
  nTraits = size(G,1)
  denom   = zeros(nTraits,nTraits)
  for i in 1:nTraits
    for j in i:nTraits
      pi_selected = filter((k,v)->k[i]==1.0 && k[j]==1.0,Pi)
      denom[i,j] = M.sum2pq*sum(values(pi_selected))
      denom[j,i] = denom[i,j]
    end
  end
  M.G = M.G ./ denom
  M.G_is_marker_variance = true
end

function genetic2marker(M::Genotypes,π::Float64)
    M.G=M.G/((1-π)*M.sum2pq)
end
