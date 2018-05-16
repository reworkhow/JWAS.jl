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
  solMean     = zeros(sol)
  return sol,solMean
end

################################################################################
# Return Output Results (Dictionary)
################################################################################
function output_result(mme,solMean,meanVare,G0Mean,output_samples_frequency,
                       meanAlpha,meanVara,estimatePi,mean_pi,output_file="MCMC_samples")
  output = Dict()
  output["Posterior mean of location parameters"] = [getNames(mme) solMean]
  output["Posterior mean of residual variance"]   = meanVare
  if mme.ped != 0
    output["Posterior mean of Polygenic effects covariance matrix"]=G0Mean
  end

  if output_samples_frequency != 0
      for i in  mme.outputSamplesVec
          trmi   = i.term
          trmStr = trmi.trmStr
          writedlm(output_file*"_"*trmStr*".txt",[transubstrarr(getNames(trmi))
                                                  i.sampleArray])
      end
  end

  if mme.M != 0
    if mme.M.markerID[1]!="NA"
        markerout=[mme.M.markerID meanAlpha]
    else
        markerout= meanAlpha
    end

    output["Posterior mean of marker effects"] = markerout
    output["Posterior mean of marker effects variance"] = meanVara
    if estimatePi == true
        output["Posterior mean of Pi"] = mean_pi
    end
  end

  return output
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
