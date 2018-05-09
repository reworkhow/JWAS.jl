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
function output_result(mme,solMean,output_samples_frequency,
                       meanAlpha=false,estimatePi=false,pi=false,mean_pi=false)
  output = Dict()
  output["Posterior mean of location parameters"] = [getNames(mme) solMean]
  if output_samples_frequency != 0
      output["MCMC samples for residual variance"]    = mme.samples4R
      if mme.ped != 0
          output["MCMC samples for polygenic effects var-cov parameters"] = mme.samples4G
      end
      for i in  mme.outputSamplesVec
          trmi   = i.term
          trmStr = trmi.trmStr
          output["MCMC samples for: "*trmStr] = [getNames(trmi) i.sampleArray]
      end
      for i in  mme.rndTrmVec
          trmi   = i.term_array[1]
          trmStr = trmi.trmStr
          output["MCMC samples for: variance of "*trmStr] = i.sampleArray
      end
  end

  if mme.M != 0 && meanAlpha != false
    if mme.M.markerID[1]!="NA"
        markerout=[mme.M.markerID meanAlpha]
    else
        markerout= meanAlpha
    end

    output["Posterior mean of marker effects"] = markerout
    if estimatePi == true
        output["Posterior mean of Pi"] = mean_pi
        if  output_samples_frequency != 0
            output["MCMC samples for: π"] = pi
        end
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
