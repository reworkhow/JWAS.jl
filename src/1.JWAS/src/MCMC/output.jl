
#######################################################
# Pre-Check
#######################################################
function pre_check(mme,sol)
  if size(mme.mmeRhs)==()
      getMME(mme,df)
  end
  #starting value for sol can be provided
  if sol == false #no starting values
      sol = zeros(size(mme.mmeLhs,1))
  else            #besure type is Float64
      sol = map(Float64,sol)
  end
  solMean     = zeros(sol)
  return sol,solMean
end


#######################################################
#MCMC Samples OUTPUT
#######################################################
function output_MCMC_samples_setup(mme,nIter,output_samples_frequency,ismarker=true)
  #initialize arrays to save MCMC samples
  num_samples = Int(floor(nIter/output_samples_frequency))
  init_sample_arrays(mme,num_samples)
  out_i = 1

  if ismarker=true #write samples for marker effects to a txt file
    file_count = 1
    file_name="MCMC_samples_for_marker_effects.txt"
    while isfile(file_name)
      file_name="MCMC_samples_for_marker_effects"*"_$(file_count)"*".txt"
      file_count += 1
    end
    outfile=open(file_name,"w")

    if mme.M.markerID[1]!="NA"
        writedlm(outfile,mme.M.markerID')
    end
    pi = zeros(num_samples)#vector to save π (for BayesC)
    return out_i,outfile,pi
  else
    return out_i
  end
end

function output_MCMC_samples(mme,out_i,pi,sol,α,vRes,G0,outfile,estimatePi)
  outputSamples(mme,sol,out_i)
  mme.samples4R[:,out_i]=vRes
  if mme.ped != 0
    mme.samples4G[:,out_i]=vec(G0)
  end
  writedlm(outfile,α')
  if estimatePi==true
    pi[out_i] = π
  end
  out_i +=1
  return out_i
end

function output_MCMC_samples(mme,out_i,sol,vRes,G0)
  outputSamples(mme,sol,out_i)
  mme.samples4R[:,out_i]=vRes
  if mme.ped != 0
    mme.samples4G[:,out_i]=vec(G0)
  end
  out_i +=1
  return out_i
end


#######################################################
# Return Output Results (Dictionary)
#######################################################
function output_result(mme,solMean,output_samples_frequency,estimatePi)
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
          trmi   = i.term
          trmStr = trmi.trmStr
          output["MCMC samples for: variance of "*trmStr] = i.sampleArray
      end
  end

  if mme.M != 0
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
