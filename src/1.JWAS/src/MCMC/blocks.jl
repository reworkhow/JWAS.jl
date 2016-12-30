
#######################################################
# Pre-Check
#######################################################
function pre_check(mme,π,sol)
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
function output_MCMC_samples_setup(mme,nIter,output_samples_frequency)
  #initialize arrays to save MCMC samples
  num_samples = Int(floor(nIter/output_samples_frequency))
  init_sample_arrays(mme,num_samples)

  if mme.M != 0 #write samples for marker effects to a txt file
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
  end
  out_i = 1
  return out_i,outfile,pi
end

function output_MCMC_samples(mme,out_i,pi,sol,α,outfile)
  outputSamples(mme,sol,out_i)
  mme.samples4R[:,out_i]=vRes
  if mme.ped != 0
    mme.samples4G[:,out_i]=vec(G0)
  end
  writedlm(outfile,α')
  pi[out_i] = π
  out_i +=1
  return out_i
end
