################################################################################
#User-interface to output MCMC Samples for specific variables
################################################################################
"""
    outputMCMCsamples(mme::MME,trmStr::AbstractString...)

Get samples for specific variables.
"""
function outputMCMCsamples(mme::MME,trmStr::AbstractString...)
    for i in trmStr
      outputSamplesFor(mme,i)
    end
end

################################################################################
#Wraps for Output MCMC Samples
################################################################################
function output_MCMC_samples_setup(mme,nIter,output_samples_frequency,ismarker=true)
  #initialize arrays to save MCMC samples
  num_samples = Int(floor(nIter/output_samples_frequency))
  init_sample_arrays(mme,num_samples)
  out_i = 1

  if ismarker==true #write samples for marker effects to a txt file
    file_count = 1
    file_name="MCMC_samples_for_marker_effects.txt"
    while isfile(file_name)
      file_name="MCMC_samples_for_marker_effects"*"_$(file_count)"*".txt"
      file_count += 1
    end
    outfile=open(file_name,"w")

    if mme.M.markerID[1]!="NA"
        writedlm(outfile,transpose(mme.M.markerID))
    end
    pi = zeros(num_samples)#vector to save π (for BayesC)
    return out_i,outfile,pi
  else
    return out_i #for conventional analyses (no markers)
  end
end

function output_MCMC_samples(mme,out_i,sol,vRes,G0,
                             α=false,pi=false,outfile=false,estimatePi=false)
  outputSamples(mme,sol,out_i)
  mme.samples4R[:,out_i]=vRes
  if mme.ped != 0
    mme.samples4G[:,out_i]=vec(G0)
  end
  if α != false
    writedlm(outfile,α')
    if estimatePi==true
      pi[out_i] = π
    end
  end
  out_i +=1
  return out_i
end
################################################################################
#function blocks
################################################################################

#define samples for WHICH location parameters to output
function outputSamplesFor(mme::MME,trmStr::AbstractString)
    #add model number => "1:age"
    res = []
    for (m,model) = enumerate(mme.modelVec)
        strVec  = split(model,['=','+'])
        strpVec = [strip(i) for i in strVec]
        if trmStr in strpVec
            res = [res;string(m)*":"*trmStr]
        end
    end #"age"->"1:age","2:age"

    for trmStr in res
        trm     = mme.modelTermDict[trmStr]
        samples = MCMCSamples(trm,Array(Float64,1,1))
        push!(mme.outputSamplesVec,samples)
    end
end

#init
function init_sample_arrays(mme::MME,niter)
    #varaince components for residual
    mme.samples4R = zeros(mme.nModels^2,niter)

    #variance components for random polygenic effects
    if mme.ped != 0
        mme.samples4G = zeros(length(mme.pedTrmVec)^2,niter)
    end

    #location parameters for fixed and random effects except markers
    for i in  mme.outputSamplesVec #resize
        trmi = i.term
        i.sampleArray = zeros(trmi.nLevels,niter)
    end

    #variance components for iid random effects
    for i in  mme.rndTrmVec #resize to be size of nTraits
        i.sampleArray = zeros(mme.nModels^2,niter)#Bug maybe many diff
    end
end

#output samples for location parameers
function outputSamples(mme::MME,sol,iter::Int64)
    for i in  mme.outputSamplesVec
        trmi = i.term
        startPosi  = trmi.startPos
        endPosi    = startPosi + trmi.nLevels - 1
        i.sampleArray[:,iter] = sol[startPosi:endPosi]
    end
    for effect in  mme.rndTrmVec
        effect.sampleArray[iter] = effect.vcNew
    end
end

#function to replace Array{SubString{String}' for issue 8
function transpose(vec::Array{SubString{String},1})
    lvec=length(vec)
    res =Array(String,1,lvec)
    for i in 1:lvec
        res[1,i]=vec[i]
    end
    return res
end
