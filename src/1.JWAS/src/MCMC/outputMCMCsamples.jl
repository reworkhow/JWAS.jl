################################################################################
#User-interface to output MCMC Samples for specific variables
################################################################################
"""
    outputMCMCsamples(mme::MME,trmStr::AbstractString...)

Get MCMC samples for specific location parameters.
"""
function outputMCMCsamples(mme::MME,trmStr::AbstractString...)
    for i in trmStr
      outputSamplesFor(mme,i)
    end
end

#define samples for WHICH location parameters to output
function outputSamplesFor(mme::MME,trmStr::AbstractString)
    #add model number => "1:age"
    #"age" may be in trait 1 but not 2
    res = []
    for (m,model) = enumerate(mme.modelVec)
        strVec  = split(model,['=','+'])
        strpVec = [strip(i) for i in strVec]
        if trmStr in strpVec || trmStr in ["J","ϵ"]
            res = [res;string(m)*":"*trmStr]
        end
    end #"age"->"1:age","2:age"

    for trmStr in res
        trm     = mme.modelTermDict[trmStr]
        samples = MCMCSamples(trm,Array{Float64}(1,1))
        push!(mme.outputSamplesVec,samples)
    end
end

################################################################################
#Set-Up files and arrays to save MCMC Samples
################################################################################
function output_MCMC_samples_setup(mme,nIter,output_samples_frequency,file_name="MCMC_samples")
  #initialize arrays to save MCMC samples
  ntraits     = size(mme.lhsVec,1)
  num_samples = Int(floor(nIter/output_samples_frequency))
  init_sample_arrays(mme,num_samples)
  out_i = 1

  outfile=Dict{String,IOStream}()

  outvar = ["residual_variance"]
  if mme.pedTrmVec != 0
      push!(outvar,"polygenic_effects_variance")
  end
  if mme.M !=0 #write samples for marker effects to a text file
      for traiti in 1:ntraits
          push!(outvar,"marker_effects_"*string(mme.lhsVec[traiti]))
      end
      push!(outvar,"marker_effects_variances")
      push!(outvar,"pi")
  end
  #non-marker random effects variances
  for i in  mme.rndTrmVec
      trmStri   = split(i.term_array[1],':')[end]
      push!(outvar,trmStri*"_variances")
  end

  for i in outvar
      file_i    = file_name*"_"*i*".txt"
      if isfile(file_i)
        warn("The file "*file_i*" already exists!!! It is overwritten by the new output.")
      else
        info("The file "*file_i*" is created to save MCMC samples for "*i*".")
      end
      outfile[i]=open(file_i,"w")
  end

  #add headers
  mytraits=map(string,mme.lhsVec)
  residual_header = repeat(mytraits,inner=length(mytraits)).*"_".*repeat(mytraits,outer=length(mytraits))
  writedlm(outfile["residual_variance"],transubstrarr(residual_header),',')

  for effect in  mme.rndTrmVec
    trmStri   = split(effect.term_array[1],':')[end]                  #x2
    thisheader= repeat(effect.term_array,inner=length(effect.term_array)).*"_".*repeat(effect.term_array,outer=length(effect.term_array))
    writedlm(outfile[trmStri*"_variances"],transubstrarr(thisheader),',') #1:x2_1:x2,1:x2_2:x2,2:x2_1:x2,2:x2_2:x2
  end

  if mme.M !=0 && mme.M.markerID[1]!="NA"
      for traiti in 1:ntraits
          writedlm(outfile["marker_effects_"*string(mme.lhsVec[traiti])],transubstrarr(mme.M.markerID),',')
      end
  end
  if mme.pedTrmVec != 0
      pedtrmvec  = mme.pedTrmVec
      thisheader = repeat(pedtrmvec,inner=length(pedtrmvec)).*"_".*repeat(pedtrmvec,outer=length(pedtrmvec))
      writedlm(outfile["polygenic_effects_variance"],transubstrarr(thisheader),',')
  end

  return out_i,outfile
end

#init sample arrays
function init_sample_arrays(mme::MME,niter)
    #location parameters for fixed and random effects except markers
    for i in  mme.outputSamplesVec #resize
        trmi = i.term
        i.sampleArray = zeros(niter,trmi.nLevels)
    end
    #WRITE to FILES
    #variance components for iid random effects
    #for i in  mme.rndTrmVec #resize to be size of nTraits
    #    i.sampleArray = zeros(niter,length(i.term_array)^2)
    #end
end

################################################################################
#Output MCMC Samples every output_samples_frequency steps
################################################################################
function output_MCMC_samples(mme,out_i,sol,vRes,G0,
                             π=false,
                             α=false,
                             locusEffectVar=false,
                             outfile=false)
  ntraits     = size(mme.lhsVec,1)
  outputSamples(mme,sol,out_i)
  #random effects variances
  for effect in  mme.rndTrmVec
    trmStri   = split(effect.term_array[1],':')[end]
    writedlm(outfile[trmStri*"_variances"],vec(inv(effect.Gi))',',')
  end

  writedlm(outfile["residual_variance"],issubtype(typeof(vRes),Number)?vRes:vec(vRes)',',')
  #mme.samples4R[out_i,:]=vRes

  if mme.pedTrmVec != 0
    #mme.samples4G[out_i,:]=vec(G0)
    writedlm(outfile["polygenic_effects_variance"],vec(G0)',',')
  end
  if α != false && outfile != false
      if ntraits == 1
          writedlm(outfile["marker_effects_"*string(mme.lhsVec[1])],α',',')
      else
          for traiti in 1:ntraits
              writedlm(outfile["marker_effects_"*string(mme.lhsVec[traiti])],α[traiti]',',')
          end
      end
      if locusEffectVar!=false
          writedlm(outfile["marker_effects_variances"],locusEffectVar',',')
      end
      writedlm(outfile["pi"],π,',')
      if !issubtype(typeof(π),Number)#add a blank line
          println(outfile["pi"])
      end
  end
  out_i +=1
  return out_i
end

#output samples for location parameers
function outputSamples(mme::MME,sol,iter::Int64)
    for i in  mme.outputSamplesVec
        trmi = i.term
        startPosi  = trmi.startPos
        endPosi    = startPosi + trmi.nLevels - 1
        i.sampleArray[iter,:] = sol[startPosi:endPosi]
    end
end


# for vec::Array{AbstractString,1} and vec::Array{String,1}
# and Array{SubString{String} for issue 8
function transubstrarr(vec)
    lvec=length(vec)
    res =Array{String}(1,lvec)
    for i in 1:lvec
        res[1,i]=vec[i]
    end
    return res
end
