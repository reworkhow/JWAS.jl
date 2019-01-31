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
        push!(mme.outputSamplesVec,trm)
    end
end

################################################################################
#Set-Up files to save MCMC Samples
################################################################################
function output_MCMC_samples_setup(mme,nIter,output_samples_frequency,file_name="MCMC_samples")
  ntraits     = size(mme.lhsVec,1)

  outfile = Dict{String,IOStream}()

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

  #non-marker location parameters
  for trmi in  mme.outputSamplesVec
      trmStr = trmi.trmStr
      push!(outvar,trmStr)
  end

  #non-marker random effects variances
  for i in  mme.rndTrmVec
      trmStri   = split(i.term_array[1],':')[end]
      push!(outvar,trmStri*"_variances")
  end

  #EBV
  if mme.MCMCinfo.outputEBV == true
      for traiti in 1:ntraits
          push!(outvar,"EBV_"*string(mme.lhsVec[traiti]))
      end
      if mme.MCMCinfo.output_heritability == true && mme.MCMCinfo.single_step_analysis == false
          push!(outvar,"genetic_variance")
          push!(outvar,"heritability")
      end
  end

  for i in outvar
      file_i    = file_name*"_"*i*".txt"
      if isfile(file_i)
        printstyled("The file "*file_i*" already exists!!! It is overwritten by the new output.\n",bold=false,color=:red)
      else
        printstyled("The file "*file_i*" is created to save MCMC samples for "*i*".\n",bold=false,color=:green)
      end
      outfile[i]=open(file_i,"w")
  end

  #add headers
  mytraits=map(string,mme.lhsVec)
  varheader = repeat(mytraits,inner=length(mytraits)).*"_".*repeat(mytraits,outer=length(mytraits))
  writedlm(outfile["residual_variance"],transubstrarr(varheader),',')

  for trmi in  mme.outputSamplesVec
      trmStr = trmi.trmStr
      writedlm(outfile[trmStr],transubstrarr(getNames(trmi)),',')
  end

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

  if mme.MCMCinfo.outputEBV == true
      for traiti in 1:ntraits
          writedlm(outfile["EBV_"*string(mme.lhsVec[traiti])],transubstrarr(mme.output_ID),',')
      end
      if mme.MCMCinfo.output_heritability == true && mme.MCMCinfo.single_step_analysis == false
          writedlm(outfile["genetic_variance"],transubstrarr(varheader),',')
          writedlm(outfile["heritability"],transubstrarr(map(string,mme.lhsVec)),',')
      end
  end

  return outfile
end

#output samples for location parameers
function output_location_parameters_samples(mme::MME,sol,outfile)
    for trmi in  mme.outputSamplesVec
        trmStr = trmi.trmStr
        startPosi  = trmi.startPos
        endPosi    = startPosi + trmi.nLevels - 1
        samples4locations = sol[startPosi:endPosi]
        writedlm(outfile[trmStr],samples4locations',',')
    end
end


################################################################################
#Output MCMC Samples every output_samples_frequency steps
################################################################################
function output_MCMC_samples(mme,sol,vRes,G0,
                             π=false,
                             α=false,
                             locusEffectVar=false,
                             outfile=false)
  ntraits     = size(mme.lhsVec,1)
  #location parameters
  output_location_parameters_samples(mme,sol,outfile)
  #random effects variances
  for effect in  mme.rndTrmVec
    trmStri   = split(effect.term_array[1],':')[end]
    writedlm(outfile[trmStri*"_variances"],vec(inv(effect.Gi))',',')
  end

  writedlm(outfile["residual_variance"],(typeof(vRes) <: Number) ? vRes : vec(vRes)' ,',')

  if mme.pedTrmVec != 0
    writedlm(outfile["polygenic_effects_variance"],vec(G0)',',')
  end
  if mme.M != 0 && outfile != false
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
      if !(typeof(π) <: Number) #add a blank line
          println(outfile["pi"])
      end
  end

  if mme.MCMCinfo.outputEBV == true
      if mme.output_ID != 0 &&  (mme.pedTrmVec != 0 || mme.M != 0 )
          if ntraits == 1
             myEBV = getEBV(mme,sol,α,1)
             writedlm(outfile["EBV_"*string(mme.lhsVec[1])],myEBV',',')
             if mme.MCMCinfo.output_heritability == true && mme.MCMCinfo.single_step_analysis == false
                 mygvar = var(myEBV)
                 writedlm(outfile["genetic_variance"],mygvar,',')
                 writedlm(outfile["heritability"],mygvar/(mygvar+vRes),',')
             end
          else
              EBVmat = myEBV = getEBV(mme,sol,α[1],1)
              writedlm(outfile["EBV_"*string(mme.lhsVec[1])],myEBV',',')
              for traiti in 2:ntraits
                  myEBV = getEBV(mme,sol,α[traiti],traiti)
                  writedlm(outfile["EBV_"*string(mme.lhsVec[traiti])],myEBV',',')
                  EBVmat = [EBVmat myEBV]
              end
              if mme.MCMCinfo.output_heritability == true && mme.MCMCinfo.single_step_analysis == false
                  mygvar = cov(EBVmat)
                  writedlm(outfile["genetic_variance"],vec(mygvar)',',')
                  writedlm(outfile["heritability"],(diag(mygvar./(mygvar+vRes)))',',')
              end
          end
       end
  end
end

# for vec::Array{AbstractString,1} and vec::Array{String,1}
# and Array{SubString{String} for issue 8
function transubstrarr(vec)
    lvec=length(vec)
    res =Array{String}(undef,1,lvec)
    for i in 1:lvec
        res[1,i]=vec[i]
    end
    return res
end
