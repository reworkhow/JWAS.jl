#This file is used for functions to output, including results (posterior mean and
#varainces) and MCMC samples.
#1.public function
#2.return point estimates as returned values (a dictionary) from runMCMC
#3.return MCMC samples as text files
################################################################################
#*******************************************************************************
#1. PUBLIC FUNCTIONS
#*******************************************************************************
################################################################################
"""
    outputEBV(model,IDs::Array)

Output estimated breeding values and prediction error variances for IDs.
"""
function outputEBV(model,IDs)
    #print("Estimated breeding values and prediction error variances will be included in the output.")
    IDs = map(string,vec(IDs)) #Array{String,1}
    model.output_ID=IDs
end

"""
    outputMCMCsamples(mme::MME,trmStrs::AbstractString...)

Get MCMC samples for specific location parameters.
"""
function outputMCMCsamples(mme::MME,trmStrs::AbstractString...)
    for trmStr in trmStrs
      res = []
      #add trait name to variables,e.g, age => "y1:age"
      #"age" may be in trait 1 but not trait 2
      for (m,model) = enumerate(mme.modelVec)
          strVec  = split(model,['=','+'])
          strpVec = [strip(i) for i in strVec]
          if trmStr in strpVec || trmStr in ["J","ϵ"]
              res = [res;string(mme.lhsVec[m])*":"*trmStr]
          end
      end
      for trmStr in res
          trm     = mme.modelTermDict[trmStr]
          push!(mme.outputSamplesVec,trm)
      end
    end
end

################################################################################
#*******************************************************************************
#2. Return Output Results (Dictionary)
#*******************************************************************************
#Posterior means and variances are calculated for all parameters in the model
#when MCMC is running; Other paramters (e.g., EBV), which is a function of those
#are calculated from files storing MCMC samples at the end of MCMC.
################################################################################
function output_result(mme,output_file,
                       solMean,meanVare,
                       G0Mean,meanAlpha,meanDelta,meanVara,estimatePi,mean_pi,estimateScale,meanScaleVar,
                       solMean2      = missing,
                       meanVare2     = missing,
                       G0Mean2       = missing,
                       meanAlpha2    = mme.nModels == 1 ? missing : fill(missing,mme.nModels),
                       meanVara2     = missing,
                       mean_pi2      = missing,
                       meanScaleVar2 = missing)
  output = Dict()
  location_parameters = reformat2dataframe([getNames(mme) solMean sqrt.(solMean2 .- solMean .^2)])
  output["location parameters"] = location_parameters
  output["residual variance"]   = matrix2dataframe(string.(mme.lhsVec),meanVare,meanVare2)


  if mme.pedTrmVec != 0
    output["polygenic effects covariance matrix"]=matrix2dataframe(mme.pedTrmVec,G0Mean,G0Mean2)
  end

  if mme.M != 0
    if mme.nModels == 1
        meanAlpha =[meanAlpha]
        meanAlpha2=[meanAlpha2]
        meanDelta = [meanDelta]
    end
    traiti      = 1
    whichtrait  = fill(string(mme.lhsVec[traiti]),length(mme.M.markerID))
    whichmarker = mme.M.markerID
    whicheffect = meanAlpha[traiti]
    whicheffectsd = sqrt.(meanAlpha2[traiti] .- meanAlpha[traiti] .^2)
    whichdelta    = meanDelta[traiti]
    for traiti in 2:mme.nModels
        whichtrait     = vcat(whichtrait,fill(string(mme.lhsVec[traiti]),length(mme.M.markerID)))
        whichmarker    = vcat(whichmarker,mme.M.markerID)
        whicheffect    = vcat(whicheffect,meanAlpha[traiti])
        whicheffectsd  = vcat(whicheffectsd,sqrt.(meanAlpha2[traiti] .- meanAlpha[traiti] .^2))
        whichdelta     = vcat(whichdelta,meanDelta[traiti])
    end
    output["marker effects"]=
    DataFrame([whichtrait whichmarker whicheffect whicheffectsd whichdelta],[:Trait,:Marker_ID,:Estimate,:Std_Error,:Model_Frequency])

    output["marker effects variance"] = matrix2dataframe(string.(mme.lhsVec),meanVara,meanVara2)
    if estimatePi == true
        output["Pi"] = dict2dataframe(mean_pi,mean_pi2)
    end
    if estimateScale == true
        output["ScaleEffectVar"] = matrix2dataframe(string.(mme.lhsVec),meanScaleVar,meanScaleVar2)
    end
  end
  #Get EBV and PEV from MCMC samples text files
  if mme.output_ID != 0 && mme.MCMCinfo.outputEBV == true
      for traiti in 1:mme.nModels
          EBVkey         = "EBV"*"_"*string(mme.lhsVec[traiti])
          EBVsamplesfile = output_file*"_"*EBVkey*".txt"
          EBVsamples,IDs = readdlm(EBVsamplesfile,',',header=true)
          EBV            = vec(mean(EBVsamples,dims=1))
          PEV            = vec(var(EBVsamples,dims=1))
          if vec(IDs) == mme.output_ID
              output[EBVkey] = DataFrame([mme.output_ID EBV PEV],[:ID,:EBV,:PEV])
          else
              error("The EBV file is wrong.")
          end
      end
      if mme.MCMCinfo.output_heritability == true  && mme.MCMCinfo.single_step_analysis == false
          for i in ["genetic_variance","heritability"]
              samplesfile = output_file*"_"*i*".txt"
              samples,names = readdlm(samplesfile,',',header=true)
              samplemean    = vec(mean(samples,dims=1))
              samplevar     = vec(std(samples,dims=1))
              output[i] = DataFrame([vec(names) samplemean samplevar],[:Covariance,:Estimate,:Std_Error])
          end
      end

  end
  return output
end

#Reformat Output Array to DataFrame
function reformat2dataframe(res::Array)
    out_names = Array{String}(undef,size(res,1),3)
    for rowi in 1:size(res,1)
        out_names[rowi,:]=[strip(i) for i in split(res[rowi,1],':',keepempty=false)]
    end

    if size(out_names,2)==1 #convert vector to matrix
        out_names = reshape(out_names,length(out_names),1)
    end
    #out_names=permutedims(out_names,[2,1]) #rotate
    out_values   = map(Float64,res[:,2])
    out_variance = convert.(Union{Missing, Float64},res[:,3])
    out =[out_names out_values out_variance]
    out = DataFrame(out, [:Trait, :Effect, :Level, :Estimate,:Std_Error])
    return out
end

function matrix2dataframe(names,meanVare,meanVare2)
    names     = repeat(names,inner=length(names)).*"_".*repeat(names,outer=length(names))
    meanVare  = (typeof(meanVare) <: Union{Number,Missing}) ? meanVare : vec(meanVare)
    meanVare2 = (typeof(meanVare2) <: Union{Number,Missing}) ? meanVare2 : vec(meanVare2)
    stdVare   = sqrt.(meanVare2 .- meanVare .^2)
    DataFrame([names meanVare stdVare],[:Covariance,:Estimate,:Std_Error])
end

function dict2dataframe(mean_pi,mean_pi2)
    if typeof(mean_pi) <: Union{Number,Missing}
        names = "π"
    else
        names = collect(keys(mean_pi))
    end
    mean_pi  = (typeof(mean_pi) <: Union{Number,Missing}) ? mean_pi : collect(values(mean_pi))
    mean_pi2 = (typeof(mean_pi2) <: Union{Number,Missing}) ? mean_pi2 : collect(values(mean_pi2))
    stdpi    = sqrt.(mean_pi2 .- mean_pi .^2)
    DataFrame([names mean_pi stdpi],[:π,:Estimate,:Std_Error])
end

"""
    getEBV(model::MME,sol,α,traiti)

(internal function) Get breeding values for individuals defined by outputEBV(),
defaulting to all genotyped individuals. This function is used inside MCMC functions for
one MCMC samples from posterior distributions.
"""
function getEBV(mme,sol,α,traiti)
    traiti = string(mme.lhsVec[traiti])
    EBV=zeros(length(mme.output_ID))

    location_parameters = reformat2dataframe([getNames(mme) sol zero(sol)])
    if mme.pedTrmVec != 0
        for pedtrm in mme.pedTrmVec
            mytrait, effect = split(pedtrm,':')
            if mytrait == traiti
                sol_pedtrm     = map(Float64,location_parameters[(location_parameters[!,:Effect].==effect).&(location_parameters[!,:Trait].==traiti),:Estimate])
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
            sol_J  = map(Float64,location_parameters[(location_parameters[!,:Effect].=="J").&(location_parameters[!,:Trait].==string(mme.lhsVec[traiti])),:Estimate])[1]
            sol_ϵ  = map(Float64,location_parameters[(location_parameters[!,:Effect].=="ϵ").&(location_parameters[!,:Trait].==string(mme.lhsVec[traiti])),:Estimate])
            EBV_J  = mme.output_X["J"]*sol_J
            EBV_ϵ  = mme.output_X["ϵ"]*sol_ϵ
            EBV   += (EBV_J+EBV_ϵ)
        end
    end

    return EBV
end

"""
    get_outputX_others(model)

(internal function) Make incidence matrices for effects involved in
calculation of EBV including J, ϵ, pedTrmVec except marker covariates.
"""
function get_outputX_others(model,single_step_analysis)
    #trick to avoid errors (PedModule.getIDs(ped) [nongeno ID;geno ID])
    if single_step_analysis == true
        model.output_X["ϵ"]=mkmat_incidence_factor(model.output_ID,
                            PedModule.getIDs(model.ped))[:,1:length(model.ped.setNG)]
        #Note that model.output_X["J"] is in SSBRrun
    end
    if model.pedTrmVec != 0
        for i in model.pedTrmVec
            model.output_X[i]=mkmat_incidence_factor(model.output_ID,
                                                     model.modelTermDict[i].names)
        end
    end
end
################################################################################
#*******************************************************************************
#3. Save MCMC Samples to Text Files
#*******************************************************************************
#MCMC samples for all hyperparameters, all marker effects, and
#location parameters defined by outputMCMCsamples() are saved
#into files every output_samples_frequency iterations.
################################################################################
"""
    output_MCMC_samples_setup(mme,nIter,output_samples_frequency,file_name="MCMC_samples")

(internal function) Set up text files to save MCMC samples.
"""
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

  if mme.M !=0
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
"""
    output_MCMC_samples(mme,sol,vRes,G0,π=false,α=false,locusEffectVar=false,outfile=false)

(internal function) Save MCMC samples every output_samples_frequency iterations to the text file.
"""
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

  if mme.MCMCinfo.outputEBV == true #add error message
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
                  myEBV = getEBV(mme,sol,α[traiti],traiti) #actually BV
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
"""
    output_location_parameters_samples(mme::MME,sol,outfile)

(internal function) Save MCMC samples for location parameers
"""
function output_location_parameters_samples(mme::MME,sol,outfile)
    for trmi in  mme.outputSamplesVec
        trmStr = trmi.trmStr
        startPosi  = trmi.startPos
        endPosi    = startPosi + trmi.nLevels - 1
        samples4locations = sol[startPosi:endPosi]
        writedlm(outfile[trmStr],samples4locations',',')
    end
end
"""
    transubstrarr(vec)

(internal function) Transpose a column vector of strings (vec' doesn't work here)
"""
function transubstrarr(vec)
    lvec=length(vec)
    res =Array{String}(undef,1,lvec)
    for i in 1:lvec
        res[1,i]=vec[i]
    end
    return res
end
