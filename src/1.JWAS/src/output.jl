#This file is used for output, including results (posterior mean and varainces)
#and MCMC samples.
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
    IDs = map(string,vec(IDs))
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
function output_result(mme,output_folder,
                       solMean,meanVare,G0Mean,
                       solMean2 = missing,meanVare2 = missing,G0Mean2 = missing)
  output = Dict()
  location_parameters = reformat2dataframe([getNames(mme) solMean sqrt.(abs.(solMean2 .- solMean .^2))])
  output["location parameters"] = location_parameters
  output["residual variance"]   = matrix2dataframe(string.(mme.lhsVec),meanVare,meanVare2)


  if mme.pedTrmVec != 0
    output["polygenic effects covariance matrix"]=matrix2dataframe(mme.pedTrmVec,G0Mean,G0Mean2)
  end

  if mme.M != 0
      for Mi in mme.M
          traiti      = 1
          whichtrait  = fill(string(mme.lhsVec[traiti]),length(Mi.markerID))
          whichmarker = Mi.markerID
          whicheffect = Mi.meanAlpha[traiti]
          whicheffectsd = sqrt.(abs.(Mi.meanAlpha2[traiti] .- Mi.meanAlpha[traiti] .^2))
          whichdelta    = Mi.meanDelta[traiti]
          for traiti in 2:mme.nModels
              whichtrait     = vcat(whichtrait,fill(string(mme.lhsVec[traiti]),length(Mi.markerID)))
              whichmarker    = vcat(whichmarker,Mi.markerID)
              whicheffect    = vcat(whicheffect,Mi.meanAlpha[traiti])
              whicheffectsd  = vcat(whicheffectsd,sqrt.(abs.(Mi.meanAlpha2[traiti] .- Mi.meanAlpha[traiti] .^2)))
              whichdelta     = vcat(whichdelta,Mi.meanDelta[traiti])
          end
          output["marker effects "*Mi.name]=DataFrame([whichtrait whichmarker whicheffect whicheffectsd whichdelta],[:Trait,:Marker_ID,:Estimate,:Std_Error,:Model_Frequency])
          #output["marker effects variance "*Mi.name] = matrix2dataframe(string.(mme.lhsVec),Mi.meanVara,Mi.meanVara2)
          if Mi.estimatePi == true
              output["pi_"*Mi.name] = dict2dataframe(Mi.mean_pi,Mi.mean_pi2)
          end
          if Mi.estimateScale == true
              output["ScaleEffectVar"*Mi.name] = matrix2dataframe(string.(mme.lhsVec),Mi.meanScaleVara,Mi.meanScaleVara2)
          end
      end
  end
  #Get EBV and PEV from MCMC samples text files
  if mme.output_ID != 0 && mme.MCMCinfo.outputEBV == true
      output_file = output_folder*"/MCMC_samples"
      EBVkeys = ["EBV"*"_"*string(mme.lhsVec[traiti]) for traiti in 1:mme.nModels]
      if mme.latent_traits == true
          push!(EBVkeys, "EBV_NonLinear")
      end
      for EBVkey in EBVkeys
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

      if mme.latent_traits == true && mme.nonlinear_function == "Neural Network"
          myvar         = "neural_networks_bias_and_weights"
          samplesfile   = output_file*"_"*myvar*".txt"
          samples       = readdlm(samplesfile,',',header=false)
          names         = ["bias";"weight".*string.(1:(size(samples,2)-1))]
          samplemean    = vec(mean(samples,dims=1))
          samplevar     = vec(std(samples,dims=1))
          output[myvar] = DataFrame([vec(names) samplemean samplevar],[:weights,:Estimate,:Std_Error])
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
    stdVare   = sqrt.(abs.(meanVare2 .- meanVare .^2))
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
    stdpi    = sqrt.(abs.(mean_pi2 .- mean_pi .^2))
    DataFrame([names mean_pi stdpi],[:π,:Estimate,:Std_Error])
end

"""
    getEBV(model::MME,traiti)

(internal function) Get breeding values for individuals defined by outputEBV(),
defaulting to all genotyped individuals. This function is used inside MCMC functions for
one MCMC samples from posterior distributions.
"""
function getEBV(mme,traiti)
    traiti_name = string(mme.lhsVec[traiti])
    EBV=zeros(length(mme.output_ID))

    location_parameters = reformat2dataframe([getNames(mme) mme.sol zero(mme.sol)])
    if mme.pedTrmVec != 0
        for pedtrm in mme.pedTrmVec
            mytrait, effect = split(pedtrm,':')
            if mytrait == traiti_name
                sol_pedtrm     = map(Float64,location_parameters[(location_parameters[!,:Effect].==effect).&(location_parameters[!,:Trait].==traiti_name),:Estimate])
                EBV_pedtrm     = mme.output_X[pedtrm]*sol_pedtrm
                EBV += EBV_pedtrm
            end
        end
    end

    if mme.M != 0
        for Mi in mme.M
            EBV += Mi.output_genotypes*Mi.α[traiti]
        end
    end

    if haskey(mme.output_X,"J") #single-step analyis
            sol_J  = map(Float64,location_parameters[(location_parameters[!,:Effect].=="J").&(location_parameters[!,:Trait].==traiti_name),:Estimate])[1]
            EBV_J  = mme.output_X["J"]*sol_J
            EBV   += EBV_J
    end
    if haskey(mme.output_X,"ϵ") #single-step analyis
            sol_ϵ  = map(Float64,location_parameters[(location_parameters[!,:Effect].=="ϵ").&(location_parameters[!,:Trait].==traiti_name),:Estimate])
            EBV_ϵ  = mme.output_X["ϵ"]*sol_ϵ
            EBV   += EBV_ϵ
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
      for Mi in mme.M
          for traiti in 1:ntraits
              push!(outvar,"marker_effects_"*Mi.name*"_"*string(mme.lhsVec[traiti]))
          end
          push!(outvar,"marker_effects_variances"*"_"*Mi.name)
          push!(outvar,"pi"*"_"*Mi.name)
      end
  end

  #non-marker location parameters
  for trmi in  mme.outputSamplesVec
      trmStr = trmi.trmStr
      push!(outvar,trmStr)
  end

  #non-marker random effects variances
  for i in  mme.rndTrmVec
      trmStri   = join(i.term_array, "_")                  #x2
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
      if mme.latent_traits == true
          push!(outvar,"EBV_NonLinear")
          if mme.nonlinear_function == "Neural Network"
              push!(outvar,"neural_networks_bias_and_weights")
          end
      end
  end
  #categorical traits
  if mme.MCMCinfo.categorical_trait == true
      push!(outvar,"threshold")
  end

  for i in outvar
      file_i    = file_name*"_"*replace(i,":"=>".")*".txt" #replace ":" by "." to avoid reserved characters in Windows
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
    trmStri   = join(effect.term_array, "_")                  #x2
    thisheader= repeat(effect.term_array,inner=length(effect.term_array)).*"_".*repeat(effect.term_array,outer=length(effect.term_array))
    writedlm(outfile[trmStri*"_variances"],transubstrarr(thisheader),',') #1:x2_1:x2,1:x2_2:x2,2:x2_1:x2,2:x2_2:x2
  end

  if mme.M !=0
      for Mi in mme.M
          for traiti in 1:ntraits
              writedlm(outfile["marker_effects_"*Mi.name*"_"*string(mme.lhsVec[traiti])],transubstrarr(Mi.markerID),',')
          end
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
      if mme.latent_traits == true
          writedlm(outfile["EBV_NonLinear"],transubstrarr(mme.output_ID),',')
      end
  end

  return outfile
end
"""
    output_MCMC_samples(mme,vRes,G0,outfile=false)

(internal function) Save MCMC samples every output_samples_frequency iterations to the text file.
"""
function output_MCMC_samples(mme,vRes,G0,
                             outfile=false)
    ntraits     = size(mme.lhsVec,1)
    #location parameters
    output_location_parameters_samples(mme,mme.sol,outfile)
    #random effects variances
    for effect in  mme.rndTrmVec
    trmStri   = join(effect.term_array, "_")
    writedlm(outfile[trmStri*"_variances"],vec(inv(effect.Gi))',',')
    end

    writedlm(outfile["residual_variance"],(typeof(vRes) <: Number) ? vRes : vec(vRes)' ,',')

    if mme.pedTrmVec != 0
    writedlm(outfile["polygenic_effects_variance"],vec(G0)',',')
    end
    if mme.M != 0 && outfile != false
      for Mi in mme.M
          for traiti in 1:ntraits
              writedlm(outfile["marker_effects_"*Mi.name*"_"*string(mme.lhsVec[traiti])],Mi.α[traiti]',',')
          end
          if Mi.G != false
              if mme.nModels == 1
                  writedlm(outfile["marker_effects_variances"*"_"*Mi.name],Mi.G',',')
              else
                  if methods in ["BayesC","RR-BLUP","BayesL","GBLUP"]
                      writedlm(outfile["marker_effects_variances"*"_"*Mi.name],vec(Mi.G),',')
                  elseif methods == "BayesB"
                      writedlm(outfile["marker_effects_variances"*"_"*Mi.name],hcat([x for x in Mi.G]...),',')
                  end
              end
          end
          writedlm(outfile["pi"*"_"*Mi.name],Mi.π,',')
          if !(typeof(Mi.π) <: Number) #add a blank line
              println(outfile["pi"*"_"*Mi.name])
          end
      end
    end

    if mme.MCMCinfo.outputEBV == true #add error message
      if mme.output_ID != 0 &&  (mme.pedTrmVec != 0 || mme.M != 0 )
          if ntraits == 1
             EBVmat = myEBV = getEBV(mme,1)
             writedlm(outfile["EBV_"*string(mme.lhsVec[1])],myEBV',',')
             if mme.MCMCinfo.output_heritability == true && mme.MCMCinfo.single_step_analysis == false
                 mygvar = var(myEBV)
                 writedlm(outfile["genetic_variance"],mygvar,',')
                 writedlm(outfile["heritability"],mygvar/(mygvar+vRes),',')
             end
          else
              EBVmat = myEBV = getEBV(mme,1)
              writedlm(outfile["EBV_"*string(mme.lhsVec[1])],myEBV',',')
              for traiti in 2:ntraits
                  myEBV = getEBV(mme,traiti) #actually BV
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
    if mme.latent_traits == true
        EBVmat = EBVmat .+ mme.sol' #mme.sol here only contains intercepts
        if mme.nonlinear_function != "Neural Network"
            BV_NN = mme.nonlinear_function.(Tuple([view(EBVmat,:,i) for i in 1:size(EBVmat,2)])...)
        else
            BV_NN = [ones(size(EBVmat,1)) tanh.(EBVmat)]*mme.weights_NN
            writedlm(outfile["neural_networks_bias_and_weights"],mme.weights_NN',',')
        end
        writedlm(outfile["EBV_NonLinear"],BV_NN',',')
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


#output mean and variance of posterior distribution of parameters of interest
function output_posterior_mean_variance(mme,nsamples)
    mme.solMean   += (mme.sol - mme.solMean)/nsamples
    mme.solMean2  += (mme.sol .^2 - mme.solMean2)/nsamples
    mme.meanVare  += (mme.R - mme.meanVare)/nsamples
    mme.meanVare2 += (mme.R .^2 - mme.meanVare2)/nsamples

    if mme.pedTrmVec != 0
        mme.G0Mean  += (inv(mme.Gi)  - mme.G0Mean )/nsamples
        mme.G0Mean2 += (inv(mme.Gi) .^2  - mme.G0Mean2 )/nsamples
    end
    if mme.M != 0
        for Mi in mme.M
            for trait in 1:Mi.ntraits
                Mi.meanAlpha[trait] += (Mi.α[trait] - Mi.meanAlpha[trait])/nsamples
                Mi.meanAlpha2[trait]+= (Mi.α[trait].^2 - Mi.meanAlpha2[trait])/nsamples
                Mi.meanDelta[trait] += (Mi.δ[trait] - Mi.meanDelta[trait])/nsamples
            end
            if Mi.estimatePi == true
                if Mi.ntraits == 1 || mme.MCMCinfo.mega_trait
                    Mi.mean_pi += (Mi.π-Mi.mean_pi)/nsamples
                    Mi.mean_pi2 += (Mi.π .^2-Mi.mean_pi2)/nsamples
                else
                    for i in keys(Mi.π)
                      Mi.mean_pi[i] += (Mi.π[i]-Mi.mean_pi[i])/nsamples
                      Mi.mean_pi2[i] += (Mi.π[i].^2-Mi.mean_pi2[i])/nsamples
                    end
                end
            end
            if Mi.method != "BayesB"
                Mi.meanVara += (Mi.G - Mi.meanVara)/nsamples
                Mi.meanVara2 += (Mi.G .^2 - Mi.meanVara2)/nsamples
            end
            if Mi.estimateScale == true
                Mi.meanScaleVara += (Mi.scale - Mi.meanScaleVara)/nsamples
                Mi.meanScaleVara2 += (Mi.scale .^2 - Mi.meanScaleVara2)/nsamples
            end
        end
    end
end
