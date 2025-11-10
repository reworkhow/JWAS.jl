# ===============================================================================================================
# This script handles the output of analysis, including results (posterior mean and variances) and MCMC samples.
# Features:
# - contains public functions
# - returns point estimates and standard deveiations (as a dictionary) from runMCMC function
# - outputs MCMC samples in text file format
# ===============================================================================================================

# =============================================================================
#                                KEY TERMS & NOTATIONS
# =============================================================================
# 
# PEV: Prediction Error Variance of a random effect, PEV=var(u-u^hat)=var(u|y), 
#      i.e., the variance of the posterior distribution of u. 
#
#

# ===============================================================================================================
#                                PUBLIC FUNCTIONS
# ===============================================================================================================
"""
    prediction_setup(mme::MME)

* (internal function) Create incidence matrices for individuals of interest based on a usere-defined
prediction equation, defaulting to genetic values including effects defined with genomic and pedigre
information. For now, genomic data is always included.
* J and ϵ are always included in single-step analysis (added in SSBR.jl)
"""
function prediction_setup(model)
    if model.MCMCinfo.prediction_equation == false
        prediction_equation = []
        if model.pedTrmVec != 0
            for i in model.pedTrmVec
                push!(prediction_equation,i)
            end
        end
    else
        prediction_equation = string.(strip.(split(model.MCMCinfo.prediction_equation,"+")))
        if model.MCMCinfo.output_heritability != false
            printstyled("User-defined prediction equation is provided. ","The heritability is the ",
            "proportion of phenotypic variance explained by the value defined by the prediction equation.\n",
            bold=false,color=:green)
        end
        for i in prediction_equation
            term_symbol = Symbol(split(i,":")[end])
            if !(haskey(model.modelTermDict,i) || (isdefined(Main,term_symbol) && typeof(getfield(Main,term_symbol)) == Genotypes))
                error("Terms $i in the prediction equation is not found.")
            end
        end
    end
    printstyled("Predicted values for individuals of interest will be obtained as the summation of ",
    prediction_equation, " (Note that genomic data is always included for now).",bold=false,color=:green)
    if length(prediction_equation) == 0 && model.M == false
        println("Default or user-defined prediction equation are not available.")
        model.MCMCinfo.outputEBV = false
    end
    filter!(e->(e in keys(model.modelTermDict)),prediction_equation) #remove "genotypes" for now
    model.MCMCinfo.prediction_equation = prediction_equation
end

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


# ===============================================================================================================
#                                PRIVATE FUNCTIONS
# ===============================================================================================================

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
  if mme.MCMCinfo.RRM == false
      output["residual variance"]   = matrix2dataframe(string.(mme.lhsVec),meanVare,meanVare2)
  else
      output["residual variance"]   = matrix2dataframe(["1"],meanVare,meanVare2)
  end

  if mme.pedTrmVec != 0
    output["polygenic effects covariance matrix"]=matrix2dataframe(mme.pedTrmVec,G0Mean,G0Mean2) 
  end

  ntraits = length(mme.lhsVec)

  if mme.M != 0
      for Mi in mme.M
         ntraits_geno = mme.MCMCinfo.RRM == false ? Mi.ntraits : length(mme.lhsVec)
         traiti      = 1
         whichtrait  = fill(string(mme.lhsVec[traiti]),length(Mi.markerID))
         whichmarker = Mi.markerID
         whicheffect = Mi.meanAlpha[traiti]
         whicheffectsd = sqrt.(abs.(Mi.meanAlpha2[traiti] .- Mi.meanAlpha[traiti] .^2))
         whichdelta    = Mi.meanDelta[traiti]
          for traiti in 2:ntraits_geno
                whichtrait     = vcat(whichtrait,fill(string(mme.lhsVec[traiti]),length(Mi.markerID)))
                whichmarker    = vcat(whichmarker,Mi.markerID)
                whicheffect    = vcat(whicheffect,Mi.meanAlpha[traiti])
                whicheffectsd  = vcat(whicheffectsd,sqrt.(abs.(Mi.meanAlpha2[traiti] .- Mi.meanAlpha[traiti] .^2)))
                whichdelta     = vcat(whichdelta,Mi.meanDelta[traiti])
            end

          output["marker effects "*Mi.name]=DataFrame([whichtrait whichmarker whicheffect whicheffectsd whichdelta],[:Trait,:Marker_ID,:Estimate,:SD,:Model_Frequency])
          #output["marker effects variance "*Mi.name] = matrix2dataframe(string.(mme.lhsVec),Mi.meanVara,Mi.meanVara2)
          if Mi.estimatePi == true
              output["pi_"*Mi.name] = dict2dataframe(Mi.mean_pi,Mi.mean_pi2)
          end
          if Mi.G.estimate_scale == true
              output["ScaleEffectVar"*Mi.name] = matrix2dataframe(string.(mme.lhsVec),Mi.meanScaleVara,Mi.meanScaleVara2)
          end
      end
  end

  #Get EBV and PEV from MCMC samples text files
  if mme.output_ID != 0 && mme.MCMCinfo.outputEBV == true
      output_file = output_folder*"/MCMC_samples"
      EBVkeys = ["EBV"*"_"*string(mme.lhsVec[traiti]) for traiti in 1:ntraits]
      if mme.nonlinear_function != false  #NNBayes
          push!(EBVkeys, "EBV_NonLinear")
          EBVkeys=[EBVkeys[end]]  #only keep "EBV_NonLinear" (remove EBV_gene1, EBV_gene2,...)
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
          if mme.MCMCinfo.RRM != false
              genetic_trm = ["genetic_variance"]
          else
              genetic_trm = ["genetic_variance","heritability"]
          end
          for i in genetic_trm
              samplesfile = output_file*"_"*i*".txt"
              samples,names = readdlm(samplesfile,',',header=true)
              samplemean    = vec(mean(samples,dims=1))
              samplevar     = vec(std(samples,dims=1))
              output[i] = DataFrame([vec(names) samplemean samplevar],[:Covariance,:Estimate,:SD])
          end
      end
      if mme.nonlinear_function != false && mme.is_activation_fcn == true  #Neural Network with activation function
          myvar         = "neural_networks_bias_and_weights"
          samplesfile   = output_file*"_"*myvar*".txt"
          samples       = readdlm(samplesfile,',',header=false)
          names         = ["bias";"weight".*string.(1:(size(samples,2)-1))]
          samplemean    = vec(mean(samples,dims=1))
          samplevar     = vec(std(samples,dims=1))
          output[myvar] = DataFrame([vec(names) samplemean samplevar],[:weights,:Estimate,:SD])
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
    out = DataFrame(out, [:Trait, :Effect, :Level, :Estimate,:SD])
    return out
end

#convert a scalar (single-trait), a matrix (multi-trait), a vector (mega-trait) to a DataFrame
function matrix2dataframe(names,meanVare,meanVare2) #also works for scalar
    if !(typeof(meanVare) <: Vector)
        names = repeat(names,inner=length(names)).*"_".*repeat(names,outer=length(names))
    end
    meanVare  = (typeof(meanVare)  <: Union{Number,Missing,Vector}) ? meanVare  : vec(meanVare)
    meanVare2 = (typeof(meanVare2) <: Union{Number,Missing,Vector}) ? meanVare2 : vec(meanVare2)
    stdVare   = sqrt.(abs.(meanVare2 .- meanVare .^2))
    DataFrame([names meanVare stdVare],[:Covariance,:Estimate,:SD])
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
    DataFrame([names mean_pi stdpi],[:π,:Estimate,:SD])
end

"""
    getEBV(model::MME,traiti)

(internal function) Get breeding values for individuals defined by outputEBV(),
defaulting to all genotyped individuals. This function is used inside MCMC functions for
one MCMC samples from posterior distributions.
e.g.,
non-NNBayes_partial (multi-classs Bayes) : y1=M1*α1[1]+M2*α2[1]+M3*α3[1]
                                           y2=M1*α1[2]+M2*α2[2]+M3*α3[2];
NNBayes_partial:     y1=M1*α1[1]
                     y2=M2*α2[1]
                     y3=M3*α3[1];
"""
function getEBV(mme,traiti)
    traiti_name = string(mme.lhsVec[traiti])
    EBV=zeros(length(mme.output_ID))

    location_parameters = reformat2dataframe([getNames(mme) mme.sol zero(mme.sol)])
    for term in keys(mme.output_X)
        mytrait, effect = split(term,':')
        if mytrait == traiti_name
            sol_term     = map(Float64,location_parameters[(location_parameters[!,:Effect].==effect).&(location_parameters[!,:Trait].==traiti_name),:Estimate])
            if length(sol_term) == 1 #1-element Array{Float64,1} doesn't work below; Will be deleted
                sol_term = sol_term[1]
            end
            EBV_term = mme.output_X[term]*sol_term
            if length(sol_term) == 1 #1-element Array{Float64,1} doesn't work below; Will be deleted
                EBV_term = vec(EBV_term)
            end
            EBV += EBV_term
        end
    end
    is_partial_connect = mme.nonlinear_function != false && mme.is_fully_connected==false
    if mme.M != 0
        for i in 1:length(mme.M) #for Mi in mme.M
            Mi=mme.M[i]
            if !is_partial_connect  #non-NNBayes_partial
                EBV += Mi.output_genotypes*Mi.α[traiti]
            else  #NNBayes_partial
                if i==traiti
                    EBV = Mi.output_genotypes*mme.M[i].α[1]
                end
            end
        end
    end
    return EBV
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
        geno_names = mme.MCMCinfo.RRM == false ? Mi.trait_names : string.(mme.lhsVec)
        for traiti in geno_names
            push!(outvar,"marker_effects_"*Mi.name*"_"*traiti)
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
          if mme.MCMCinfo.RRM == false
              push!(outvar,"heritability")
          else
              printstyled("heritability is not computed for Random Regression Model. \n",bold=false,color=:green)
          end
      end
      if mme.nonlinear_function != false  #NNBayes
          push!(outvar,"EBV_NonLinear")
          if mme.is_activation_fcn == true #Neural Network with activation function
              push!(outvar,"neural_networks_bias_and_weights")
          end
      end
  end
  #categorical/censored traits
  for t in 1:mme.nModels
      if mme.traits_type[t] ∈ ["categorical","categorical(binary)","censored"] #liability are sampled
          push!(outvar,"liabilities_"*string(mme.lhsVec[t]))
          if mme.traits_type[t] ∈ ["categorical","categorical(binary)"] #thresholds will be saved
              push!(outvar,"threshold_"*string(mme.lhsVec[t]))
          end
      end
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
  if mme.R.constraint == false
      varheader = repeat(mytraits,inner=length(mytraits)).*"_".*repeat(mytraits,outer=length(mytraits))
  else
      varheader = transubstrarr(map(string,mme.lhsVec))
  end
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
        geno_names = mme.MCMCinfo.RRM == false ? Mi.trait_names : string.(mme.lhsVec)
        for traiti in geno_names
            writedlm(outfile["marker_effects_"*Mi.name*"_"*traiti],transubstrarr(Mi.markerID),',')
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
          varheader = repeat(mytraits,inner=length(mytraits)).*"_".*repeat(mytraits,outer=length(mytraits))

          writedlm(outfile["genetic_variance"],transubstrarr(varheader),',')
          if mme.MCMCinfo.RRM == false
              writedlm(outfile["heritability"],transubstrarr(map(string,mme.lhsVec)),',')
          end
      end
      if mme.nonlinear_function != false #NNBayes
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
        writedlm(outfile[trmStri*"_variances"],vec(inv(effect.Gi.val))',',')
    end

    if mme.R.constraint == true
        vRes=diag(vRes)
    end
    writedlm(outfile["residual_variance"],(typeof(vRes) <: Number) ? vRes : vec(vRes)' ,',')

    if mme.pedTrmVec != 0
        writedlm(outfile["polygenic_effects_variance"],vec(G0)',',')
    end
    is_partial_connect = mme.nonlinear_function != false && mme.is_fully_connected==false
    if mme.M != 0 && outfile != false
      for Mi in mme.M
         ntraits_geno = mme.MCMCinfo.RRM == false ? Mi.ntraits : length(mme.lhsVec)
         geno_names = mme.MCMCinfo.RRM == false ? Mi.trait_names : string.(mme.lhsVec)
         for traiti in 1:ntraits_geno
            writedlm(outfile["marker_effects_"*Mi.name*"_"*geno_names[traiti]],Mi.α[traiti]',',')
         end
          
         if Mi.G.val != false && mme.nonlinear_function == false #do not save marker effect variances in NNMM to save space
              if mme.nModels == 1
                  writedlm(outfile["marker_effects_variances"*"_"*Mi.name],Mi.G.val',',')
              else
                  if Mi.method == "BayesB"
                      writedlm(outfile["marker_effects_variances"*"_"*Mi.name],hcat([x for x in Mi.G.val]...),',')
                  else
                      writedlm(outfile["marker_effects_variances"*"_"*Mi.name],Mi.G.val,',')
                  end
              end
          end
          writedlm(outfile["pi"*"_"*Mi.name],Mi.π,',')
          if !(typeof(Mi.π) <: Number) #add a blank line
              println(outfile["pi"*"_"*Mi.name])
          end
      end
    end

    if mme.MCMCinfo.outputEBV == true
         EBVmat = myEBV = getEBV(mme,1)
         writedlm(outfile["EBV_"*string(mme.lhsVec[1])],myEBV',',')
         for traiti in 2:ntraits
             myEBV = getEBV(mme,traiti) #actually BV
             trait_name = is_partial_connect ? mme.M[traiti].trait_names[1] : string(mme.lhsVec[traiti])
             writedlm(outfile["EBV_"*trait_name],myEBV',',')
             EBVmat = [EBVmat myEBV]
         end

         if mme.MCMCinfo.output_heritability == true && mme.MCMCinfo.single_step_analysis == false
             #single-trait: a scalar ;  multi-trait: a matrix; mega-trait: a vector
             if mme.M != 0 && mme.M[1].G.constraint==true
                mygvar = Diagonal(vec(var(EBVmat,dims=1)))
             else
                mygvar = cov(EBVmat)
             end
             genetic_variance = (ntraits == 1 ? mygvar : vec(mygvar)')
             if mme.MCMCinfo.RRM == false
                 vRes = mme.R.constraint==true ? Diagonal(vRes) : vRes #change to diagonal matrix to avoid error
                 heritability = (ntraits == 1 ? mygvar/(mygvar+vRes) : (diag(mygvar)./(diag(mygvar)+diag(vRes)))')
                 writedlm(outfile["heritability"],heritability,',')
             end
             writedlm(outfile["genetic_variance"],genetic_variance,',')
         end
    end
    if mme.nonlinear_function != false #NNBayes
        # EBVmat = EBVmat .+ mme.sol' #mme.sol here only contains intercepts
        # do not include intercepts in EBVmat
        if mme.is_activation_fcn == false  #user-defined nonlinear function
            BV_NN = mme.nonlinear_function.(Tuple([view(EBVmat,:,i) for i in 1:size(EBVmat,2)])...)
        else  #activation function
            # BV_NN = [ones(size(EBVmat,1)) mme.nonlinear_function.(EBVmat)]*mme.weights_NN
            BV_NN = mme.nonlinear_function.(EBVmat)*mme.weights_NN #do not include intercepts
            writedlm(outfile["neural_networks_bias_and_weights"],mme.weights_NN',',')
        end
        writedlm(outfile["EBV_NonLinear"],BV_NN',',')
    end
    #categorical/binary/censored traits
    if !isempty(intersect(mme.traits_type, ["categorical","categorical(binary)","censored"]))
        ySparse = reshape(mme.ySparse,:,ntraits) #liability (=mme.ySparse)
        for t in 1:mme.nModels
            if mme.traits_type[t] ∈ ["categorical","categorical(binary)","censored"] #save liability
                writedlm(outfile["liabilities_"*string(mme.lhsVec[t])], ySparse[:,t]', ',')
                if mme.traits_type[t] ∈ ["categorical","categorical(binary)"] #save thresholds
                    writedlm(outfile["threshold_"*string(mme.lhsVec[t])], mme.thresholds[t]', ',')
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

#output mean and variance of posterior distribution of parameters of interest
function output_posterior_mean_variance(mme,nsamples)
    mme.solMean   += (mme.sol - mme.solMean)/nsamples
    mme.solMean2  += (mme.sol .^2 - mme.solMean2)/nsamples
    mme.meanVare  += (mme.R.val - mme.meanVare)/nsamples
    mme.meanVare2 += (mme.R.val .^2 - mme.meanVare2)/nsamples

    if mme.pedTrmVec != 0
        polygenic_pos = findfirst(i -> i.randomType=="A", mme.rndTrmVec)
        mme.G0Mean  += (inv(mme.rndTrmVec[polygenic_pos].Gi.val)  - mme.G0Mean )/nsamples
        mme.G0Mean2 += (inv(mme.rndTrmVec[polygenic_pos].Gi.val) .^2  - mme.G0Mean2 )/nsamples
    end
    if mme.M != 0
        for Mi in mme.M
            for trait in 1:Mi.ntraits
                Mi.meanAlpha[trait] += (Mi.α[trait] - Mi.meanAlpha[trait])/nsamples
                Mi.meanAlpha2[trait]+= (Mi.α[trait].^2 - Mi.meanAlpha2[trait])/nsamples
                Mi.meanDelta[trait] += (Mi.δ[trait] - Mi.meanDelta[trait])/nsamples
            end
            if Mi.estimatePi == true
                if Mi.ntraits == 1 || mme.M[1].G.constraint==true #may need to change for multiple M
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
                Mi.meanVara += (Mi.G.val - Mi.meanVara)/nsamples
                Mi.meanVara2 += (Mi.G.val .^2 - Mi.meanVara2)/nsamples
            end
            if Mi.G.estimate_scale == true
                Mi.meanScaleVara += (Mi.G.scale - Mi.meanScaleVara)/nsamples
                Mi.meanScaleVara2 += (Mi.G.scale .^2 - Mi.meanScaleVara2)/nsamples
            end
        end
    end
end
