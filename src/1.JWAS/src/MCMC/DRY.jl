################################################################################
# Pre-Check
################################################################################
function errors_args(mme,methods,Pi)
    if methods == "conventional (no markers)" && mme.M!=0
        error("Conventional analysis runs without genotypes!")
    end
    if mme.M!=0 && methods=="GBLUP" && mme.M.G_is_marker_variance==true
        error("Please provide values for the genetic variance for GBLUP analysis")
    end
    if mme.nModels > 1 && mme.M!=0
        if Pi != 0.0 && round(sum(values(Pi)),digits=2)!=1.0
          error("Summation of probabilities of Pi is not equal to one.")
        end
    end
end

function check_pedigree(mme,df,pedigree)
    if mme.ped == 0 && pedigree == false
        return
    end
    if pedigree!=false
        pedID=map(string,collect(keys(pedigree.idMap)))
    else
        pedID=map(string,collect(keys(mme.ped.idMap)))
    end

    if mme.M!=0 && !issubset(mme.M.genoID,pedID)
        error("Not all genotyped individuals are found in pedigree!")
    end

    phenoID = map(String,df[1])
    if !issubset(phenoID,pedID)
        error("Not all phenotyped individuals are found in pedigree!")
    end
end

function check_outputID(outputEBV,mme)
    #Genotyped individuals are usaully not many, and are used in GWAS (complete
    #and incomplete), thus are used as default output_ID if not provided
    if outputEBV == false
        mme.output_ID = 0
    elseif mme.output_ID == 0 && mme.M != 0
        mme.output_ID = mme.M.genoID
    elseif mme.output_ID == 0 && mme.M == 0 && mme.pedTrmVec != 0
        #output EBV for all individuals in the pedigree for PBLUP
        pedID=map(String,collect(keys(mme.ped.idMap)))
        mme.output_ID = pedID
    end
end

function check_phenotypes(mme,df,single_step_analysis=false)
    #mme.M.obsID is
    #IDs for M after imputation, i.e., all individuals in pedigree, in SSBR
    #IDs for all genotyped individuals in complete genomic data
    phenoID = map(string,df[1])#same to df[:,1] in deprecated CSV
    if mme.M == 0
        return
    end
    if single_step_analysis == false
        if !issubset(phenoID,mme.M.genoID)
            printstyled("Phenotyped individuals are not a subset of\n",
            "genotyped individuals (complete genomic data,non-single-step).\n",
            "Only use phenotype information for genotyped individuals.\n",bold=false,color=:red)
            index = [phenoID[i] in mme.M.genoID for i=1:length(phenoID)]
            df    = df[index,:]
        elseif mme.output_ID!=0 && !issubset(mme.output_ID,mme.M.genoID)
            printstyled("Testing individuals are not a subset of \n",
            "genotyped individuals (complete genomic data,non-single-step).\n",
            "Only output EBV for tesing individuals with genotypes.\n",bold=false,color=:red)
            mme.output_ID = intersect(mme.output_ID,mme.M.genoID)
        end
    else
        pedID = map(string,collect(keys(ped.idMap)))
        if !issubset(phenoID,pedID)
            printstyled("Phenotyped individuals are not a subset of\n",
            "individuals in pedigree (incomplete genomic data (single-step) or PBLUP).\n",
            "Only use phenotype information for individuals in the pedigree.\n",bold=false,color=:red)
            index = [phenoID[i] in pedID for i=1:length(phenoID)]
            df    = df[index,:]
        elseif mme.output_ID!=0 && !issubset(mme.output_ID,pedID)
            printstyled("Testing individuals are not a subset of \n",
            "individuals in pedigree (incomplete genomic data (single-step) or PBLUP).\n",
            "Only output EBV for tesing individuals in the pedigree.\n",bold=false,color=:red)
            mme.output_ID = intersect(mme.output_ID,pedID)
        end
    end
end

function pre_check(mme,df,sol)
    getMME(mme,df)
    #starting value for sol can be provided
    if sol == false #no starting values
        sol = zeros(size(mme.mmeLhs,1))
    else            #besure type is Float64
        sol = map(Float64,sol)
    end
    return sol,df
end

################################################################################
# Return Output Results (Dictionary)
################################################################################
function output_result(mme,solMean,meanVare,G0Mean,output_samples_frequency,
                       meanAlpha,meanVara,estimatePi,mean_pi,output_file="MCMC_samples")
  output = Dict()
  if mme.output_ID != 0 &&  (mme.pedTrmVec != 0 || mme.M != 0)
      for traiti in 1:mme.nModels
          output["EBV"*"_"*string(mme.lhsVec[traiti])]=zeros(length(mme.output_ID))
      end
  end

  location_parameters = reformat2DataFrame([getNames(mme) solMean])
  output["Posterior mean of location parameters"] = location_parameters
  output["Posterior mean of residual variance"]   = meanVare
  if mme.pedTrmVec != 0
    output["Posterior mean of polygenic effects covariance matrix"]=G0Mean
    if mme.output_ID != 0
        for pedtrm in mme.pedTrmVec
            traiti, effect = split(pedtrm,':')
            sol_pedtrm     = map(Float64,location_parameters[(location_parameters[:Effect].==effect).&(location_parameters[:Trait].==traiti),:Estimate])
            EBV_pedtrm     = mme.output_X[pedtrm]*sol_pedtrm
            output["EBV"*"_"*string(mme.lhsVec[parse(Int64,traiti)])] += EBV_pedtrm
        end
    end
  end

  #samples for non-marker effects
  if output_samples_frequency != 0
      for i in  mme.outputSamplesVec
          trmi   = i.term
          trmStr = trmi.trmStr
          writedlm(output_file*"_"*trmStr*".txt",[transubstrarr(getNames(trmi))
                                                  i.sampleArray],',')
      end
  end

  if mme.M != 0
    if mme.nModels == 1
        meanAlpha=[meanAlpha] # make st array of array
    end
    markerout        = []
    if mme.M.markerID[1]!="NA"
        for markerArray in meanAlpha
          push!(markerout,[mme.M.markerID markerArray])
        end
    else
        for markerArray in meanAlpha
          push!(markerout,markerArray)
        end
    end

    output["Posterior mean of marker effects"] = (mme.nModels==1) ? markerout[1] : markerout
    output["Posterior mean of marker effects variance"] = meanVara
    if estimatePi == true
        output["Posterior mean of Pi"] = mean_pi
    end

    if mme.output_ID != 0
        for traiti in 1:mme.nModels
            EBV_markers  = mme.output_genotypes*meanAlpha[traiti] #fixed for mt
            output["EBV"*"_"*string(mme.lhsVec[traiti])] += EBV_markers
        end
    end
  end

  if mme.output_ID != 0 && haskey(mme.output_X,"J") #single-step analyis
      for traiti in 1:mme.nModels
          sol_J        = map(Float64,location_parameters[(location_parameters[:Effect].=="J").&(location_parameters[:Trait].==string(traiti)),:Estimate])[1]
          sol_ϵ        = map(Float64,location_parameters[(location_parameters[:Effect].=="ϵ").&(location_parameters[:Trait].==string(traiti)),:Estimate])
          EBV_J        = mme.output_X["J"]*sol_J
          EBV_ϵ        = mme.output_X["ϵ"]*sol_ϵ
          output["EBV"*"_"*string(mme.lhsVec[traiti])] += (EBV_J+EBV_ϵ)
      end
  end

  if mme.output_ID != 0 &&  (mme.pedTrmVec != 0 || mme.M != 0)
      for traiti in 1:mme.nModels
          EBV = output["EBV"*"_"*string(mme.lhsVec[traiti])]
          if EBV != zeros(length(mme.output_ID))
              output["EBV"*"_"*string(mme.lhsVec[traiti])]= [mme.output_ID EBV]
          end
      end
  end
  return output
end

################################################################################
# Reformat Output Array to DataFrame
################################################################################
function reformat2DataFrame(res::Array)
    ##SLOW, 130s
    #out_names=[strip(i) for i in split(res[1,1],':',keepempty=false)]
    #for rowi in 2:size(res,1)
    #    out_names=[out_names [strip(i) for i in split(res[rowi,1],':',keepempty=false)]]#hcat two vectors
    #end
    out_names = Array{String}(undef,size(res,1),3)
    for rowi in 1:size(res,1)
        out_names[rowi,:]=[strip(i) for i in split(res[rowi,1],':',keepempty=false)]
    end

    if size(out_names,2)==1 #convert vector to matrix
        out_names = reshape(out_names,length(out_names),1)
    end
    #out_names=permutedims(out_names,[2,1]) #rotate
    out_values=map(Float64,res[:,2])
    out=[out_names out_values]
    out = DataFrame(out, [:Trait, :Effect, :Level, :Estimate])
    return out
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
      #pi_selected = filter((k,v)->k[i]==1.0 && k[j]==1.0,Pi)
      pi_selected = filter(d->d.first[i]==1.0 && d.first[j]==1.0,Pi)

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
