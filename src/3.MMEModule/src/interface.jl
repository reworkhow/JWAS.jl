#build model form model equations
function build_model(models::AbstractString,R::Float64)
  initMME(models,R)
end

#set as covariate
function set_covariate(mme::MME,covStr::AbstractString)##string seprated by space
    covList(mme,covStr)
end

#set random effects
function set_random(mme::MME,randomStr::AbstractString,ped::PedModule.Pedigree, G)
    setAsRandom(mme,randomStr,ped, G)
end

function set_random(mme::MME,randomStr::AbstractString, vc::Float64, df::Float64)
    setAsRandom(mme,randomStr, vc, df)
end

#add pedigree information
function get_pedigree(pedfile::AbstractString)
  PedModule.mkPed(pedfile)
end

#add marker information
function add_markers(mme::MME,file,G::Float64;separator=' ',header=true)
    addMarkers(mme,file,G,separator=separator,header=header)
end

#get MME without MCMC_markers
function showMME(mme::MME,df::DataFrame)
   if size(mme.mmeRhs)==()
     getMME(mme,df)
   end
    return mme.mmeLhs,mme.mmeRhs
end

#get samples for specific variables
function outputMCMCsamples(mme::MME,trmStr::AbstractString...)
    for i in trmStr
      outputSamplesFor(mme,i)
    end
end


#solve equations without estimate variance components
function solve(mme::MME,
                df::DataFrame;
                solver="Jacobi",
                printout_frequency=100,
                tolerance = 0.000001,
                niterations = 5000
                )
    if size(mme.mmeRhs)==()
        getMME(mme,df)
    end
    p = size(mme.mmeRhs,1)
    if solver=="Jacobi"
        return [getNames(mme) Jacobi(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,0.3,tolerance=tolerance,output=printout_frequency)]
    elseif solver=="GaussSeidel"
        return [getNames(mme) GaussSeidel(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,tolerance=tolerance,output=printout_frequency)]
    elseif solver=="Gibbs"
        return [getNames(mme) Gibbs(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,mme.RNew,niterations,outFreq=printout_frequency)]
    else
        error("No this solver\n")
    end
end


function runMCMC(mme,df;
                Pi                =0.0,
                chain_length      =1000,
                starting_value    =false,
                printout_frequency=100,
                estimatePi        =false,
                methods           ="no markers", #BayesC0,BayesC,BayesCπ,BayesB
                output_marker_effects_frequency::Int64 = 0 # 0=>save samples to a file
                )
  if mme.M ==0
    res=MCMC_conventional(chain_length,mme,df,
                          sol=starting_value,
                          outFreq=printout_frequency,
                          thin=thin)
  elseif methods=="BayesC" && Pi == 0.0
    res=MCMC_BayesC0(chain_length,mme,df,
                     sol=starting_value,
                     outFreq=printout_frequency,
                     output_marker_effects_frequency =output_marker_effects_frequency)
  elseif methods=="BayesC" && Pi != 0.0
    res=MCMC_BayesC(chain_length,mme,df,Pi,
                     estimatePi = estimatePi,
                     sol=starting_value,
                     outFreq=printout_frequency,
                     output_marker_effects_frequency =output_marker_effects_frequency)
  elseif methods=="BayesB"
    res=MCMC_BayesB(chain_length,mme,df,Pi,
                     estimatePi = false,
                     sol=starting_value,
                     outFreq=printout_frequency,
                     output_marker_effects_frequency =output_marker_effects_frequency)
  else
    error("No options!!!")
  end
  res
end
