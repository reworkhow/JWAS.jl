include("samplePi.jl")
include("MCMC_conventional.jl")
include("MCMC_BayesC0.jl")
include("MCMC_BayesC.jl")
include("MCMC_BayesCC.jl")


#=
function runMCMC(mme,df;nIter=1000,sol=false,outFreq=100,thin=100)
  if mme.M !=0
    res=MCMC_markers(nIter,mme,df,sol=sol,outFreq=outFreq,thin=thin)
  else
    res=MCMC_conventional(nIter,mme,df,sol=sol,outFreq=outFreq,thin=thin)
  end
  res
end
=#
