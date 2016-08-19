"""
    set_random(mme::MME,randomStr::AbstractString,vc::Float64, df::Float64))

set variables as iid random effects
"""
function set_random(mme::MME,randomStr::AbstractString, vc::Float64, df::Float64)
    setAsRandom(mme,randomStr, vc, df)
end

"""
    add_markers(mme::MME,file,G::Float64;separator=' ',header=true)

* Get marker informtion from a genotype file (same order as the phenotype file).
* G is the additive genetic variance.
* File format:

```
Animal,marker1,marker2,marker3,marker4,marker5
S1,1,0,1,1,1
D1,2,0,2,2,1
O1,1,2,0,1,0
O3,0,0,2,1,1
```
"""
function add_markers(mme::MME,file,G::Float64;separator=' ',header=true,G_is_marker_variance=false)
    addMarkers(mme,file,G,separator=separator,header=header,G_is_marker_variance=G_is_marker_variance)
end

"""
    outputMCMCsamples(mme::MME,trmStr::AbstractString...)

Get samples for specific variables.
"""
function outputMCMCsamples(mme::MME,trmStr::AbstractString...)
    for i in trmStr
      outputSamplesFor(mme,i)
    end
end

"""
    runMCMC(mme,df;Pi=0.0,chain_length=1000,starting_value=false,printout_frequency=100,estimatePi=false,methods="no markers",output_marker_effects_frequency::Int64 = 0)

Run MCMC (marker information included or not) with sampling of variance components.
Available methods include "no markers", "BayesB", "BayesC".
"""
function runMCMC(mme,df;
                Pi                =0.0,
                chain_length      =1000,
                starting_value    =false,
                printout_frequency=chain_length+1,
                estimatePi        =false,
                methods           ="conventional (no markers)", #BayesC0,BayesC,BayesCπ,BayesB
                output_marker_effects_frequency::Int64 = 0 # 0=>save samples to a file
                )
  if mme.M ==0
    #res=MCMC_conventional(chain_length,mme,df,
    res=MCMC_BayesC(chain_length,mme,df,
                          sol=starting_value,
                          outFreq=printout_frequency)
  elseif methods=="BayesC" || methods=="BayesC0"
    res=MCMC_BayesC(chain_length,mme,df,
                    π      =Pi,
                    methods=methods,
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
