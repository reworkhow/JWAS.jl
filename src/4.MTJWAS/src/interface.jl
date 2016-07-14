"""
    build_model(model_equations::AbstractString,R::Array{Float64,2})

build models from **model equations** with a residual covaraince matrix **R**

```julia
model_equations = "BW = intercept + age + sex;
                   CW = intercept + age + sex";
R               = [6.72   24.84
                   24.84  708.41]
models          = MT.build_model(model_equations,R);
```
"""
function build_model(model_equations::AbstractString,R::Array{Float64,2})
  initMME(model_equations,R)
end

"""
    set_covariate(mme::MME,covStr::AbstractString)

set variables as covariates; covStr is a string of varaibles separated by spaces

```julia
model_equations = "BW = intercept + age + sex;
                   CW = intercept + age + sex";
R               = [6.72   24.84
                   24.84  708.41]
models          = MT.build_model(model_equations,R)

#set the variable age as covariates
MT.set_covariate(models,"age")
```
"""
function set_covariate(mme::MME,covStr::AbstractString)
    covList(mme,covStr)
end

"""
    set_random(mme::MME,randomStr::AbstractString,ped::PedModule.Pedigree, G::Array{Float64,2})

set variables as random polygenic effects

```julia
model_equations = "BW = intercept + age + sex + Animal;
                   CW = intercept + age + sex + Animal";
model           = MT.build_model(model_equations,R);
ped             = MT.get_pedigree(pedfile);

setAsRandom(model,"Animal", ped,G)
```
"""
function set_random(mme::MME,randomStr::AbstractString,ped::PedModule.Pedigree, G::Array{Float64,2})
    setAsRandom(mme,randomStr,ped, G)
end

"""
    get_pedigree(pedfile::AbstractString)
get pedigree informtion from a pedigree file
"""
function get_pedigree(pedfile::AbstractString)
  PedModule.mkPed(pedfile)
end

"""
    add_markers(mme::MME,file,G::Array{Float64,2});separator=' ',header=true,G_is_genetic_variance=true)
Get marker informtion from a genotype file (same order as the phenotype file).\\
File format:
Animal,marker1,marker2,marker3,marker4,marker5\\
S1,1,0,1,1,1\\
D1,2,0,2,2,1\\
O1,1,2,0,1,0\\
O3,0,0,2,1,1\\
"""
function add_markers(mme::MME,file,G::Array{Float64,2};separator=' ',header=true,G_is_genetic_variance=true)
    addMarkers(mme,file,G,separator=separator,header=header,G_is_genetic_variance=G_is_genetic_variance)
end


"""
    solve(mme::MME,df::DataFrame;solver="Jacobi",printout_frequency=100,tolerance = 0.000001,niterations = 5000)

Solve the mixed model equations (no marker information) without estimating variance components.
Available solvers includes `Jacobi`,`GaussSeidel`,`Gibbs sampler`.
"""
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

"""
    runMCMC(mme,df;Pi=0.0,chain_length=1000,starting_value=false,printout_frequency=100,missing_phenotypes= false,methods="no markers",output_marker_effects_frequency::Int64 = 0)

Run MCMC (marker information included or not) with sampling of variance components.

* available **methods** include "no markers", "BayesC0", "BayesC", "BayesCC".
* **missing_phenotypes**
* **Pi** is a dictionary such as `Pi=Dict([1.0; 1.0]=>0.7,[1.0; 0.0]=>0.1,[0.0; 1.0]=>0.1,[0.0; 0.0]=>0.1)`
* save samples of marker effects every **output_marker_effects_frequency** iterations to files
* **starting_value** can be provided as vector for all location parameteres except marker effects.
* print out the monte carlo mean in REPL with **printout_frequency**,
"""

function runMCMC(mme,df;
                Pi                =0.0,   #Dict{Array{Float64,1},Float64}()
                chain_length      =1000,
                starting_value    =false,
                printout_frequency=100,
                missing_phenotypes= false,
                methods           = "no markers", #BayesC0,BayesC,BayesCC
                output_marker_effects_frequency::Int64 = 0)
  if mme.M ==0
    res=MCMC_conventional(chain_length,mme,df,
                          sol=starting_value,
                          outFreq=printout_frequency,
                          missing_phenotypes=missing_phenotypes)
  elseif methods=="BayesC0"
    res=MCMC_BayesC0(chain_length,mme,df,
                     sol=starting_value,
                     outFreq=printout_frequency,
                     missing_phenotypes=missing_phenotypes,
                     output_marker_effects_frequency=output_marker_effects_frequency)
  elseif methods=="BayesC"
    if Pi == 0.0
      error("Pi is not provided!!")
    end
    res=MCMC_BayesC(chain_length,mme,df,Pi,
                     sol=starting_value,
                     outFreq=printout_frequency,
                     missing_phenotypes=missing_phenotypes,
                     output_marker_effects_frequency=output_marker_effects_frequency)
  elseif methods=="BayesCC"
    if Pi == 0.0
      error("Pi is not provided!!")
    end
    res=MCMC_BayesCC(chain_length,mme,df,Pi,
                     sol=starting_value,
                     outFreq=printout_frequency,
                     missing_phenotypes=missing_phenotypes,
                     output_marker_effects_frequency=output_marker_effects_frequency)
  elseif methods=="BayesB"
    if Pi == 0.0
      error("Pi is not provided!!")
    end
    res=MCMC_BayesB(chain_length,mme,df,Pi,
                     sol=starting_value,
                     outFreq=printout_frequency,
                     missing_phenotypes=missing_phenotypes,
                     output_marker_effects_frequency=output_marker_effects_frequency)
  else
    error("No methods options!!!")
  end
  res
end
