"""
    build_model(model_equations::AbstractString,R)

* build models from **model equations** with residual varainces **R**

```julia
#single-trait
model_equations = "BW = intercept + age + sex"
R               = 6.72
models          = build_model(model_equations,R);

#multi-trait
model_equations = "BW = intercept + age + sex;
                   CW = intercept + age + sex";
R               = [6.72   24.84
                   24.84  708.41]
models          = build_model(model_equations,R);
```
"""
function build_model(model_equations::AbstractString,R)
  initMME(model_equations,R)
end

"""
    set_covariate(mme::MME,covStr::AbstractString)

* set variables as covariates; covStr is a string of varaibles separated by spaces


```julia
model_equations = "BW = intercept + age + sex;
                   CW = intercept + age + sex";
R               = [6.72   24.84
                   24.84  708.41]
models          = build_model(model_equations,R)
#set the variable age as covariates
set_covariate(models,"age")
```
"""
function set_covariate(mme::MME,covStr::AbstractString)
    covList(mme,covStr)
end

"""
    set_random(mme::MME,randomStr::AbstractString,ped::Pedigree, G)

* set variables as random polygenic effects

```julia
model_equations = "BW = intercept + age + sex + Animal;
                   CW = intercept + age + sex + Animal";
model           = build_model(model_equations,R);
ped             = get_pedigree(pedfile);
G               = [6.72   2.84
                   2.84  8.41]

setAsRandom(model,"Animal", ped,G)
```
"""
function set_random(mme::MME,randomStr::AbstractString,ped::PedModule.Pedigree, G)
    setAsRandom(mme,randomStr,ped, G)
end

"""
    set_random(mme::MME,randomStr::AbstractString,vc::Float64, df::Float64))

set variables as iid random effects
"""
function set_random(mme::MME,randomStr::AbstractString, vc::Float64, df::Float64)
    setAsRandom(mme,randomStr, vc, df)
end

"""
    get_pedigree(pedfile::AbstractString)
* Get pedigree informtion from a pedigree file.
* File format:

```
a 0 0
b 0 0
c a b
d a c
```
"""
function get_pedigree(pedfile::AbstractString)
  PedModule.mkPed(pedfile)
end

"""
    add_markers(mme::MME,file,G::Array{Float64,2});separator=' ',header=true,G_is_marker_variance=false)
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
function add_markers(mme::MME,file,G::Array{Float64,2};separator=' ',header=true,G_is_marker_variance=false)
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
    solve(mme::MME,df::DataFrame;solver="default",printout_frequency=100,tolerance = 0.000001,maxiter = 5000)

* Solve the mixed model equations (no marker information) without estimating variance components.
Available solvers includes `default`,`Jacobi`,`GaussSeidel`,`Gibbs sampler`.
"""
function solve(mme::MME,
                df::DataFrame;
                solver="default",
                printout_frequency=100,
                tolerance = 0.000001,
                maxiter = 5000)
    if size(mme.mmeRhs)==()
        getMME(mme,df)
    end
    p = size(mme.mmeRhs,1)
    if solver=="Jacobi"
        return [getNames(mme) Jacobi(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,0.3,
                                    tolerance=tolerance,outFreq=printout_frequency,maxiter=maxiter)]
    elseif solver=="GaussSeidel"
        return [getNames(mme) GaussSeidel(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,
                              tolerance=tolerance,outFreq=printout_frequency,maxiter=maxiter)]
    elseif solver=="Gibbs" && mme.nModels !=1
        return [getNames(mme) Gibbs(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,maxiter,
                              outFreq=printout_frequency)]
    elseif solver=="Gibbs" && mme.nModels==1
        return [getNames(mme) Gibbs(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,mme.RNew,maxiter,outFreq=printout_frequency)]
    elseif solver=="default"
        return [getNames(mme) mme.mmeLhs\mme.mmeRhs]
    else
        error("No this solver\n")
    end
end

"""
    runMCMC(mme,df;Pi=0.0,estimatePi=false,chain_length=1000,starting_value=false,printout_frequency=100,missing_phenotypes=false,constraint=false,methods="conventional (no markers)",output_samples_frequency::Int64 = 0)

Run MCMC (marker information included or not) with sampling of variance components.

* available **methods** include "conventional (no markers)", "BayesC0", "BayesC", "BayesCC","BayesB".
* **missing_phenotypes**
* **Pi** is a dictionary such as `Pi=Dict([1.0; 1.0]=>0.7,[1.0; 0.0]=>0.1,[0.0; 1.0]=>0.1,[0.0; 0.0]=>0.1)`
* save MCMC samples of variance components and marker effects every **output_samples_frequency** iterations to files
* **starting_value** can be provided as a vector for all location parameteres except marker effects.
* print out the monte carlo mean in REPL with **printout_frequency**
* **constraint**=true if constrain residual covariances between traits to be zero.
"""

function runMCMC(mme,df;
                Pi                = 0.0,   #Dict{Array{Float64,1},Float64}()
                chain_length      = 100,
                starting_value    = false,
                printout_frequency= chain_length+1,
                missing_phenotypes= false,
                constraint        = false,
                estimatePi        = false,
                methods           = "conventional (no markers)",
                output_samples_frequency::Int64 = 0)

  if mme.M != 0 && mme.M.G_is_marker_variance==false
    genetic2marker(mme.M,Pi)
    println("Priors for marker effects covariance matrix were calculated from genetic covariance matrix and π.")
    if !isposdef(mme.M.G)
      error("Marker effects covariance matrix is not postive definite! Please modify the argument: Pi.")
    end
    println("Marker effects covariance matrix is ")
    println(round(mme.M.G,6),".\n\n")
  end

  if mme.nModels ==1
      if methods in ["BayesC","BayesC0","conventional (no markers)"]
          res=MCMC_BayesC(chain_length,mme,df,
                          π          =Pi,
                          methods    =methods,
                          estimatePi =estimatePi,
                          sol        =starting_value,
                          outFreq    =printout_frequency,
                          output_marker_effects_frequency =output_marker_effects_frequency)
       elseif methods =="BayesB"
           res=MCMC_BayesB(chain_length,mme,df,Pi,
                            estimatePi =false,
                            sol        =starting_value,
                            outFreq    =printout_frequency,
                            output_marker_effects_frequency =output_marker_effects_frequency)
        else
            error("No options!!!")
        end
    elseif mme.nModels > 1
        if Pi != 0.0 && sum(values(Pi))!=1.0
          error("Summation of probabilities of Pi is not equal to one.")
        end
        if methods in ["BayesC","BayesCC","BayesC0","conventional (no markers)"]
          res=MT_MCMC_BayesC(chain_length,mme,df,
                          Pi     = Pi,
                          sol    = starting_value,
                          outFreq= printout_frequency,
                          missing_phenotypes=missing_phenotypes,
                          constraint = constraint,
                          estimatePi = estimatePi,
                          methods    = methods,
                          output_samples_frequency=output_samples_frequency)
        elseif methods=="BayesB"
            if Pi == 0.0
                error("Pi is not provided!!")
            end
            res=MT_MCMC_BayesB(chain_length,mme,df,Pi,
                            sol=starting_value,
                            outFreq=printout_frequency,
                            missing_phenotypes=missing_phenotypes,
                            constraint=constraint,
                            output_marker_effects_frequency=output_samples_frequency)
        else
            error("No methods options!!!")
        end
    else
        error("No options!")
    end
  res
end

"""
    showMME(mme::MME,df::DataFrame)

* Show left-hand side and right-hand side of mixed model equations (no markers).
"""
function showMME(mme::MME,df::DataFrame)
   if size(mme.mmeRhs)==()
     getMME(mme,df)
   end
    return mme.mmeLhs,mme.mmeRhs
end
