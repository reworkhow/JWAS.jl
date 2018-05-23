
################################################################################
# model__equations = "y1 = A + B
#                     y2 = A + B + A*B"
#
# the class ModelTrem is shown below for one term in model_equations, such as
# terms: 1:A,1:B,2:A,2:B,2:A*B
#
################################################################################
type ModelTerm
    iModel::Int64                  # 1st or 2nd model_equation

                                   # | trmStr | nFactors | factors |
                                   # |--------|----------|---------|
    trmStr::AbstractString         # | "1:A"  |    1     | :A      |
    nFactors::Int64                # | "2:A"  |    1     | :A      |
    factors::Array{Symbol,1}       # | "1:A*B"|    2     | :A,:B   |

                                   #DATA             |          str               |     val       |
                                   #                :|----------------------------|---------------|
    str::Array{AbstractString,1}   #covariate^2     :|["A x B", "A X B", ...]     | df[:A].*df[:B]|
    val::Array{Float64,1}          #factor^2        :|["A1 x B1", "A2 X B2", ...] | [1.0,1.0,...] |
                                   #factor*covariate:|["A1 x B","A2 X B", ...]    | 1.0.*df[:B]   |

                                   #OUTPUT           | nLevels |     names        |
                                   #                 |---------|------------------|
    nLevels::Int64                 #covariate   :    | 1       | "A"              |
    names::Array{Any,1}#for OUTPUT #factor      :    | nLevels | "A1", "A2", ...  |
                                   #animal (ped):    | nAnimals| ids              |
                                   #animal(ped)*age: | nAnimals| "A1*age","A2*age"|
                                   #factor*covariate:| nLevels | "A1*age","A2*age"|

    startPos::Int64                   #start postion for this term in incidence matrix
    X::SparseMatrixCSC{Float64,Int64} #incidence matrix

    function ModelTerm(trmStr,m)
        iModel    = m
        trmStr    = strip(trmStr)
        factorVec = split(trmStr,"*")
        nFactors  = length(factorVec)
        factors   = [Symbol(strip(f)) for f in factorVec]
        trmStr    = string(m)*":"*trmStr
        new(iModel,trmStr,nFactors,factors,[],[],0,[],0,spzeros(0,0))
    end
    #e.g. animal*age, nLevels
end

################################################################################
#A class for residual covariance matrix for all observations of size (nob*nModel)
#where Ri is modified based on missing pattern (number of Ri= 2^nTrait-1)
#It allows using the same incidence matrix X for all traits in multi-trait analyses
#In JWAS, ONLY used when residual variance is constant
#or missing phenotypes are not imputed at each step of MCMC (no marker effects).
################################################################################
type ResVar
    R0::Array{Float64,2}
    RiDict::Dict{BitArray{1},Array{Float64,2}}
end

################################################################################
#General (iid) random effects
#NOW assuming different variables are independent cov(A,B)=0 #NO
#single-trait e.g. termarray: [ModelTerm(1:A)]
#multi-trait  e.g. termarray: [ModelTerm(1:A), ModelTerm(2:A)]
################################################################################
type RandomEffect
    term_array::Array{ModelTerm,1} #term_array::::Array{AbstractString,1}
    G::Array{Float64,2}      #covariance matrix (multi-trait)
    GiOld::Array{Float64,2}  #specific for lambda version of MME (single-trait)
    GiNew::Array{Float64,2}  #specific for lambda version of MME (single-trait)
    #vcOld #::Array{Float64,2}
    #vcNew #::Array{Float64,2}
    df::Float64
    scale #::Array{Float64,2}
end

# type RandomEffect
#     pedTrmVec::Array{AbstractString,1}            #polygenic effects(pedigree): "1:Animal","1:Mat","2:Animal"
#     Vi                                            #inverse of covariance matrix
#     Gi::Array{Float64,2}                          #inverse of genetic covariance matrix for pedTrmVec (multi-trait)
#     GiOld::Array{Float64,2}                       #specific for lambda version of MME (single-trait)
#     GiNew::Array{Float64,2}                       #specific for lambda version of MME (single-trait)
#
#     df::Float64
#     scale #::Array{Float64,2}
# end

#MCMC samplers for location parameters
type MCMCSamples
    term::ModelTerm
    sampleArray::Array{Float64,2}
end

type Genotypes
  obsID::Array{String,1}    #row ID of genotypes
  markerID
  nObs::Int64
  nMarkers::Int64
  alleleFreq::Array{Float64,2}
  sum2pq::Float64
  centered::Bool
  genotypes::Array{Float64,2}
  G  #ST->Float64;MT->Array{Float64,2}
  G_is_marker_variance::Bool
  Genotypes(a1,a2,a3,a4,a5,a6,a7,a8)=new(a1,a2,a3,a4,a5,a6,a7,a8,0.0,false)
end

type DF
    residual::Float64
    polygenic::Float64
    marker::Float64
    random::Float64
end

################################################################################
# the class MME is shown below with members for models, mixed model equations...
#
# single-trait analyses: lambda version of MME
#  multi-trait analysis: formal version of MME
################################################################################
type MME
    nModels::Int64                                #number of model equations
    modelVec::Array{AbstractString,1}             #["y1 = A + B","y2 = A + B + A*B"]
    modelTerms::Array{ModelTerm,1}                #ModelTerms for "1:intercept","1:A","2:intercept","2:A","2:A*B"...;
    modelTermDict::Dict{AbstractString,ModelTerm} #key: "1:A*B" value: ModelTerm; convert modelTerms above to dictionary
    lhsVec::Array{Symbol,1}                       #phenotypes: [:y1,:y2]
    covVec::Array{Symbol,1}                       #variables those are covariates

                                                  #MIXED MODEL EQUATIONS
    X                                             #incidence matrix
    ySparse                                       #phenotypes
    IDs                                           #IDs for phenotypes
    mmeLhs                                        #Lhs of Mixed Model Equations
    mmeRhs                                        #Rhs of Mixed Model Equations

                                                  #RANDOM EFFCTS
    pedTrmVec::Array{AbstractString,1}            #polygenic effects(pedigree): "1:Animal","1:Mat","2:Animal"
    ped                                           #PedModule.Pedigree
    Ai                                            #inverse of numerator relationship matrix
    Gi::Array{Float64,2}                          #inverse of genetic covariance matrix for pedTrmVec (multi-trait)
    GiOld::Array{Float64,2}                       #specific for lambda version of MME (single-trait)
    GiNew::Array{Float64,2}                       #specific for lambda version of MME (single-trait)

    rndTrmVec::Array{RandomEffect,1}              #iid random effects
    #should have one for pedigree randomEffect
                                                  #RESIDUAL EFFECTS
    R::Array{Float64,2}                           #residual covariance matrix (multi-trait)
    missingPattern                                #for impuation of missing residual
    resVar                                        #for impuation of missing residual
    ROld::Float64                                 #residual variance (single-trait)
    RNew::Float64                                 #for lambda version of MME (single-trait)

    M                                             #GENOTYPES

    mmePos::Int64                                 #temporary value to record term position (start from 1)

                                                  #MCMC SAMPLES
    samples4R::Array{Float64,2}                   #residual variance  (ndim^2 * niter)
    samples4G::Array{Float64,2}                   #polygenic variance (matrix -> vector)
    outputSamplesVec::Array{MCMCSamples,1}        #location parameters
    #samples4π                                    #samples4π is not used in MME type since π is BayesC-specific

    df::DF                                        #prior degree of freedom

    function MME(nModels,modelVec,modelTerms,dict,lhsVec,R,ν)
      if nModels==1 && typeof(R)==Float64             #single-trait
        return new(nModels,modelVec,modelTerms,dict,lhsVec,[],
                   0,0,0,0,0,
                   [],0,0,zeros(1,1),zeros(1,1),zeros(1,1),
                   [],
                   zeros(1,1),0,0,R,R,
                   0,
                   1,
                   zeros(1,1),zeros(1,1),[],
                   DF(ν,4,4,4))
      elseif nModels>1 && typeof(R)==Array{Float64,2} #multi-trait
        return new(nModels,modelVec,modelTerms,dict,lhsVec,[],
                   0,0,0,0,0,
                   [],0,0,zeros(1,1),zeros(1,1),zeros(1,1),
                   [],
                   R,0,0,0.0,0.0,
                   0,
                   1,
                   zeros(1,1),zeros(1,1),[],
                   DF(ν,4,4,4))
      else
        error("Residual variance R should be a scalar for single-trait analyses or a matrix for multi-trait analyses.")
      end
    end
end
