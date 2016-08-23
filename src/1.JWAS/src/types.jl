#Type for one term in model equations, e.g. age or comtempary group
type ModelTerm
    iModel::Int64                  # 1st, 2nd model...

    trmStr::AbstractString         #"1:A"  ;   1; :A
    nFactors::Int64                #"1:A*B";   2; :A,:B
    factors::Array{Symbol,1}       #:2:A   ;   1; :A

    str::Array{AbstractString,1}   #covariate       : str->["A x B", "A X B", ...];     val -> df[:A].*df[:B]
    val::Array{Float64,1}          #factor          : str->["A1 x B1", "A2 X B2", ...]; val -> [1.0,1.0,...]
                                   #factor&covariate: str->["A x B1","A X B2", ...];    val -> 1.0.*df[:B]

    nLevels::Int64                 #covariate   : nLevels ->1;       names -> "A x B";
    names::Array{Any,1}            #factor      : nLevels ->nLevels; names -> "A1 x B1", "A2 X B2", ...;
                                   #animal (ped): nLevels ->nAnimals;names -> ids

    startPos::Int64                   #start postion for this term in incidence matrix
    X::SparseMatrixCSC{Float64,Int64} #incidence matrix

    function ModelTerm(trmStr,m)
        iModel    = m
        trmStr    = strip(trmStr)
        factorVec = split(trmStr,"*")
        nFactors  = length(factorVec)
        factors   = [symbol(strip(f)) for f in factorVec]
        trmStr    = string(m)*":"*trmStr
        new(iModel,trmStr,nFactors,factors,[],[],0,[],0,spzeros(0,0))
    end
end

#using same incidence matrix for all traits
#Modify Ri based on missing pattern (number of Ri= 2^nTrait-1)
type ResVar
    R0::Array{Float64,2}
    RiDict::Dict{BitArray{1},Array{Float64,2}}
end

#general (iid) random effects
#single-trait
#multi-trait
type RandomEffect
    term::ModelTerm
    vcOld::Float64
    vcNew::Float64
    df::Float64
    scale::Float64
    sampleArray::Array{Float64,2} #for variance components
end

#MCMC samplers for location parameters
type MCMCSamples
    term::ModelTerm
    sampleArray::Array{Float64,2}
end

type Genotypes
  obsID::Array{UTF8String,1}    #row ID of genotypes
  #obsID  #now maybe string or int
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

#type for MME
#single-trait: lambda version of mixed model equations
#multi-trait : formal version
type MME
    nModels::Int64
    modelVec::Array{AbstractString,1}             #vector of model equations
    modelTerms::Array{ModelTerm,1}                #"1:intercept","1:A","2:intercept","2:A","2:A*B"...;
    modelTermDict::Dict{AbstractString,ModelTerm} #key: "1:A*B" value: ModelTerm
    lhsVec::Array{Symbol,1}                       #[:y1; :y2; ...]
    covVec::Array{Symbol,1}                       #variables those are covariates

    X                                             #Mixed Model Equations
    ySparse
    mmeLhs
    mmeRhs

    pedTrmVec::Array{AbstractString,1}      #random variables(pedigree);"1:Animal","1:Mat","2:Animal"
    ped                                     #PedModule.Pedigree
    Ai                                      #inv of numerator relationship matrix
    Gi::Array{Float64,2}                    #inv of genetic covariance matrix for pedTrmVec (multi-trait)
    GiOld::Array{Float64,2}                 #single-trait
    GiNew::Array{Float64,2}                 #(specific for lambda version of MME)

    rndTrmVec::Array{RandomEffect,1}        #iid random effects

    R::Array{Float64,2}                     #residual covariance matrix (multi-trait)
    missingPattern
    resVar
    ROld::Float64                           #residual variance (single-trait)
    RNew::Float64                           #for lambda MME

    M                                       #Genotypes

    mmePos::Int64                           #temporary value to record term position
                                            #(starting from 1)

    samples4R::Array{Float64,2}             #residual variance  (ndim^2 * niter)
    samples4G::Array{Float64,2}             #polygenic variance (matrix -> vector)
    outputSamplesVec::Array{MCMCSamples,1}  #location parameters

    df::DF

    function MME(nModels,modelVec,modelTerms,dict,lhsVec,R,ν)
      if nModels==1 && typeof(R)==Float64 #single-trait
        return new(nModels,modelVec,modelTerms,dict,lhsVec,[],0,0,0,0,[],0,0,zeros(1,1),zeros(1,1),zeros(1,1),[],zeros(1,1),0,0,0.0,R,0,1,zeros(1,1),zeros(1,1),[],DF(ν,4,4,4))
    elseif nModels>1 && typeof(R)==Array{Float64,2} #multi-trait
        return new(nModels,modelVec,modelTerms,dict,lhsVec,[],0,0,0,0,[],0,0,zeros(1,1),zeros(1,1),zeros(1,1),[],R,0,0,0.0,0.0,0,1,zeros(1,1),zeros(1,1),[],DF(ν,4,4,4))
    else
        error("Residual variance R should be a scalar for single-trait analyses or a matrix for multi-trait analyses.")
      end
    end
end
