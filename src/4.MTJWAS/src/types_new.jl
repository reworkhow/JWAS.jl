type ModelTerm #Type for all terms in model equations
    iModel::Int64                  # 1st, 2nd model...
    trmStr::AbstractString         #"1:A" ; "1:A*B"

    nFactors::Int64                # 1    ;  2
    factors::Array{Symbol,1}       #:A    ; :A , :B

    str::Array{AbstractString,1}   #covariate       : str->["A x B", "A X B", ...];     val -> df[:A].*df[:B]
    val::Array{Float64,1}          #factor          : str->["A1 x B1", "A2 X B2", ...]; val -> [1.0,1.0,...]
                                   #factor&covariate: str->["A x B1","A X B2", ...];    val -> 1.0.*df[:B]

    startPos::Int64                #start postion for this term in incidence matrix
    nLevels::Int64                 #covariate   : nLevels=1,       "A x B";
                                   #factor      : nLevels=nLevels, "A1 x B1", "A2 X B2", ...;
                                   #animal (ped): nLevels=nAnimals

    X::SparseMatrixCSC{Float64,Int64}
    names::Array{Any,1}            #names for this variable

    ModelTerm(x1,x2)=new(x1,x2,0,[],[],[],0,0,spzeros(0,0),[])
end

type MME #Difference in single-trait and Multi-trait is whether using lambda version of MME
    modelVec::Array{AbstractString,1}             #vector of model equations
    modelTerms::Array{ModelTerm,1}                #"1:intercept","1:A","2:intercept","2:A","2:A*B"...;
    modelTermDict::Dict{AbstractString,ModelTerm} #key: "1:A*B" value: ModelTerm
    lhsVec::Array{Symbol,1}                       #[:y1; :y2; ...]
    covVec::Array{Symbol,1}                       #variables those are covariates

    X                                             #Mixed Model Equations
    ySparse
    mmeLhs
    mmeRhs

    pedTrmVec::Array{AbstractString,1}           #random variables based on (pedigree);"1:Animal","1:Mat","2:Animal"
    ped                                          #PedModule.Pedigree
    Ai                                           #inverse of numerator relationship matrix
    Gi::Array{Float64,2}                         #inverse of genetic covariance matrix for pedTrmVec (multu-trait)
    GiOld::Array{Float64,2}                      #single-trait
    GiNew::Array{Float64,2}                      #(specific for lambda version of MME)


    R::Array{Float64,2}                          #residual(multi-trait)
    missingPattern
    resVar
    ROld::Float64                                #single-trait
    RNew::Float64                                #(specific for lambda version of MME)


    M                                            #Genotype

    mmePos::Int64                                #temporary value to record term position (starting from 1)


    function MME(modelVec,modelTerms,dict,lhsVec,R)
      if length(modelVec)!=1 #multi-trait
        return new(modelVec,modelTerms,dict,lhsVec,[],0,0,0,0,[],0,0,zeros(1,1),zeros(1,1),zeros(1,1),
                   R,0,0,0.0,0.0,0,1)
      else #single-trait
        return new(modelVec,modelTerms,dict,lhsVec,[],0,0,0,0,[],0,0,zeros(1,1),zeros(1,1),zeros(1,1),
                   zeros(1,1),0,0,0.0,R,0,1)
      end
    end
end


#using same incidence matrix for all traits
#Modify Ri based on missing pattern (number of Ri= 2^nTrait-1)
type ResVar
    R0::Array{Float64,2}
    RiDict::Dict{BitArray{1},Array{Float64,2}}
end

type MCMCsamples
  samples4R::Array{Array{Float64,2},1} #residual variance
  samples4G::Array{Array{Float64,2},1} #polygenic variance
end

#general (iid) random effects; should also make a specific type for BV ped effects
type RandomEffect
    term::ModelTerm
    vcOld::Float64
    vcNew::Float64
    df::Float64
    scale::Float64
    sampleArray::Array{Float64,1}
end
