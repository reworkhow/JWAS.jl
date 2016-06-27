type ModelTerm
    trmStr::AbstractString         #"A" ; "A*B"
    nFactors::Int64                # 1  ;  2
    factors::Array{Symbol,1}       #:A  ; :A,:B
    
    str::Array{AbstractString,1}   #covariate: str-> ["A x B", "A X B", ...]; val -> df[:A].*df[:B]
    val::Array{Float64,1}          #factor: str->["A1 x B1", "A2 X B2", ...]; val -> [1.0,1.0,...]
                                   #factor&covariate: str->["A x B1","A X B2", ...]; val->1.0.*df[:B]
    
    startPos::Int64                #start postion for this term in incidence matrix                       
    nLevels::Int64
    X::SparseMatrixCSC{Float64,Int64}
    names::Array{Any,1}            #names for this variable: 
                                   #covariate:     nLevels=1,       "A x B"; 
                                   #factor:        nLevels=nLevels, "A1 x B1", "A2 X B2", ...;
                                   #animal (ped) : nLevels=nAnimals       
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

type MCMCSamples
    term::ModelTerm
    sampleArray::Array{Float64,2}
end

type MME
    modelEquation::AbstractString      #"y= A + B + A*B"
    modelTerms::Array{ModelTerm,1}     # "A" ; "B"; "A*B";
    modelTermDict::Dict{AbstractString,ModelTerm}
    lhs::Symbol                        #:y
    covVec::Array{Symbol,1}            #variables those are covariates
    pedTrmVec::Array{AbstractString,1} #random variables related to A (pedigree);"Animal","Animal*Age"

    rndTrmVec::Array{RandomEffect,1}
    outputSamplesVec::Array{MCMCSamples,1}
    genVarSampleArray::Array{Float64,2}
    resVarSampleArray::Array{Float64,1}    
    
    
    X                                  #
    ySparse                            #
    mmeLhs                             #
    mmeRhs                             #
    
    ped                                #PedModule.Pedigree
    
    GiOld::Array{Float64,2}            #inverse of genetic covariance matrix for pedTrmVec
    GiNew::Array{Float64,2}            
    ROld::Float64                      #residual variance
    RNew::Float64                      
    
    
    Ai                                 #inverse of numerator relationship matrix
    mmePos::Int64                      #temporary value to record term positions
    M                                  #marker genotypes
    
    function MME(modelEquation,modelTerms,dict,lhs,RNew)
        new(modelEquation,modelTerms,dict,lhs,[],[],[],[],
            Array(Float64,1,1),Array(Float64,1),0,0,0,0,0,
        Array(Float64,1,1),Array(Float64,1,1),0.0,RNew,0,1,0)
    end
end