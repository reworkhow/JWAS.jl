type ModelTerm
    iModel::Int64                  #1 ; 2; 3... (1st, 2nd model...)
    trmStr::AbstractString         #"1:A" ; "1:A*B" **
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

type MME
    modelVec::Array{AbstractString,1}
    modelTerms::Array{ModelTerm,1}                #"1:intercept","1:A","2:intercept","2:A"...;
    modelTermDict::Dict{AbstractString,ModelTerm}
    lhsVec::Array{Symbol,1}
    covVec::Array{Symbol,1}
    pedTrmVec::Array{AbstractString,1}            #"1:Animal","1:Mat","2:Animal"
    X
    ySparse
    mmeLhs
    mmeRhs
    ped
    Gi::Array{Float64,2}
    R::Array{Float64,2}
    Ai
    mmePos::Int64
    
    missingPattern
    resVar
    
    M #genotype
end


#using same incidence matrix for all traits
#Modify Ri based on missing pattern (number of Ri= 2^nTrait-1)
type ResVar 
    R0::Array{Float64,2}
    RiDict::Dict{BitArray{1},Array{Float64,2}}
end

