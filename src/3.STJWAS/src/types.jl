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
