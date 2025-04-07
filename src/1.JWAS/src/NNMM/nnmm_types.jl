mutable struct Layer
    layer_name::String  
    data_path::String
    separator
    header
    data
    quality_control
    MAF
    missing_value
    center   
    function Layer(;layer_name, 
                    data_path,
                    separator=",",
                    header=true, 
                    data=false,
                    quality_control=true,
                    MAF=0.01,
                    missing_value=9.0,
                    center=true)
        new(layer_name, 
            data_path, 
            separator, 
            header, 
            data, 
            quality_control, 
            MAF, 
            missing_value,
            center)
    end
end


# mutable struct nnmmGenotypes
#     name                            #name for this category, eg. "geno1"
#     trait_names                     #names for the corresponding traits, eg.["y1","y2"]
  
#     obsID::Array{AbstractString,1}  #row ID for (imputed) genotyped and phenotyped inds (finally)
#     markerID
#     nObs::Int64                     #length of obsID
#     nMarkers::Int64
#     alleleFreq
#     sum2pq::AbstractFloat
#     centered::Bool
#     data::Union{Array{Float64,2},Array{Float32,2}}
#     nLoci             #number of markers included in the model
#     ntraits           #number of traits included in the model
#     nFeatures         #number of features included in the model, for genotype, same as nMarkers
  
#     genetic_variance  #genetic variance, type: Variance struct
#     G       #marker effect variance; the "Variance" object, for Variance.val: ST->Float64;MT->Array{Float64,2}

#     method            #prior for marker effects (Bayesian ALphabet, GBLUP ...)
#     estimatePi

#     mArray            #a collection of matrices used in Bayesian Alphabet
#     mRinvArray        #a collection of matrices used in Bayesian Alphabet
#     mpRinvm           #a collection of matrices used in Bayesian Alphabet
#     mΦΦArray          # a collection of matrices used in RRM
#     D                 #eigen values used in GBLUP
#     gammaArray        #array used in Bayesian LASSO
  
#     α                 #array of current MCMC samples
#     β
#     δ
#     π
  
#     meanAlpha         #arrays of results
#     meanAlpha2
#     meanDelta
#     mean_pi
#     mean_pi2
#     meanVara
#     meanVara2
#     meanScaleVara
#     meanScaleVara2
  
#     output_genotypes #output genotypes
  
#     isGRM  #whether genotypes or relationship matirx is provided
  
#     nnmmGenotypes(a1,a2,a3,a4,a5,a6,a7,a8,a9)=new(false,false,
#                                            a1,a2,a3,a4,a5,a6,a7,a8,a4,false,a4,
#                                            Variance(false,false,false,true,false,false),Variance(false,false,false,true,false,false), #false,false,
#                                            false,true, #true,false,
#                                            false,false,false,false,false,false,
#                                            false,false,false,false,
#                                            false,false,false,false,false,false,false,false,false,
#                                            false,a9)
# end


mutable struct Omics
    name                            #name for this category, eg. "geno1"
    trait_names                     #names for the corresponding traits, eg.["y1","y2"]
  
    obsID::Array{AbstractString,1}  #row ID for (imputed) genotyped and phenotyped inds (finally)
    featureID # = omicsID
    nObs::Int64                     #length of obsID
    nFeatures::Int64
    centered::Bool
    data::Union{Array{Float64,2},Array{Float32,2}}
    ntraits           #number of traits included in the model
  
    genetic_variance  #genetic variance, type: Variance struct
    G       #marker effect variance; the "Variance" object, for Variance.val: ST->Float64;MT->Array{Float64,2}
  
    method            #prior for marker effects (Bayesian ALphabet, GBLUP ...)
    estimatePi
  
    mArray            #a collection of matrices used in Bayesian Alphabet
    mRinvArray        #a collection of matrices used in Bayesian Alphabet
    mpRinvm           #a collection of matrices used in Bayesian Alphabet
    mΦΦArray          # a collection of matrices used in RRM
    D                 #eigen values used in GBLUP
    gammaArray        #array used in Bayesian LASSO
  
    α                 #array of current MCMC samples
    β
    δ
    π
  
    meanAlpha         #arrays of results
    meanAlpha2
    meanDelta
    mean_pi
    mean_pi2
    meanVara
    meanVara2
    meanScaleVara
    meanScaleVara2
  
    output_genotypes #output genotypes
  
    isGRM
    Omics(obsID,featureID,nObs,nFeatures,centered,data) = 
        new(false,false,
            obsID,featureID,nObs,nFeatures,centered,data,false,
            Variance(false,false,false,true,false,false),Variance(false,false,false,true,false,false), 
            false,true,
            false,false,false,false,false,false,
            false,false,false,false,
            false,false,false,false,false,false,false,false,false,
            false,false)
end

mutable struct Phenotypes
    obsID::Array{AbstractString,1}  #row ID for (imputed) genotyped and phenotyped inds (finally)
    featureID::Array{AbstractString,1} # = obsID
    nObs::Int64                     #length of obsID
    nPheno::Int64
    data
    nFeatures # = nPheno for Phenotypes = ntraits
    Phenotypes(obsID,featureID,nObs,nPheno,data) = 
        new(obsID,featureID,nObs,nPheno,data,nPheno)
end


mutable struct Equation
    from_layer_name::String 
    to_layer_name::String
    equation::String       
    covariate
    random
    activation_function
    partial_connect_structure
    ## method:
    method
    Pi
    estimatePi
    ## variance:
    G
    G_is_marker_variance
    df
    estimate_variance
    estimate_scale
    constraint
    #below function is used for genotypes data
    function Equation(;from_layer_name, 
                      to_layer_name, 
                      equation,
                      covariate=false,
                      random=false,
                      activation_function="linear",
                      partial_connect_structure=false,
                      method="BayesC", 
                      Pi=0.0,
                      estimatePi=true,
                      G=false,
                      G_is_marker_variance=false,
                      df = 4.0,
                      estimate_variance=true,
                      estimate_scale=false,
                      constraint=true)
        
        new(from_layer_name::String, 
            to_layer_name::String,
            equation::String,       
            covariate,
            random,
            activation_function,
            partial_connect_structure,
            method,
            Pi,
            estimatePi,
            G,
            G_is_marker_variance,
            df,
            estimate_variance,
            estimate_scale,
            constraint)    
    end
end


#    ## method:
#    method = "BayesC", Pi = 0.0, estimatePi = true, 
#    ## variance:
#    G_is_marker_variance = false, df = 4.0,
#    estimate_variance = true, estimate_scale = false,
#    constraint = false, #for multi-trait only, constraint=true means no genetic covariance among traits
