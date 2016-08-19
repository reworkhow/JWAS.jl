include("tools.jl")
include("genotype_type.jl")
include("readgenotypes.jl")
include("BayesianAlphabet/BayesC0.jl")
include("BayesianAlphabet/BayesC.jl")
include("BayesianAlphabet/BayesB.jl")
include("BayesianAlphabet/MTBayesC0.jl")
include("BayesianAlphabet/MTBayesC.jl")
include("BayesianAlphabet/MTBayesCC.jl")
include("BayesianAlphabet/MTBayesB.jl")


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


function addMarkers(mme::MME,file,G;
                    separator=' ',header=true,
                    G_is_marker_variance=false)
    mme.M   = readgenotypes(file;separator=separator,header=header,center=true)
    mme.M.G = G
    mme.M.G_is_marker_variance = G_is_marker_variance
end
