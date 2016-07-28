type Genotypes
  obsID       #now maybe string or int
  markerID
  nObs::Int64
  nMarkers::Int64
  alleleFreq::Array{Float64,2}
  sum2pq::Float64
  centered::Bool
  genotypes::Array{Float64,2}
  G::Float64
  G_is_marker_variance::Bool
  Genotypes(a1,a2,a3,a4,a5,a6,a7,a8)=new(a1,a2,a3,a4,a5,a6,a7,a8,0.0,false)
end

function addMarkers(mme::MME,file,G::Float64;
                    separator=' ',header=true,
                    G_is_marker_variance=false)
    mme.M   = readgenotypes(file;separator=separator,header=header,center=true)
    mme.M.G = G
    mme.M.G_is_marker_variance = G_is_marker_variance
end

include("tools.jl")
include("BayesC0.jl")
include("BayesC.jl")
include("BayesB.jl")
include("readgenotypes.jl")
