type Genotypes
  obsID       #now maybe string or int
  markerID
  nObs::Int64
  nMarkers::Int64
  alleleFreq::Array{Float64,2}
  sum2pq::Float64
  centered::Bool
  genotypes::Array{Float64,2}
  G::Float64    ##genetic vairance (prior)
  Genotypes(a,b,c,d,e,f,g,h)=new(a,b,c,d,e,f,g,h,0.0)
end

function addMarkers(mme::MME,file,G::Float64;separator=' ',header=true)
    mme.M   = readgenotypes(file;separator=separator,header=header,center=true)
    mme.M.G = G
end

#function addMarkers(mme::MME,df,G::Float64)
#    M = convert(Array,df)
#    mme.M = MarkerMatrix(M,G)
#end

include("tools.jl")
include("BayesC0.jl")
include("BayesC.jl")
include("BayesB.jl")
include("readgenotypes.jl")
