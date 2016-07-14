type Genotypes
  #obsID::Array{UTF8String,1}    #row ID of genotypes
  obsID  #now maybe string or int
  markerID
  nObs::Int64
  nMarkers::Int64
  alleleFreq::Array{Float64,2}
  sum2pq::Float64
  centered::Bool
  genotypes::Array{Float64,2}
  G::Array{Float64,2} ##marker effects covariance matrix
  Genotypes(a,b,c,d,e,f,g,h)=new(a,b,c,d,e,f,g,h,zeros(2,2))
end

function addMarkers(mme::MME,file,G::Array{Float64,2};separator=' ',
                    header=true,G_is_genetic_variance=G_is_genetic_variance)
    mme.M   = readgenotypes(file;separator=separator,header=header,center=true)
    if G_is_genetic_variance==true
      mme.M.G = G/mme.M.sum2pq
    else
      mme.M.G = G
    end
end

include("tools.jl")
include("MTBayesC0.jl")
include("MTBayesC.jl")
include("MTBayesCC.jl")
include("MTBayesB.jl")
include("readgenotypes.jl")
