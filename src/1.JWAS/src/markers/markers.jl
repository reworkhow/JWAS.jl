include("tools4genotypes.jl")
include("readgenotypes.jl")
include("BayesianAlphabet/BayesC0.jl")
include("BayesianAlphabet/BayesC.jl")
include("BayesianAlphabet/BayesB.jl")
include("BayesianAlphabet/MTBayesC.jl")
include("BayesianAlphabet/MTBayesCC.jl")
include("BayesianAlphabet/MTBayesB.jl")
include("BayesianAlphabet/MTBayesC0L.jl")

################################################################################
#Set Marker Effcets Hyperparameters: Variances and Pi
################################################################################
function set_marker_hyperparameters_variances_and_pi(mme::MME,Pi,methods)
  #multi-trait (Pi)
  if mme.nModels !=1 && Pi==0.0
      printstyled("Pi (Π) is not provided.\n",bold=false)
      printstyled("Pi (Π) is generated assuming all markers have effects on all traits.\n",bold=false)
      mykey=Array{Float64}(undef,0)
      ntraits=mme.nModels
      Pi=Dict{Array{Float64,1},Float64}()
      #for i in [ bin(n,ntraits) for n in 0:2^ntraits-1 ] `bin(n, pad)` is deprecated, use `string(n, base=2, pad=pad)
      for i in [ string(n,base=2,pad=ntraits) for n in 0:2^ntraits-1 ]
        Pi[parse.(Float64, split(i,""))]=0.0
      end
      Pi[ones(ntraits)]=1.0
  end
  #marker effect variances
  if mme.M.G == false && methods!="GBLUP"
      genetic2marker(mme.M,Pi)
      println()
      if mme.nModels != 1
        if !isposdef(mme.M.G) #also work for scalar
          error("Marker effects covariance matrix is not postive definite! Please modify the argument: Pi.")
        end
        println("The prior for marker effects covariance matrix is calculated from genetic covariance matrix and Π.")
        println("The mean of the prior for the marker effects covariance matrix is:")
        Base.print_matrix(stdout,round.(mme.M.G,digits=6))
      else
        if !isposdef(mme.M.G) #positive scalar (>0)
          error("Marker effects variance is negative!")
        end
        println("The prior for marker effects variance is calculated from the genetic variance and π.")
        print("The mean of the prior for the marker effects variance is: ")
        print(round.(mme.M.G,digits=6))
      end
      print("\n\n\n")
  end
  return Pi
end

function genetic2marker(M::Genotypes,Pi::Dict)
  nTraits = size(M.genetic_variance,1)
  denom   = zeros(nTraits,nTraits)
  for i in 1:nTraits
    for j in i:nTraits
      #pi_selected = filter((k,v)->k[i]==1.0 && k[j]==1.0,Pi)
      pi_selected = filter(d->d.first[i]==1.0 && d.first[j]==1.0,Pi)

      denom[i,j] = M.sum2pq*sum(values(pi_selected))
      denom[j,i] = denom[i,j]
    end
  end
  M.G = M.genetic_variance ./ denom
end

function genetic2marker(M::Genotypes,π::Float64)
    M.G = M.genetic_variance/((1-π)*M.sum2pq)
end
