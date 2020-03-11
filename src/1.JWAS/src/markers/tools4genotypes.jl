#the function below using pointers is deprecated
function get_column(X,j)
    nrow,ncol = size(X)
    if j>ncol||j<0
        error("column number is wrong!")
    end
    indx = 1 + (j-1)*nrow
    ptr = pointer(X,indx)
    unsafe_wrap(Array,ptr,nrow) #pointer_to_array(ptr,nrow) #deprected in Julia 0.5
end

function get_column_ref(X)
    ncol = size(X)[2]
    xArray = Array{Union{Array{Float64,1},Array{Float32,1}}}(undef,ncol)
    for i=1:ncol
        xArray[i] = get_column(X,i)
    end
    return xArray
end

function center!(X)
    nrow,ncol = size(X)
    colMeans = mean(X,dims=1)
    BLAS.axpy!(-1,ones(nrow)*colMeans,X)
    return colMeans
end

function center(X)
    nrow,ncol = size(X)
    colMeans = mean(X,dims=1)
    return colMeans
end

function getXpRinvX(X, Rinv)
    XpRinvX = [((X[:,i].*Rinv)'X[:,i]) for i=1:size(X,2)]
    return XpRinvX
end

function getXpRinvX(X)
    XpRinvX = [dot(X[:,i],X[:,i]) for i=1:size(X,2)]
    return XpRinvX
end

mutable struct GibbsMats
    X::Union{Array{Float64,2},Array{Float32,2}}
    nrows::Int64
    ncols::Int64
    xArray::Array{Union{Array{Float64,1},Array{Float32,1}},1}
    xRinvArray::Array{Union{Array{Float64,1},Array{Float32,1}},1}
    xpRinvx::Union{Array{Float64,1},Array{Float32,1}}
    function GibbsMats(X::Union{Array{Float64,2},Array{Float32,2}},Rinv)
        nrows,ncols = size(X)
        xArray      = get_column_ref(X)
        xpRinvx     = getXpRinvX(X,Rinv)
        if Rinv == ones(length(Rinv))
            xRinvArray = xArray  #avoid using extra memory for xRinvArray
        else
            xRinvArray = [x.*Rinv for x in xArray]
        end
        new(X,nrows,ncols,xArray,xRinvArray,xpRinvx)
    end
end

################################################################################
#align genotypes with phenotyps given IDs #MAY USE TOO MUCH MEMORY
################################################################################
#input: mme.M.genotypes (1) complete genomic data  -> all genotyped individuals
#                       (2)* incomplete genomic data -> all pedigree individuals
#output:mme.M.genotypes (1) subset
#
# mme.M.genotypes is modified to same IDs in phenotypes for coding convenience
#
#Note that Aligning genotypes with phenotypes and output_ID in single-step analysis
#(incomplete genomic data) has been moved to SSBR.jl file.
function align_genotypes(mme::MME,output_heritability=false,single_step_analysis=false)
    if mme.output_ID != 0 && single_step_analysis==false
        if mme.output_ID != mme.M.obsID
            Zo  = map(Float32,mkmat_incidence_factor(mme.output_ID,mme.M.obsID))
            mme.output_genotypes =  Zo*mme.M.genotypes
        else
            mme.output_genotypes =  mme.M.genotypes
        end
    end
    #***************************************************************************
    #Align genotypes with phenotypes
    #Change mme.M.genotypes to phenotyped individuals to accommodate
    #individuals with repeated records or individuals without records
    #
    #**********CENTERING?*******************************************************
    if mme.obsID != mme.M.obsID && single_step_analysis==false
        Z  = map(Float32,mkmat_incidence_factor(mme.obsID,mme.M.obsID))
        mme.M.genotypes = Z*mme.M.genotypes
        mme.M.obsID     = mme.obsID
        mme.M.nObs      = length(mme.M.obsID)
    end
end

#get an incidence matrix Z to reorder uID to yID by yID = Z*uID
function mkmat_incidence_factor(yID,uID)
    Z = spzeros(length(yID),length(uID))

    uIDdict = Dict()
    for (index,id) in enumerate(uID)
        uIDdict[id]=index
    end

    rowi = 1
    for id in yID
        if haskey(uIDdict,id)
            index = uIDdict[id]
        else
            error(id, " is not found!")
        end
        Z[rowi,index]=1
        rowi = rowi+1
    end
    return Z
end

################################################################################
#Set Marker Effcets Hyperparameters: Variances and Pi
################################################################################
function set_marker_hyperparameters_variances_and_pi(mme::MME,Pi,methods)
  #multi-trait (Pi)
  if mme.nModels !=1 && Pi==0.0 && !(methods in ["RR-BLUP","BayesL","BayesA"])
      println()
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
  if mme.nModels == 1
      mme.M.scale = mme.M.G*(mme.df.marker-2)/mme.df.marker #scale parameter for marker effect variance
  else
      mme.M.scale = mme.M.G*(mme.df.marker - mme.nModels - 1)
  end
  return Pi
end

function genetic2marker(M::Genotypes,Pi::Dict)
  nTraits = size(M.genetic_variance,1)
  denom   = zeros(nTraits,nTraits)
  for i in 1:nTraits
    for j in i:nTraits
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
