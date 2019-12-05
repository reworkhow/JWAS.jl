using DataFrames, LinearAlgebra, Statistics, Distributions, Random, DelimitedFiles, GSL;
##Phenotypes adjusted for all effects (assuming zeros now)
ycorr  = vec(Matrix(mme.ySparse))
function Transform_DM(IDvec, ys, tpvec)
    # input:
    # IDvec - a vector of ID for all the observed data
    # ys: a vector of observations from the data
    # tpvec: a vector of TimePoints for all observed data
    time_points = sort(unique(tpvec)) # vector of whole timepoints sorted
    nt = length(time_points) # number of timepoints measured in the data
    IDs = unique(IDvec) # vector of unique individuals
    nind = length(IDs) # number of individuals having records
    nobs = length(ys) # number of observations from the data
    observed_IDbyDays = string.(IDvec) .* string.(tpvec) # an extra vector as “ID x TimePoint” for all observations
    all_IDbyDays = vcat([ID .* string.(time_points)  for ID in string.(IDs)]...) # an extra vector as “ID x TimePoint” for full y
    Is = Array{Int64}(undef, nobs)
    Js = Array{Int64}(undef, nobs)
    Vs = Array{Float64}(undef, nobs)
    indicator = Array{Float64, 1}(undef, nind*nt)
    rows = 1
    for i in 1:nind*nt
        if all_IDbyDays[i] in observed_IDbyDays
            Is[rows] = rows
            Js[rows] = i
            Vs[rows] = 1 # put 1 at each column corresponding to the observations
            indicator[i] = false
            rows += 1
        else
            indicator[i] = true
        end
    end
    DesignMatrix = sparse(Is, Js, Vs,nobs,nind*nt)
    return DesignMatrix, indicator
end


function GenrateFullPhi(tpvec, nCoeff)
    # input:
    # tpvec: a vector of TimePoints for all observed data
    # nCoeff = number of polynomial coefficients to be estimated
    ti = sort(unique(tpvec)) # vector of whole timepoints sorted
    tmin = minimum(ti) # minimal timepoint
    tmax = maximum(ti) # maximal timepoint
    # standardized time points
    qi = 2*(ti .- tmin)./(tmax .- tmin) .- 1
    Phi = Matrix{Float64}(undef, length(ti), nCoeff) # a matrix for the whole Phi corresponding to whole timepoints
    for i in 1:nCoeff
        n = i - 1
        ## Given the normalized function from L. R. Schaeffer
        Phi[:,i] = sqrt((2*n+1)/2)*sf_legendre_Pl.(n,qi) # construct the matrix
    end
    return Phi
end


#Φ =  GenrateFullPhi(tpvec, nCoeff) # complete Φ for all time points

function get_mΦΦArray(Φ,T,xArray)
    nind = length(xArray[1])
    nMarkers = length(xArray)
    nTimes = size(Φ)[1]
    Φw   = T*repeat(Φ; outer=nind) # whole Φ for observed data
    #### Construct an array (of lengh number of markers) of matrices (of order nCoeff x nCoeff)
    mΦΦArray = Array{Array{Float64,2}}(undef,nMarkers)
    for marker in 1:nMarkers
        x = xArray[marker]
        Mj = T*repeat(x;inner = nTimes) # vector of genetypes corresponding to ycor (observed y)
        term = Mj .* Φw
        mΦΦArray[marker] = term'term
    end
    return mΦΦArray
end

#MTBayesC requires the support for prior for delta for d is the set
#of all 2^n trait outcomes of dj:
function BayesABCRRM!(xArray,xpx,wArray,
                      betaArray,deltaArray,alphaArray,
                      vare,varEffects,BigPi,labels,
                      Φ, whichzeros, mΦΦArray)
#For example of n individuals, 5 time points and 3 RR coefficient
#wArray is an array (of length number of timepoints)
#                    of vectors (of length number of individuals)
#betaArray,deltaArray,alphaArray is an array (of length number of RR coefficients)
#                                   of vectors (of length of number of markers)
#xArray is an array (of length number of markers)
#                                   of vector (of length number of individuals)
#invG0 is the inverse of the covariance matrix among RR coefficients
#        [G11 G12 G13
#         G21 G22 G23
#         G31 G32 G33]
#invR0 is the inverse of the covariance matrix among time pointes
#        Diag(σ2,σ2,σ2,σ2,σ2)

    nMarkers = length(xArray)
    nind     = length(xArray[1])
    nCoeff   = length(alphaArray)
    nTimes   = length(wArray)
    nlabels  = length(labels)
    Rinv     = inv(vare)           #Do Not Use inv.(): elementwise
    Ginv     = inv.(varEffects)
    β        = zeros(typeof(betaArray[1][1]),nCoeff)
    newα     = zeros(typeof(alphaArray[1][1]),nCoeff)
    oldα     = zeros(typeof(alphaArray[1][1]),nCoeff)
    δ        = zeros(typeof(deltaArray[1][1]),nCoeff)
    # xw       = zeros(nTimes) # for rhs # construct directly
    xw_first = zeros(nTimes)
    αcandidate = Array{Array{typeof(alphaArray[1][1]),1}}(undef,nlabels) ###
    βcandidate = Array{Array{typeof(betaArray[1][1]),1}}(undef,nlabels) ### an array (of length number of labels)
#                                                                  of vector (of length number of RR coefficients)
    for marker=1:nMarkers
        x    = xArray[marker]
        for coeffi = 1:nCoeff
            β[coeffi]  = betaArray[coeffi][marker]  # we don't need this
         oldα[coeffi]  = newα[coeffi] = alphaArray[coeffi][marker]
            δ[coeffi]  = deltaArray[coeffi][marker] # we don't need this
        end
        ## calculate sum(m_ij * Phi_i' w_ij)
        ## first part (see derivation)
        for timei = 1:nTimes
            xw_first[timei] = dot(x,wArray[timei])
        end
        mΦΦ = mΦΦArray[marker]
        xw = Φ'xw_first + mΦΦ*oldα
        #########################################

        for label in 1:nlabels
            δj = labels[label]
            Dj   = Diagonal(δj)
            lhs  = Dj'mΦΦ*Dj/vare + Ginv
            rhs  = Dj'xw/vare
            Σ    = inv(lhs)
            μ    = Σ *rhs
            logDelta[label]=-0.5*(log(det(lhs))-(rhs*μ)[1,1])+log(BigPi[label])
            beta[label] = rand(MvNormal(μ, Σ))
            α[label]  = Dj*beta[label]
        end
        for label in 1:nlabels
          deno =0.0
          for i in 1:nlabels
            deno += exp(logDelta[i]-logDelta[label])
          end
          probDelta[label]=1/deno
        end

        whichlabel = rand((Categorical(probDelta)))
        newα = α[whichlabel]
        β = beta[whichlabel]
        δ = labels[whichlabel]

        # adjust for locus j
        for coeffi = 1:nCoeff
            betaArray[coeffi][marker]       = β[coeffi]
            deltaArray[coeffi][marker]      = δ[coeffi]
            alphaArray[coeffi][marker]      = newα[coeffi]
        end

        #### when we change wArray, we change yfull correspondingly because we use pointer to generate wArray
        ## indicator for yfull
        Φαold = Φ*oldα
        Φαnew = Φ*newα
        for timei = 1:nTimes
            BLAS.axpy!(Φαold[timei]-Φαnew[timei],x,wArray[nTimes])
        end
        # renew wArray
        yfull[whichzeros] .= 0
    end
end
