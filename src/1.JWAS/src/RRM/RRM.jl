#ID, time, obs
#a1, 1, y11
#a1, 3, y13
#a2, 1, y21
#a2, 2, y22
#a2, 3, y23
#a3, 2, y32
#
#yobs order: a1, a1, a2, a2, a2, a3
#yfull order: a1, a1, a1, a2, a2, a2, a3, a3, a3 (fill missing values)
function matrix_yfull_to_yobs(IDvec, yvec, timevec)
    times  = sort(unique(timevec))      # vector of whole timepoints sorted
    ntimes = length(times)           # number of timepoints measured in the data
    IDs    = unique(IDvec)           # vector of unique individuals
    nIDs   = length(IDs)             # number of individuals from data
    nobs   = length(yvec)            # number of observations from the data
    yobs = string.(IDvec) .* string.(timevec) # an extra vector as “ID x TimePoint” for all observations
    yfull = vcat([ID .* string.(times)  for ID in string.(IDs)]...) # an extra vector as “ID x TimePoint” for full y
    Is = Array{Int64}(undef, nobs)
    Js = Array{Int64}(undef, nobs)
    zeros_indicator = fill(true, length(yfull))
    whichrow = 1
    for i in 1:length(yfull)
        if yfull[i] in yobs
            Is[whichrow] = whichrow
            Js[whichrow] = i
            zeros_indicator[i] = false
            whichrow += 1
        end
    end
    matrix_yfull2yobs = sparse(Is, Js, ones(nobs),length(yobs),length(yfull))
    return matrix_yfull2yobs, zeros_indicator
end

#using Statistics, GSL
#Φ =  generatefullPhi(timevec)
function generatefullPhi(timevec, ncoeff=3)
    # timevec: a vector of time points for all observations
    # ncoeff : the number of polynomial coefficients to be estimated
    times = sort(unique(timevec)) # vector of whole timepoints sorted
    tmin  = minimum(times)
    tmax  = maximum(times)
    # standardized time points
    qi = 2*(times .- tmin)./(tmax .- tmin) .- 1
    Φ  = Matrix{Float64}(undef, length(times), ncoeff) # Phi of size ntimes x ncoeff
    for i in 1:ncoeff
        n = i - 1
        ## Given the normalized function from L. R. Schaeffer
        Φ[:,i] = sqrt((2*n+1)/2)*sf_legendre_Pl.(n,qi) # construct the matrix
    end
    return Φ
end

function get_mΦΦarray(Φ,T,mArray)
    ninds    = length(mArray[1])
    nmarkers = length(mArray)
    ntimes   = size(Φ,1)
    Φyobs    = T*repeat(Φ; outer=ninds) # whole Φ for all observations of size #observations X #coeff
    #### Construct an array (of lengh number of markers) of matrices (of order nCoeff x nCoeff)
    mΦΦarray = Array{Array{Float64,2}}(undef,nmarkers)
    for i in 1:nmarkers
        m  = mArray[i]
        mj = T*repeat(m;inner = ntimes) # vector of genotypes corresponding to yobs
        term = mj .* Φyobs
        mΦΦarray[i] = term'term
    end
    return mΦΦarray
end

function BayesABCRRM!(xArray,xpx,wArray,yfull,
                      betaArray,deltaArray,alphaArray,
                      vare,varEffects,BigPi,
                      Φ, zeros_indicator, mΦΦArray)
#For example of n individuals, 5 time points and 3 RR coefficient
#wArray is an array (of length number of timepoints)
#                    of vectors (of length number of individuals)
#betaArray,deltaArray,alphaArray is an array (of length number of RR coefficients)
#                                   of vectors (of length of number of markers)
#xArray is an array (of length number of markers)
#                                   of vector (of length number of individuals)
#varEffects is an array of the covariance matrix among RR coefficients
#         [G11 G12 G13
#         G21 G22 G23
#         G31 G32 G33]
#vare is a scallar (residual variance)

    nMarkers = length(xArray)
    nind     = length(xArray[1])
    nCoeff   = length(alphaArray)
    nTimes   = length(wArray)
    Rinv     = inv(vare)           #Do Not Use inv.(): elementwise
    Ginv     = inv.(varEffects)
    β        = zeros(typeof(betaArray[1][1]),nCoeff)
    newα     = zeros(typeof(alphaArray[1][1]),nCoeff)
    oldα     = zeros(typeof(alphaArray[1][1]),nCoeff)
    δ        = zeros(typeof(deltaArray[1][1]),nCoeff)
    # xw       = zeros(nTimes) # for rhs # construct directly
    xw_first   = zeros(nTimes)

    logDelta   = Dict{Any,Float64}() #of length number of labels
    probDelta  = Dict{Any,Float64}() #of length number of labels
    αcandidate = Dict()
    βcandidate = Dict()

    for marker=1:nMarkers
        x    = xArray[marker]
        for coeffi = 1:nCoeff
          oldα[coeffi]  = newα[coeffi] = alphaArray[coeffi][marker]
        end
        ## calculate sum(m_ij * Phi_i' w_ij)
        ## first part (see derivation)
        for timei = 1:nTimes
            xw_first[timei] = dot(x,wArray[timei])
        end
        mΦΦ = mΦΦArray[marker]
        xw = Φ'xw_first + mΦΦ*oldα
        #########################################

        for label in keys(BigPi)
            δj   = label
            Dj   = Diagonal(δj)
            lhs  = Dj'mΦΦ*Dj/vare + Ginv[marker]
            rhs  = Dj'xw/vare
            Σ    = inv(lhs)
            μ    = Σ *rhs
            logDelta[label]=-0.5*(log(det(lhs))-(rhs'μ)[1,1])+log(BigPi[label])
            beta = rand(MvNormal(μ, convert(Array,Symmetric(Σ))))
            βcandidate[label] = beta
            αcandidate[label] = Dj*beta
        end
        for i in keys(BigPi)
            if logDelta[i] == -Inf
                deno = Inf
            else
                deno =0.0
                for j in keys(BigPi)
                    deno += exp(logDelta[j]-logDelta[i])
                end
            end
            probDelta[i]=1/deno
        end
        whichlabel = rand((Categorical(collect(values(probDelta)))))
        newα = collect(values(αcandidate))[whichlabel]
        β    = collect(values(βcandidate))[whichlabel]
        δ    = collect(keys(BigPi))[whichlabel]

        ## adjust for locus j
        #when we change wArray, yfull is changed correspondingly
        #because we use pointer to generate wArray
        Φαold = Φ*oldα
        Φαnew = Φ*newα
        for timei = 1:nTimes
            BLAS.axpy!(Φαold[timei]-Φαnew[timei],x,wArray[timei])
        end
        yfull[zeros_indicator] .= 0         # renew wArray

        for coeffi = 1:nCoeff
            betaArray[coeffi][marker]       = β[coeffi]
            deltaArray[coeffi][marker]      = δ[coeffi]
            alphaArray[coeffi][marker]      = newα[coeffi]
        end
    end
end
