## assume ys is a nobs x 1 observations from the data
function Transform_DM(IDvec, ys, tpvec)
    time_points = sort(unique(tpvec))
    nt = length(time_points)
    IDs = unique(IDvec)
    nind = length(IDs)
    nobs = length(ys)
    DesignMatrix = zeros(Int8,nobs, nind*nt)
    observed_IDbyDays = string.(IDvec) .* string.(tpvec)
    all_IDbyDays = vcat([ID .* string.(time_points)  for ID in string.(IDs)]...)
    rows = 1
    for i in 1:nind*nt
        if all_IDbyDays[i] in observed_IDbyDays
            DesignMatrix[rows,i] = 1
            rows += 1
        end
    end
    return DesignMatrix
end
function GenrateFullPhi(tpvec, nCoeff) # nCoeff = number of polynomial coefficients to be estimated
    ti = sort(unique(tpvec))
    tmin = minimum(ti)
    tmax = maximum(ti)
    # standardized time points
    qi = 2*(ti .- tmin)./(tmax .- tmin) .- 1
    Phi = Matrix{Float64}(undef, length(ti), nCoeff)
    for i in 1:nCoeff
        n = i - 1
        ## Given the normalized function from L. R. Schaeffer
        Phi[:,i] = sqrt((2*n+1)/2)*sf_legendre_Pl.(n,qi)
    end
    return Phi
end

##Phenotypes adjusted for all effects (assuming zeros now)
ycorr          = vec(Matrix(mme.ySparse))
DM = Transform_DM(IDvec, ys, tpvec)
yfull = DM'ycorr
nObs = div(length(yfull),nTimes)
wArray         = Array{Array{Float64,1}}(undef,nTraits)
for ti = 1:nTimes
            startPosi              = (ti-1)*nObs  + 1
            ptr                    = pointer(yfull,startPosi)
            wArray[ti]         = unsafe_wrap(Array,ptr,nObs)
            betaArray[ti]      = copy(α[(ti-1)*nMarkers+1:ti*nMarkers])
            deltaArray[ti]     = ones(nMarkers)
            meandeltaArray[ti] = zeros(nMarkers)
            alphaArray[ti]         = copy(α[(ti-1)*nMarkers+1:ti*nMarkers])
            meanalphaArray[ti]     = zeros(nMarkers)
            meanalphaArray2[ti]    = zeros(nMarkers)
        end
    end
DM = Transform_DM(IDvec, ys, tpvec)
Φ =  GenrateFullPhi(tpvec, 3)

#MTBayesC requires the support for prior for delta for d is the set
#of all 2^n trait outcomes of dj:
function RRM!(xArray,xpx,wArray,betaArray,
              deltaArray,alphaArray,
              invR0,invG0,BigPi,
              Φ,vare,labels, DM)
#For example of n individuals, 5 time points and 3 RR coefficient
#wArray is an array (of length number of timepoints)
#                    of vectors (of length number of individuals)
#betaArray,deltaArray,alphaArray is an array (of length number of RR coefficients)
#                                   of vectors (of length of number of markers)
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
    Ginv     = invG0
    Rinv     = invR0

    β        = zeros(nCoeff)
    newα     = zeros(nCoeff)
    oldα     = zeros(nCoeff)
    δ        = zeros(nCoeff)
    xw       = zeros(nTimes) #for rhs
    α        = Array{Array{Float64,1}}(undef,nlabels) ###
    beta     = Array{Array{Float64,1}}(undef,nlabels) ###

    T = DM'DM

    for i in 1:length(labels)
      δi = labels[i]
      D  = diagm(δi)
      RinvLhs[i] = D*Rinv*D #split better
      RinvRhs[i] = Rinv*D
    end

    for marker=1:nMarkers
        x    = xArray[marker]
        for coeffi = 1:nCoeff
            β[coeffi]  = betaArray[coeffi][marker]  # we don't need this
         oldα[coeffi]  = newα[coeffi] = alphaArray[coeffi][marker]
            δ[coeffi]  = deltaArray[coeffi][marker] # we don't need this
        end

        oldΦα =  Φ*oldα

        for timei = 1:nTimes
            if
            xw[timei]  = dot(x,wArray[timei])+xpx[marker]*oldΦα[timei] #x_{ij}*w_{ij} #WRONG!!!!
        end

        for label in 1:length(labels)
            δi   = labels[label]
            Dj   = Diagonal(δj)
            temp = Φ*Dj
            lhs  = temp'temp * xpx[marker]/vare + Ginv
            rhs  = Dj'Φ'xw/vare
            Σ    = inv(lhs)
            μ    = Σ *rhs
            logDelta[label]=-0.5*(log(det(lhs))-(rhs*μ)[1,1])+log(BigPi[label])

            beta[label] = rand(MvNormal(μ, Σ))
            α[label]  = Dj*beta[label]
        end
        for label in 1:length(labels)
          deno =0.0
          for i in 1:length(labels)
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
        Φαold = Φ*oldα #individual phi is needed
        Φαnew = Φ*newα
        for timei = 1:nTimes
            BLAS.axpy!(Φαold[timei]-Φαnew[timei],x,wArray[nTimes])
        end

        ###################################################
        # renew wArray
        yfull = [wArray[i][j] for i=1:nTimes, j = 1:nind] #[wArray[1] wArray[2] wArray[3]]
        yfull = reshape(yfull, 1, :)
        # updated yfull and wArray
        yfull = T*yfull
        startPosi          = [(ti-1)*nind  + 1 for ti in 1:nTimes]
        ptr                = pointer(yfull,startPosi)
        wArray[ti]         = unsafe_wrap(Array,ptr,nind)
    end
end
