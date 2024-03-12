# Reference:
#Cheng et al. 2018 Genomic prediction from multiple-trait bayesian regression methods using mixture priors. Genetics209:  89–103.
# Gibbs sampler I: only one of the t indicator lables is sampled at a time.
# Gibbs sampler II: all indicator labels are sampled jointly, can be used for the restrictive model assuming any particular locus affects all traits or none of them.

function MTBayesABC!(genotypes,ycorr_array,vare,locus_effect_variances,nModels)
    if length(genotypes.π) == 2^nModels  #Gibbs sampler I
        MTBayesABC!(genotypes.mArray,genotypes.mRinvArray,genotypes.mpRinvm,
                    ycorr_array,genotypes.β,genotypes.δ,genotypes.α,vare,
                    locus_effect_variances,genotypes.π)
    else  #Gibbs sampler II
        MTBayesABC_samplerII!(genotypes.mArray,genotypes.mRinvArray,genotypes.mpRinvm,
                    ycorr_array,genotypes.β,genotypes.δ,genotypes.α,vare,
                    locus_effect_variances,genotypes.π)
    end
end

#Gibbs sampler I
function MTBayesABC!(xArray,xRinvArray,xpRinvx,
                     wArray,betaArray,
                     deltaArray,
                     alphaArray,
                     vare,varEffects,
                     BigPi)
    nMarkers = length(xArray)
    ntraits  = length(alphaArray)

    Rinv     = inv(vare) #Do Not Use inv.(): elementwise inversion
    Ginv     = inv.(varEffects)

    β        = zeros(typeof(betaArray[1][1]),ntraits)
    newα     = zeros(typeof(alphaArray[1][1]),ntraits)
    oldα     = zeros(typeof(alphaArray[1][1]),ntraits)
    δ        = zeros(typeof(deltaArray[1][1]),ntraits)
    w        = zeros(typeof(wArray[1][1]),ntraits) #for rhs

    for marker=1:nMarkers
        x, xRinv = xArray[marker], xRinvArray[marker]

        for trait = 1:ntraits
            β[trait]  = betaArray[trait][marker]
         oldα[trait]  = newα[trait] = alphaArray[trait][marker]
            δ[trait]  = deltaArray[trait][marker]
            w[trait]  = dot(xRinv,wArray[trait])+xpRinvx[marker]*oldα[trait]
        end

        for k=1:ntraits
            Ginv11 = Ginv[marker][k,k]
            nok    = deleteat!(collect(1:ntraits),k)
            Ginv12 = Ginv[marker][k,nok]
            C11    = Ginv11+Rinv[k,k]*xpRinvx[marker]
            C12    = Ginv12+xpRinvx[marker]*Matrix(Diagonal(δ[nok]))*Rinv[k,nok]

            invLhs0  = 1/Ginv11
            rhs0     = - Ginv12'β[nok]
            gHat0    = (rhs0*invLhs0)[1,1]
            invLhs1  = 1/C11
            rhs1     = w'*Rinv[:,k]-C12'β[nok]
            gHat1    = (rhs1*invLhs1)[1,1]

            d0 = copy(δ)
            d1 = copy(δ)
            d0[k] = 0.0
            d1[k] = 1.0

            logDelta0  = -0.5*(log(Ginv11)- gHat0^2*Ginv11) + log(BigPi[d0]) #logPi
            logDelta1  = -0.5*(log(C11)-gHat1^2*C11) + log(BigPi[d1]) #logPiComp

            probDelta1 =  1.0/(1.0+exp(logDelta0-logDelta1))
            if(rand()<probDelta1)
                δ[k] = 1
                β[k] = newα[k] = gHat1 + randn()*sqrt(invLhs1)
                BLAS.axpy!(oldα[k]-newα[k],x,wArray[k])
            else
                β[k] = gHat0 + randn()*sqrt(invLhs0)
                δ[k] = 0
                newα[k] = 0
                if oldα[k] != 0
                    BLAS.axpy!(oldα[k],x,wArray[k])
                end
            end
        end
        for trait = 1:ntraits
            betaArray[trait][marker]       = β[trait]
            deltaArray[trait][marker]      = δ[trait]
            alphaArray[trait][marker]      = newα[trait]
        end
    end
end


#Gibbs sampler II
function MTBayesABC_samplerII!(xArray,
                    xRinvArray,
                    xpRinvx,
                    wArray, #ycorr
                    betaArray,
                    deltaArray,
                    alphaArray,
                    vare, #mme.R.val, t-by-t
                    varEffects, # vector of length #SNPs, each element is a t-by-t covariance matrix
                    BigPi) #genotypes.π

    nMarkers = length(xArray)
    ntraits  = length(alphaArray)

    Rinv     = inv(vare) #inv(mme.R.val)
    Ginv     = inv.(varEffects)

    β        = zeros(typeof(betaArray[1][1]),ntraits)
    newα     = zeros(typeof(alphaArray[1][1]),ntraits)
    oldα     = zeros(typeof(alphaArray[1][1]),ntraits)
    δ        = zeros(typeof(deltaArray[1][1]),ntraits)
    w        = zeros(typeof(wArray[1][1]),ntraits) #for rhs

    nlable    = length(BigPi) #e.g.,length(Dict([1.0; 1.0]=>0.3,[0.0; 0.0]=>0.7))=2
    probDelta = Array{Float64}(undef, nlable)
    logDelta  = Array{Float64}(undef, nlable)
    βeta      = Array{Array{Float64,1}}(undef, nlable)
    RinvLhs   = Array{Array{Float64,2}}(undef, nlable) # D*inv(R)*D
    RinvRhs   = Array{Array{Float64,2}}(undef, nlable) # inv(R)*D


    iloci=1
    for i in keys(BigPi)
        D  = diagm(i)
        RinvLhs[iloci] = D*Rinv*D # D*inv(R)*D
        RinvRhs[iloci] = Rinv*D   # inv(R)*D
        iloci += 1
    end

    for marker=1:nMarkers
        x, xRinv = xArray[marker], xRinvArray[marker]

        for trait = 1:ntraits
            β[trait]  = betaArray[trait][marker]
         oldα[trait]  = newα[trait] = alphaArray[trait][marker]
            δ[trait]  = deltaArray[trait][marker]
            w[trait]  = dot(xRinv,wArray[trait])+xpRinvx[marker]*oldα[trait] #w=xj'(ycorr+xj) t-by-1
        end

        #full conditional distribution of β
        stdnorm = randn(ntraits)
        iloci=1
        for pi_t in values(BigPi)
            lhs       = RinvLhs[iloci]*xpRinvx[marker]+Ginv[marker]  #C: t-by-t
            rhs       = RinvRhs[iloci]'w  #t-by-1
            invLhs    = inv(lhs)
            invLhsC   = cholesky(Hermitian(invLhs)).L # L, where LL'=invLhs Q:why Hermitian
            gHat      = invLhs*rhs  #t-by-1
            logDelta[iloci] = -0.5*(log(det(lhs))-rhs'gHat)+log(pi_t)
            βeta[iloci]     = gHat + invLhsC*stdnorm  #var(Lz)=LL'=invLhs=inv(C)
            iloci += 1
        end

        #marginal full conditional probability of δ
        for l in 1:nlable
          denominator_l = 0.0  #denominator
          for i in 1:nlable
            denominator_l += exp(logDelta[i]-logDelta[l])
          end
          probDelta[l]=1/denominator_l
        end



        #choose label
        whichlabel = rand(Categorical(probDelta)) #e.g., 1 means the 1st label
        δ           = copy(collect(keys(BigPi))[whichlabel]) ##copy required,otherwise modifiying keys(BigPi)
        β           = βeta[whichlabel]
        newα        = diagm(δ)*β


        # adjust for locus j
        for trait in 1:ntraits
            BLAS.axpy!(oldα[trait]-newα[trait],x,wArray[trait])  #update wArray (ycorr)
            betaArray[trait][marker]       = β[trait]
            deltaArray[trait][marker]      = δ[trait]
            alphaArray[trait][marker]      = newα[trait]
        end

    end
end

#block
function MTBayesABC_block!(genotypes,ycorr_array,vare,locus_effect_variances)
    MTBayesABC_block!(genotypes.MArray,genotypes.MRinvArray,genotypes.mpRinvm,
                genotypes.genotypes,genotypes.MpRinvM,
                ycorr_array,genotypes.β,genotypes.δ,genotypes.α,vare,
                locus_effect_variances,genotypes.π)
end

function MTBayesABC_block!(XArray,XRinvArray,xpRinvx,
                     X, XpRinvX,
                     wArray,
                     betaArray,deltaArray,alphaArray,
                     vare,varEffects,
                     BigPi)
    nMarkers = length(xpRinvx)
    ntraits  = length(alphaArray)

    Rinv     = inv(vare) #Do Not Use inv.(): elementwise inversion
    Ginv     = inv.(varEffects) # vector of length p, each element is a t-by-t covariance matrix

    #vector of length t, since we only focus one SNP at a time
    β        = zeros(typeof(betaArray[1][1]),ntraits)
    newα     = zeros(typeof(alphaArray[1][1]),ntraits)
    oldα     = zeros(typeof(alphaArray[1][1]),ntraits)
    δ        = zeros(typeof(deltaArray[1][1]),ntraits)
    w        = zeros(typeof(wArray[1][1]),ntraits) #for rhs

    XpRinvX       = XpRinvX
    nblocks       = length(XpRinvX)
    start_pos     = 0
    for i in 1:nblocks
        XpRinvycorr = [XRinvArray[i]*wArray[trait] for trait = 1:ntraits]
        oldalphaArray = deepcopy(alphaArray)
        block_size  = size(XpRinvX[i],1)
        nreps       = block_size #user-defined nreps=block_size, outer_niter =niter/block_size
        for reps = 1:nreps
            for j = 1:block_size #additional code to save all sampled αs is needed
                locus_j = start_pos + j
                for trait = 1:ntraits
                    β[trait]  = betaArray[trait][locus_j]
                 oldα[trait]  = newα[trait] = alphaArray[trait][locus_j]
                    δ[trait]  = deltaArray[trait][locus_j]
                    w[trait]  = XpRinvycorr[trait][j]+xpRinvx[locus_j]*oldα[trait]
                end
                for k=1:ntraits
                    Ginv11 = Ginv[locus_j][k,k]
                    nok    = deleteat!(collect(1:ntraits),k)
                    Ginv12 = Ginv[locus_j][k,nok]
                    C11    = Ginv11+Rinv[k,k]*xpRinvx[locus_j]
                    C12    = Ginv12+xpRinvx[locus_j]*Matrix(Diagonal(δ[nok]))*Rinv[k,nok]

                    invLhs0  = 1/Ginv11
                    rhs0     = - Ginv12'β[nok]
                    gHat0    = (rhs0*invLhs0)[1,1]
                    invLhs1  = 1/C11
                    rhs1     = w'*Rinv[:,k]-C12'β[nok]
                    gHat1    = (rhs1*invLhs1)[1,1]

                    d0 = copy(δ)
                    d1 = copy(δ)
                    d0[k] = 0.0
                    d1[k] = 1.0

                    logDelta0  = -0.5*(log(Ginv11)- gHat0^2*Ginv11) + log(BigPi[d0]) #logPi
                    logDelta1  = -0.5*(log(C11)-gHat1^2*C11) + log(BigPi[d1]) #logPiComp

                    probDelta1 =  1.0/(1.0+exp(logDelta0-logDelta1))
                    if(rand()<probDelta1)
                        δ[k] = 1
                        β[k] = newα[k] = gHat1 + randn()*sqrt(invLhs1)
                        BLAS.axpy!(oldα[k]-newα[k],view(XpRinvX[i],:,j),XpRinvycorr[k])
                    else
                        β[k] = gHat0 + randn()*sqrt(invLhs0)
                        δ[k] = 0
                        newα[k] = 0
                        if oldα[k] != 0
                            BLAS.axpy!(oldα[k],view(XpRinvX[i],:,j),XpRinvycorr[k])
                        end
                    end
                end
                for trait = 1:ntraits
                    betaArray[trait][locus_j]       = β[trait]
                    deltaArray[trait][locus_j]      = δ[trait]
                    alphaArray[trait][locus_j]      = newα[trait]
                end
            end
        end
        for trait = 1:ntraits
            wArray[trait][:]  = wArray[trait] + XArray[i]*(oldalphaArray[trait]-alphaArray[trait])[(start_pos+1):(start_pos+block_size)]
        end
        start_pos += block_size
    end
end
