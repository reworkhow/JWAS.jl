# Reference:
#Cheng et al. 2018 Genomic prediction from multiple-trait bayesian regression methods using mixture priors. Genetics209:  89–103.
# Gibbs sampler I: only one of the t indicator lables is sampled at a time.
# Gibbs sampler II: all indicator labels are sampled jointly, can be used for the restrictive model assuming any particular locus affects all traits or none of them.

struct GlobalPiPrior{P}
    big_pi::P
end

struct MarkerSpecificPiPrior{T}
    snp_pi::T
end

@inline log_prior(prior::GlobalPiPrior, marker::Integer, state::AbstractVector{<:Real}) =
    log(prior.big_pi[state])

@inline log_prior(prior::MarkerSpecificPiPrior, marker::Integer, state::AbstractVector{<:Real}) =
    log(prior.snp_pi[marker, annotated_bayesc_mt_state_index(state)])

function mt_bayesc_sampler_mode(genotypes, nModels)
    sampler = getproperty(genotypes, :multi_trait_sampler)
    sampler == :auto && return length(genotypes.π) == 2^nModels ? :I : :II
    sampler in (:I, :II) || error("multi_trait_sampler must be one of :auto, :I, or :II.")
    return sampler
end

function mt_bayesc_block_context(genotypes, nModels)
    prior = if has_marker_annotations(genotypes) && genotypes.method == "BayesC"
        MarkerSpecificPiPrior(genotypes.annotations.snp_pi)
    else
        GlobalPiPrior(genotypes.π)
    end
    sampler_mode = mt_bayesc_sampler_mode(genotypes, nModels)
    return prior, sampler_mode
end

function MTBayesABC!(genotypes,ycorr_array,vare,locus_effect_variances,nModels)
    prior = if has_marker_annotations(genotypes) && genotypes.method == "BayesC" && genotypes.ntraits == 2
        MarkerSpecificPiPrior(genotypes.annotations.snp_pi)
    else
        GlobalPiPrior(genotypes.π)
    end
    sampler_mode = mt_bayesc_sampler_mode(genotypes, nModels)
    if sampler_mode == :I  # Gibbs sampler I
        _MTBayesABC_samplerI!(genotypes.mArray,genotypes.mRinvArray,genotypes.mpRinvm,
                             ycorr_array,genotypes.β,genotypes.δ,genotypes.α,vare,
                             locus_effect_variances,prior)
    else  # Gibbs sampler II
        state_labels = collect(keys(genotypes.π))
        _MTBayesABC_samplerII!(genotypes.mArray,genotypes.mRinvArray,genotypes.mpRinvm,
                               ycorr_array,genotypes.β,genotypes.δ,genotypes.α,vare,
                               locus_effect_variances,state_labels,prior)
    end
end

#Gibbs sampler I
function _MTBayesABC_samplerI!(xArray,xRinvArray,xpRinvx,
                               wArray,betaArray,
                               deltaArray,
                               alphaArray,
                               vare,varEffects,
                               prior)
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

            logDelta0  = -0.5*(log(Ginv11)- gHat0^2*Ginv11) + log_prior(prior, marker, d0)
            logDelta1  = -0.5*(log(C11)-gHat1^2*C11) + log_prior(prior, marker, d1)

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

function _MTBayesABC_samplerII!(xArray,
                                xRinvArray,
                                xpRinvx,
                                wArray,
                                betaArray,
                                deltaArray,
                                alphaArray,
                                vare,
                                varEffects,
                                state_labels,
                                prior)
    nMarkers = length(xArray)
    ntraits  = length(alphaArray)

    Rinv     = inv(vare) #inv(mme.R.val)
    Ginv     = inv.(varEffects)

    β        = zeros(typeof(betaArray[1][1]),ntraits)
    newα     = zeros(typeof(alphaArray[1][1]),ntraits)
    oldα     = zeros(typeof(alphaArray[1][1]),ntraits)
    δ        = zeros(typeof(deltaArray[1][1]),ntraits)
    w        = zeros(typeof(wArray[1][1]),ntraits) #for rhs

    nlable    = length(state_labels)
    probDelta = Array{Float64}(undef, nlable)
    logDelta  = Array{Float64}(undef, nlable)
    βeta      = Array{Array{Float64,1}}(undef, nlable)
    RinvLhs   = Array{Array{Float64,2}}(undef, nlable) # D*inv(R)*D
    RinvRhs   = Array{Array{Float64,2}}(undef, nlable) # inv(R)*D

    for (iloci, state) in enumerate(state_labels)
        D  = diagm(state)
        RinvLhs[iloci] = D*Rinv*D # D*inv(R)*D
        RinvRhs[iloci] = Rinv*D   # inv(R)*D
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
        for (iloci, state) in enumerate(state_labels)
            lhs       = RinvLhs[iloci]*xpRinvx[marker]+Ginv[marker]  #C: t-by-t
            rhs       = RinvRhs[iloci]'w  #t-by-1
            invLhs    = inv(lhs)
            invLhsC   = cholesky(Hermitian(invLhs)).L # L, where LL'=invLhs
            gHat      = invLhs*rhs  #t-by-1
            logDelta[iloci] = -0.5*(log(det(lhs))-rhs'gHat)+log_prior(prior, marker, state)
            βeta[iloci]     = gHat + invLhsC*stdnorm  #var(Lz)=LL'=invLhs=inv(C)
        end

        #marginal full conditional probability of δ
        max_logDelta = maximum(logDelta)
        isfinite(max_logDelta) || error("All MTBayesABC sampler II state probabilities are zero or invalid.")
        denominator = 0.0
        for l in 1:nlable
            probDelta[l] = exp(logDelta[l] - max_logDelta)
            denominator += probDelta[l]
        end
        probDelta ./= denominator

        #choose label
        whichlabel = rand(Categorical(probDelta))
        δ           = copy(state_labels[whichlabel])
        β           = βeta[whichlabel]
        newα        = diagm(δ)*β

        for trait in 1:ntraits
            BLAS.axpy!(oldα[trait]-newα[trait],x,wArray[trait])
            betaArray[trait][marker]       = β[trait]
            deltaArray[trait][marker]      = δ[trait]
            alphaArray[trait][marker]      = newα[trait]
        end
    end
end

#block
function MTBayesABC_block!(genotypes,ycorr_array,vare,locus_effect_variances,Rinv=ones(eltype(ycorr_array[1]),length(ycorr_array[1])),independent_blocks=false)
    prior, sampler_mode = mt_bayesc_block_context(genotypes, genotypes.ntraits)
    if sampler_mode == :I
        if independent_blocks == true
            return _MTBayesABC_block_samplerI_independent!(genotypes.MArray,genotypes.mpRinvm,
                                                           genotypes.genotypes,genotypes.MpRinvM,
                                                           ycorr_array,genotypes.β,genotypes.δ,genotypes.α,vare,
                                                           locus_effect_variances,prior,Rinv)
        end
        return _MTBayesABC_block_samplerI!(genotypes.MArray,genotypes.mpRinvm,
                                           genotypes.genotypes,genotypes.MpRinvM,
                                           ycorr_array,genotypes.β,genotypes.δ,genotypes.α,vare,
                                           locus_effect_variances,prior,Rinv)
    elseif sampler_mode == :II
        state_labels = collect(keys(genotypes.π))
        if independent_blocks == true
            return _MTBayesABC_block_samplerII_independent!(genotypes.MArray,genotypes.mpRinvm,
                                                            genotypes.genotypes,genotypes.MpRinvM,
                                                            ycorr_array,genotypes.β,genotypes.δ,genotypes.α,vare,
                                                            locus_effect_variances,state_labels,prior,Rinv)
        end
        return _MTBayesABC_block_samplerII!(genotypes.MArray,genotypes.mpRinvm,
                                            genotypes.genotypes,genotypes.MpRinvM,
                                            ycorr_array,genotypes.β,genotypes.δ,genotypes.α,vare,
                                            locus_effect_variances,state_labels,prior,Rinv)
    else
        error("multi_trait_sampler must resolve to :I or :II.")
    end
end

function _MTBayesABC_block_samplerI!(XArray,xpRinvx,
                                     X, XpRinvX,
                                     wArray,
                                     betaArray,deltaArray,alphaArray,
                                     vare,varEffects,
                                     prior,Rinvw)
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
    unit_weights  = is_unit_weights(Rinvw)
    max_block_size = maximum(size(XpRinvX[i],1) for i in 1:nblocks)
    rhs_cache = [zeros(eltype(wArray[trait]),max_block_size) for trait = 1:ntraits]
    for i in 1:nblocks
        block_size  = size(XpRinvX[i],1)
        XpRinvycorr = [view(rhs_cache[trait],1:block_size) for trait = 1:ntraits]
        for trait = 1:ntraits
            block_rhs!(XpRinvycorr[trait],XArray[i],wArray[trait],Rinvw,unit_weights)
        end
        oldalphaArray = deepcopy(alphaArray)
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

                    logDelta0  = -0.5*(log(Ginv11)- gHat0^2*Ginv11) + log_prior(prior, locus_j, d0) #logPi
                    logDelta1  = -0.5*(log(C11)-gHat1^2*C11) + log_prior(prior, locus_j, d1) #logPiComp

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
                    betaArray[trait][locus_j]  = β[trait]
                    deltaArray[trait][locus_j] = δ[trait]
                    alphaArray[trait][locus_j] = newα[trait]
                end
            end
        end
        for trait = 1:ntraits
            wArray[trait][:] = wArray[trait] + XArray[i] * (oldalphaArray[trait] - alphaArray[trait])[(start_pos+1):(start_pos+block_size)]
        end
        start_pos += block_size
    end
end

function _MTBayesABC_block_samplerI_independent!(XArray,xpRinvx,
                                                 X, XpRinvX,
                                                 wArray,
                                                 betaArray,deltaArray,alphaArray,
                                                 vare,varEffects,
                                                 prior,Rinvw)
    ntraits  = length(alphaArray)

    Rinv     = inv(vare)
    Ginv     = inv.(varEffects)

    nblocks       = length(XpRinvX)
    block_sizes   = [size(XpRinvX[i],1) for i in 1:nblocks]
    block_offsets = cumsum([0; block_sizes[1:end-1]])
    unit_weights  = is_unit_weights(Rinvw)
    w_snapshot    = [copy(wArray[trait]) for trait = 1:ntraits]
    block_deltas  = Vector{Vector{Vector{eltype(alphaArray[1])}}}(undef,nblocks)
    block_seeds   = rand(UInt,nblocks)

    Threads.@threads for i in 1:nblocks
        rng = MersenneTwister(block_seeds[i])
        block_size  = block_sizes[i]
        start_pos   = block_offsets[i]
        block_range = start_pos + 1:start_pos + block_size
        XpRinvycorr = [zeros(eltype(wArray[trait]),block_size) for trait = 1:ntraits]
        for trait = 1:ntraits
            block_rhs!(XpRinvycorr[trait],XArray[i],w_snapshot[trait],Rinvw,unit_weights)
        end
        block_old_alpha = [copy(view(alphaArray[trait], block_range)) for trait in 1:ntraits]

        β        = zeros(typeof(betaArray[1][1]),ntraits)
        newα     = zeros(typeof(alphaArray[1][1]),ntraits)
        oldα     = zeros(typeof(alphaArray[1][1]),ntraits)
        δ        = zeros(typeof(deltaArray[1][1]),ntraits)
        w        = zeros(typeof(wArray[1][1]),ntraits)

        nreps = block_size
        for reps = 1:nreps
            for j = 1:block_size
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

                    logDelta0  = -0.5*(log(Ginv11)- gHat0^2*Ginv11) + log_prior(prior, locus_j, d0)
                    logDelta1  = -0.5*(log(C11)-gHat1^2*C11) + log_prior(prior, locus_j, d1)

                    probDelta1 =  1.0/(1.0+exp(logDelta0-logDelta1))
                    if(rand(rng)<probDelta1)
                        δ[k] = 1
                        β[k] = newα[k] = gHat1 + randn(rng)*sqrt(invLhs1)
                        BLAS.axpy!(oldα[k]-newα[k],view(XpRinvX[i],:,j),XpRinvycorr[k])
                    else
                        β[k] = gHat0 + randn(rng)*sqrt(invLhs0)
                        δ[k] = 0
                        newα[k] = 0
                        if oldα[k] != 0
                            BLAS.axpy!(oldα[k],view(XpRinvX[i],:,j),XpRinvycorr[k])
                        end
                    end
                end
                for trait = 1:ntraits
                    betaArray[trait][locus_j]  = β[trait]
                    deltaArray[trait][locus_j] = δ[trait]
                    alphaArray[trait][locus_j] = newα[trait]
                end
            end
        end

        for trait = 1:ntraits
            block_old_alpha[trait] .-= alphaArray[trait][block_range]
        end
        block_deltas[i] = block_old_alpha
    end

    for i in 1:nblocks
        for trait = 1:ntraits
            mul!(wArray[trait], XArray[i], block_deltas[i][trait], 1.0, 1.0)
        end
    end
    return nothing
end

function _MTBayesABC_block_samplerII!(XArray,xpRinvx,
                                      X, XpRinvX,
                                      wArray,
                                      betaArray,deltaArray,alphaArray,
                                      vare,varEffects,
                                      state_labels,
                                      prior,Rinvw)
    nMarkers = length(xpRinvx)
    ntraits  = length(alphaArray)

    Rinv     = inv(vare) #Do Not Use inv.(): elementwise inversion
    Ginv     = inv.(varEffects)

    β        = zeros(typeof(betaArray[1][1]),ntraits)
    newα     = zeros(typeof(alphaArray[1][1]),ntraits)
    oldα     = zeros(typeof(alphaArray[1][1]),ntraits)
    δ        = zeros(typeof(deltaArray[1][1]),ntraits)
    w        = zeros(typeof(wArray[1][1]),ntraits)

    nlable    = length(state_labels)
    probDelta = Array{Float64}(undef, nlable)
    logDelta  = Array{Float64}(undef, nlable)
    βeta      = Array{Array{Float64,1}}(undef, nlable)
    RinvLhs   = Array{Array{Float64,2}}(undef, nlable)
    RinvRhs   = Array{Array{Float64,2}}(undef, nlable)

    for (iloci, state) in enumerate(state_labels)
        D  = diagm(state)
        RinvLhs[iloci] = D*Rinv*D
        RinvRhs[iloci] = Rinv*D
    end

    nblocks       = length(XpRinvX)
    start_pos     = 0
    unit_weights  = is_unit_weights(Rinvw)
    max_block_size = maximum(size(XpRinvX[i],1) for i in 1:nblocks)
    rhs_cache = [zeros(eltype(wArray[trait]),max_block_size) for trait = 1:ntraits]

    for i in 1:nblocks
        block_size  = size(XpRinvX[i],1)
        block_range = start_pos + 1:start_pos + block_size
        XpRinvycorr = [view(rhs_cache[trait],1:block_size) for trait = 1:ntraits]
        for trait = 1:ntraits
            block_rhs!(XpRinvycorr[trait],XArray[i],wArray[trait],Rinvw,unit_weights)
        end

        block_old_alpha = [copy(view(alphaArray[trait], block_range)) for trait in 1:ntraits]
        nreps = block_size

        for reps = 1:nreps
            for j = 1:block_size
                locus_j = start_pos + j
                for trait = 1:ntraits
                    β[trait]  = betaArray[trait][locus_j]
                    oldα[trait]  = newα[trait] = alphaArray[trait][locus_j]
                    δ[trait]  = deltaArray[trait][locus_j]
                    w[trait]  = XpRinvycorr[trait][j] + xpRinvx[locus_j]*oldα[trait]
                end

                stdnorm = randn(ntraits)
                for (iloci, state) in enumerate(state_labels)
                    lhs       = RinvLhs[iloci]*xpRinvx[locus_j] + Ginv[locus_j]
                    rhs       = RinvRhs[iloci]'*w
                    invLhs    = inv(lhs)
                    invLhsC   = cholesky(Hermitian(invLhs)).L
                    gHat      = invLhs*rhs
                    logDelta[iloci] = -0.5*(log(det(lhs))-rhs'gHat) + log_prior(prior, locus_j, state)
                    βeta[iloci]     = gHat + invLhsC*stdnorm
                end

                max_logDelta = maximum(logDelta)
                isfinite(max_logDelta) || error("All MTBayesABC block sampler II state probabilities are zero or invalid.")
                denominator = 0.0
                for l in 1:nlable
                    probDelta[l] = exp(logDelta[l] - max_logDelta)
                    denominator += probDelta[l]
                end
                probDelta ./= denominator

                whichlabel = rand(Categorical(probDelta))
                δ = copy(state_labels[whichlabel])
                β = βeta[whichlabel]
                newα = diagm(δ)*β

                for trait = 1:ntraits
                    BLAS.axpy!(oldα[trait]-newα[trait], view(XpRinvX[i], :, j), XpRinvycorr[trait])
                    betaArray[trait][locus_j]  = β[trait]
                    deltaArray[trait][locus_j] = δ[trait]
                    alphaArray[trait][locus_j] = newα[trait]
                end
            end
        end

        for trait = 1:ntraits
            wArray[trait][:] = wArray[trait] + XArray[i] * (block_old_alpha[trait] - alphaArray[trait][block_range])
        end
        start_pos += block_size
    end
end

function _MTBayesABC_block_samplerII_independent!(XArray,xpRinvx,
                                                  X, XpRinvX,
                                                  wArray,
                                                  betaArray,deltaArray,alphaArray,
                                                  vare,varEffects,
                                                  state_labels,
                                                  prior,Rinvw)
    ntraits  = length(alphaArray)

    Rinv     = inv(vare)
    Ginv     = inv.(varEffects)

    nlable    = length(state_labels)
    RinvLhs   = Array{Array{Float64,2}}(undef, nlable)
    RinvRhs   = Array{Array{Float64,2}}(undef, nlable)

    for (iloci, state) in enumerate(state_labels)
        D  = diagm(state)
        RinvLhs[iloci] = D*Rinv*D
        RinvRhs[iloci] = Rinv*D
    end

    nblocks       = length(XpRinvX)
    block_sizes   = [size(XpRinvX[i],1) for i in 1:nblocks]
    block_offsets = cumsum([0; block_sizes[1:end-1]])
    unit_weights  = is_unit_weights(Rinvw)
    w_snapshot    = [copy(wArray[trait]) for trait = 1:ntraits]
    block_deltas  = Vector{Vector{Vector{eltype(alphaArray[1])}}}(undef,nblocks)
    block_seeds   = rand(UInt,nblocks)

    Threads.@threads for i in 1:nblocks
        rng = MersenneTwister(block_seeds[i])
        block_size  = block_sizes[i]
        start_pos   = block_offsets[i]
        block_range = start_pos + 1:start_pos + block_size
        XpRinvycorr = [zeros(eltype(wArray[trait]),block_size) for trait = 1:ntraits]
        for trait = 1:ntraits
            block_rhs!(XpRinvycorr[trait],XArray[i],w_snapshot[trait],Rinvw,unit_weights)
        end
        block_old_alpha = [copy(view(alphaArray[trait], block_range)) for trait in 1:ntraits]

        β        = zeros(typeof(betaArray[1][1]),ntraits)
        newα     = zeros(typeof(alphaArray[1][1]),ntraits)
        oldα     = zeros(typeof(alphaArray[1][1]),ntraits)
        δ        = zeros(typeof(deltaArray[1][1]),ntraits)
        w        = zeros(typeof(wArray[1][1]),ntraits)
        probDelta = Array{Float64}(undef, nlable)
        logDelta  = Array{Float64}(undef, nlable)
        βeta      = Array{Array{Float64,1}}(undef, nlable)

        nreps = block_size
        for reps = 1:nreps
            for j = 1:block_size
                locus_j = start_pos + j
                for trait = 1:ntraits
                    β[trait]  = betaArray[trait][locus_j]
                    oldα[trait]  = newα[trait] = alphaArray[trait][locus_j]
                    δ[trait]  = deltaArray[trait][locus_j]
                    w[trait]  = XpRinvycorr[trait][j] + xpRinvx[locus_j]*oldα[trait]
                end

                stdnorm = randn(rng,ntraits)
                for (iloci, state) in enumerate(state_labels)
                    lhs       = RinvLhs[iloci]*xpRinvx[locus_j] + Ginv[locus_j]
                    rhs       = RinvRhs[iloci]'*w
                    invLhs    = inv(lhs)
                    invLhsC   = cholesky(Hermitian(invLhs)).L
                    gHat      = invLhs*rhs
                    logDelta[iloci] = -0.5*(log(det(lhs))-rhs'gHat) + log_prior(prior, locus_j, state)
                    βeta[iloci]     = gHat + invLhsC*stdnorm
                end

                max_logDelta = maximum(logDelta)
                isfinite(max_logDelta) || error("All MTBayesABC independent block sampler II state probabilities are zero or invalid.")
                denominator = 0.0
                for l in 1:nlable
                    probDelta[l] = exp(logDelta[l] - max_logDelta)
                    denominator += probDelta[l]
                end
                probDelta ./= denominator

                whichlabel = rand(rng,Categorical(probDelta))
                δ = copy(state_labels[whichlabel])
                β = βeta[whichlabel]
                newα = diagm(δ)*β

                for trait = 1:ntraits
                    BLAS.axpy!(oldα[trait]-newα[trait], view(XpRinvX[i], :, j), XpRinvycorr[trait])
                    betaArray[trait][locus_j]  = β[trait]
                    deltaArray[trait][locus_j] = δ[trait]
                    alphaArray[trait][locus_j] = newα[trait]
                end
            end
        end

        for trait = 1:ntraits
            block_old_alpha[trait] .-= alphaArray[trait][block_range]
        end
        block_deltas[i] = block_old_alpha
    end

    for i in 1:nblocks
        for trait = 1:ntraits
            mul!(wArray[trait], XArray[i], block_deltas[i][trait], 1.0, 1.0)
        end
    end
    return nothing
end
