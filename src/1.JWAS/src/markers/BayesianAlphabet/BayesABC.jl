function megaBayesABC!(genotypes,wArray,vare,locus_effect_variances)
    Threads.@threads for i in 1:length(wArray) #ntraits
         BayesABC!(genotypes.mArray,genotypes.mRinvArray,genotypes.mpRinvm,
                    wArray[i],genotypes.α[i],genotypes.β[i],genotypes.δ[i],vare[i,i],
                    [vari[i,i] for vari in locus_effect_variances],genotypes.π[i])
    end
end


function BayesABC!(genotypes,ycorr,vare,locus_effect_variances)
    BayesABC!(genotypes.mArray,genotypes.mRinvArray,genotypes.mpRinvm,
              ycorr,genotypes.α[1],genotypes.β[1],genotypes.δ[1],vare,
              locus_effect_variances,genotypes.π)
end

@inline function bayesabc_update_marker!(x,
                                         yCorr,
                                         α,β,δ,
                                         j::Int,
                                         xRinvy,
                                         xpRinvx_j,
                                         invVarRes,
                                         invVarEffect_j,
                                         logVarEffect_j,
                                         varEffect_j,
                                         logDelta0,
                                         logPiComp)
    rhs = (xRinvy + xpRinvx_j*α[j])*invVarRes
    lhs = xpRinvx_j*invVarRes + invVarEffect_j
    invLhs = 1/lhs
    gHat   = rhs*invLhs
    logDelta1  = -0.5*(log(lhs) + logVarEffect_j - gHat*rhs) + logPiComp
    probDelta1 = 1/(1+ exp(logDelta0 - logDelta1))
    oldAlpha = α[j]

    if rand() < probDelta1
        δ[j] = 1
        β[j] = gHat + randn()*sqrt(invLhs)
        α[j] = β[j]
        BLAS.axpy!(oldAlpha-α[j],x,yCorr)
    else
        if oldAlpha != 0
            BLAS.axpy!(oldAlpha,x,yCorr)
        end
        δ[j] = 0
        β[j] = randn()*sqrt(varEffect_j)
        α[j] = 0
    end
    return nothing
end

function BayesABC!(xArray,xRinvArray,xpRinvx,
                   yCorr,
                   α,β,δ,
                   vare,varEffects,π)

    logPi         = log(π)
    logPiComp     = log(1-π)
    logDelta0     = logPi
    invVarRes     = 1/vare
    invVarEffects = 1 ./  varEffects
    logVarEffects = log.(varEffects)
    nMarkers      = length(α)

    for j=1:nMarkers
        x, xRinv = xArray[j], xRinvArray[j]
        bayesabc_update_marker!(x,yCorr,α,β,δ,j,
                                dot(xRinv,yCorr),xpRinvx[j],
                                invVarRes,invVarEffects[j],logVarEffects[j],varEffects[j],
                                logDelta0,logPiComp)
    end
end

function BayesABC_streaming!(genotypes,ycorr,vare,locus_effect_variances)
    BayesABC_streaming!(genotypes.stream_backend,genotypes.mpRinvm,
                        ycorr,genotypes.α[1],genotypes.β[1],genotypes.δ[1],vare,
                        locus_effect_variances,genotypes.π)
end

function BayesABC_streaming!(backend,xpRinvx,
                             yCorr,
                             α,β,δ,
                             vare,varEffects,π)
    logPi         = log(π)
    logPiComp     = log(1-π)
    logDelta0     = logPi
    invVarRes     = 1/vare
    invVarEffects = 1 ./ varEffects
    logVarEffects = log.(varEffects)
    nMarkers      = length(α)
    marker_buffer = Vector{eltype(yCorr)}(undef, backend.nObs)

    for j in 1:nMarkers
        decode_marker!(marker_buffer,backend,j)
        bayesabc_update_marker!(marker_buffer,yCorr,α,β,δ,j,
                                dot(marker_buffer,yCorr),xpRinvx[j],
                                invVarRes,invVarEffects[j],logVarEffects[j],varEffects[j],
                                logDelta0,logPiComp)
    end
end


function BayesABC_block!(genotypes,ycorr,vare,locus_effect_variances,Rinv=ones(eltype(ycorr),length(ycorr)))
    BayesABC_block!(genotypes.MArray,genotypes.mpRinvm,
              genotypes.genotypes,genotypes.MpRinvM,
              ycorr,genotypes.α[1],genotypes.β[1],genotypes.δ[1],vare,
              locus_effect_variances,genotypes.π,Rinv)
end

function BayesABC_block!(XArray,xpRinvx,
                   X, XpRinvX,
                   yCorr,
                   α,β,δ,
                   vare,varEffects,π,Rinv)

    logPi         = log(π)
    logPiComp     = log(1-π)
    logDelta0     = logPi
    invVarRes     = 1/vare
    invVarEffects = 1 ./  varEffects
    logVarEffects = log.(varEffects)
    nMarkers      = length(α)
    XpRinvX       = XpRinvX
    nblocks       = length(XpRinvX)
    start_pos     = 0
    max_block_size = maximum(size(XpRinvX[i],1) for i in 1:nblocks)
    αold_block = Vector{eltype(α)}(undef,max_block_size)
    rhs_block = Vector{eltype(yCorr)}(undef,max_block_size)
    unit_weights = is_unit_weights(Rinv)
    for i in 1:nblocks
        block_size  = size(XpRinvX[i],1)
        block_start = start_pos + 1
        block_end   = start_pos + block_size
        block_range = block_start:block_end
        copyto!(view(αold_block,1:block_size),view(α,block_range))
        XpRinvycorr = view(rhs_block,1:block_size)
        block_rhs!(XpRinvycorr,XArray[i],yCorr,Rinv,unit_weights)
        nreps       = block_size #user-defined nreps=block_size, outer_niter =niter/block_size
        for reps = 1:nreps
            for j=1:block_size #additional code to save all sampled αs is needed
                locus_j    = start_pos+j
                rhs        = (XpRinvycorr[j] + xpRinvx[locus_j]*α[locus_j])*invVarRes
                lhs        = xpRinvx[locus_j]*invVarRes + invVarEffects[locus_j]
                invLhs     = 1/lhs
                gHat       = rhs*invLhs
                logDelta1  = -0.5*(log(lhs) + logVarEffects[locus_j] - gHat*rhs) + logPiComp
                probDelta1 = 1/(1+ exp(logDelta0 - logDelta1))
                oldAlpha   = α[locus_j]

                if rand() < probDelta1
                    δ[locus_j] = 1
                    β[locus_j] = gHat + randn()*sqrt(invLhs)
                    α[locus_j] = β[locus_j]
                    BLAS.axpy!(oldAlpha-α[locus_j],view(XpRinvX[i],:,j),XpRinvycorr)
                else
                    if oldAlpha != 0
                        BLAS.axpy!(oldAlpha,view(XpRinvX[i],:,j),XpRinvycorr)
                    end
                    δ[locus_j] = 0
                    β[locus_j] = randn()*sqrt(varEffects[locus_j])
                    α[locus_j] = 0
                end
            end
        end

        @views begin
            αold_view = view(αold_block,1:block_size)
            αold_view .-= α[block_range]
            mul!(yCorr,XArray[i],αold_view,1.0,1.0)
        end
        start_pos += block_size
    end
end
