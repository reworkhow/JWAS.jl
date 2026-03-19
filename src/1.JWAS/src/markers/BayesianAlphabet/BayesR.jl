@inline function bayesr_logsumexp(log_probs)
    max_log = maximum(log_probs)
    return max_log + log(sum(exp.(log_probs .- max_log)))
end

function BayesR!(genotypes, ycorr, vare)
    BayesR!(genotypes.mArray, genotypes.mRinvArray, genotypes.mpRinvm,
            ycorr, genotypes.α[1], genotypes.δ[1], vare, genotypes.G.val, genotypes.π, BAYESR_GAMMA)
end

function BayesR!(xArray, xRinvArray, xpRinvx,
                 yCorr,
                 α, δ,
                 vare, sigmaSq, π, gamma)
    nclasses = length(gamma)
    length(π) == nclasses || error("BayesR pi vector length $(length(π)) must match the number of mixture classes ($nclasses).")
    isapprox(sum(π), 1.0; atol=1e-8) || error("BayesR pi must sum to 1.")
    sigmaSq > 0 || error("BayesR sigmaSq must be positive.")

    log_probs = zeros(Float64, nclasses)
    probs = zeros(Float64, nclasses)
    logPi = log.(π)
    invVarRes = 1 / vare

    for j in eachindex(α)
        x = xArray[j]
        xRinv = xRinvArray[j]
        rhs = (dot(xRinv, yCorr) + xpRinvx[j] * α[j]) * invVarRes
        oldAlpha = α[j]

        log_probs[1] = logPi[1]
        for k in 2:nclasses
            varEffect = gamma[k] * sigmaSq
            invVarEffect = 1 / varEffect
            lhs = xpRinvx[j] * invVarRes + invVarEffect
            invLhs = 1 / lhs
            betaHat = invLhs * rhs
            log_probs[k] = 0.5 * (log(invLhs) - log(varEffect) + betaHat * rhs) + logPi[k]
        end

        log_norm = bayesr_logsumexp(log_probs)
        for k in 1:nclasses
            probs[k] = exp(log_probs[k] - log_norm)
        end

        sampled_class = rand(Categorical(probs))
        δ[j] = sampled_class

        if sampled_class == 1
            if oldAlpha != 0
                BLAS.axpy!(oldAlpha, x, yCorr)
            end
            α[j] = zero(eltype(α))
        else
            varEffect = gamma[sampled_class] * sigmaSq
            invVarEffect = 1 / varEffect
            lhs = xpRinvx[j] * invVarRes + invVarEffect
            invLhs = 1 / lhs
            betaHat = invLhs * rhs
            α[j] = betaHat + randn() * sqrt(invLhs)
            BLAS.axpy!(oldAlpha - α[j], x, yCorr)
        end
    end
end
