@inline bayesr_mt_prior_row(π::AbstractVector, ::Integer) = π
@inline bayesr_mt_prior_row(π::AbstractMatrix, marker::Integer) = view(π, marker, :)

function bayesr_mt_validate_priors(π::AbstractVector, nmarkers::Integer)
    length(π) == length(ANNOTATED_BAYESR_MT_STATES) || error("Multi-trait BayesR pi vector must have $(length(ANNOTATED_BAYESR_MT_STATES)) entries.")
    all(x -> x >= 0, π) || error("Multi-trait BayesR pi entries must be nonnegative.")
    isapprox(sum(π), 1.0; atol=1e-8) || error("Multi-trait BayesR pi must sum to 1.")
    return nothing
end

function bayesr_mt_validate_priors(π::AbstractMatrix, nmarkers::Integer)
    size(π, 1) == nmarkers || error("Multi-trait BayesR per-marker pi must have one row per marker.")
    size(π, 2) == length(ANNOTATED_BAYESR_MT_STATES) || error("Multi-trait BayesR per-marker pi must have $(length(ANNOTATED_BAYESR_MT_STATES)) columns.")
    return nothing
end

function mt_bayesr_trait_class_conditionals(trait::Integer,
                                            marker::Integer,
                                            w::AbstractVector,
                                            beta::AbstractVector,
                                            delta::AbstractVector,
                                            xpRinvx::Real,
                                            Rinv::AbstractMatrix,
                                            Ginv::AbstractMatrix,
                                            prior_row::AbstractVector,
                                            gamma::AbstractVector)
    nclasses = length(gamma)
    log_weights = zeros(Float64, nclasses)
    means = zeros(Float64, nclasses)
    vars = zeros(Float64, nclasses)
    others = [idx for idx in eachindex(beta) if idx != trait]

    for class in 1:nclasses
        candidate_delta = Int.(delta)
        candidate_delta[trait] = class
        prior_prob = Float64(prior_row[annotated_bayesr_mt_state_index(candidate_delta)])
        scale = sqrt(Float64(gamma[class]))

        lhs = Float64(Ginv[trait, trait]) + Float64(xpRinvx) * scale^2 * Float64(Rinv[trait, trait])
        lhs > 0.0 || error("Multi-trait BayesR conditional precision must be positive.")
        cross = [
            Float64(Ginv[trait, other]) +
            Float64(xpRinvx) * scale * sqrt(Float64(gamma[Int(delta[other])])) * Float64(Rinv[trait, other])
            for other in others
        ]
        rhs = scale * dot(w, view(Rinv, :, trait)) - dot(cross, view(beta, others))
        inv_lhs = 1.0 / lhs

        means[class] = inv_lhs * rhs
        vars[class] = inv_lhs
        log_weights[class] = prior_prob <= 0.0 ? -Inf : -0.5 * (log(lhs) - means[class] * rhs) + log(prior_prob)
    end

    return log_weights, means, vars
end

function MTBayesR!(genotypes, ycorr_array, vare)
    priors = genotypes.annotations === false ? genotypes.π : genotypes.annotations.snp_pi
    MTBayesR!(genotypes.mArray, genotypes.mRinvArray, genotypes.mpRinvm,
              ycorr_array, genotypes.β, genotypes.δ, genotypes.α,
              vare, genotypes.G.val, priors, BAYESR_GAMMA)
end

function MTBayesR!(xArray,
                   xRinvArray,
                   xpRinvx,
                   wArray,
                   betaArray,
                   deltaArray,
                   alphaArray,
                   vare,
                   G,
                   priors,
                   gamma)
    nmarkers = length(xArray)
    ntraits = length(alphaArray)
    nclasses = length(gamma)
    ntraits == 2 || error("Multi-trait BayesR currently supports exactly 2 traits.")
    bayesr_mt_validate_priors(priors, nmarkers)

    Rinv = inv(vare)
    Ginv = inv(G)
    beta = zeros(Float64, ntraits)
    delta = ones(Int, ntraits)
    old_alpha = zeros(Float64, ntraits)
    new_alpha = zeros(Float64, ntraits)
    w = zeros(Float64, ntraits)
    probs = zeros(Float64, nclasses)

    for marker in 1:nmarkers
        x = xArray[marker]
        xRinv = xRinvArray[marker]
        marker_xpRinvx = xpRinvx[marker]

        for trait in 1:ntraits
            beta[trait] = betaArray[trait][marker]
            delta[trait] = Int(deltaArray[trait][marker])
            old_alpha[trait] = alphaArray[trait][marker]
            new_alpha[trait] = old_alpha[trait]
            w[trait] = dot(xRinv, wArray[trait]) + marker_xpRinvx * old_alpha[trait]
        end

        prior_row = bayesr_mt_prior_row(priors, marker)
        for trait in 1:ntraits
            log_weights, means, vars = mt_bayesr_trait_class_conditionals(
                trait, marker, w, beta, delta, marker_xpRinvx, Rinv, Ginv, prior_row, gamma,
            )
            isfinite(maximum(log_weights)) || error("All multi-trait BayesR class probabilities are zero or invalid.")
            log_norm = bayesr_logsumexp(log_weights)
            for class in 1:nclasses
                probs[class] = exp(log_weights[class] - log_norm)
            end

            sampled_class = rand(Categorical(probs))
            previous_alpha = new_alpha[trait]
            delta[trait] = sampled_class
            beta[trait] = means[sampled_class] + randn() * sqrt(vars[sampled_class])
            new_alpha[trait] = sqrt(Float64(gamma[sampled_class])) * beta[trait]
            BLAS.axpy!(previous_alpha - new_alpha[trait], x, wArray[trait])
        end

        for trait in 1:ntraits
            betaArray[trait][marker] = beta[trait]
            deltaArray[trait][marker] = delta[trait]
            alphaArray[trait][marker] = new_alpha[trait]
        end
    end
    return nothing
end
