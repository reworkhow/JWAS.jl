function has_marker_annotations(Mi)
    return Mi.annotations !== false
end

"""
    annotation_binary_bounds!(lower, upper, response)

Set the truncation bounds for a binary probit step.

For binary indicator `z_i ∈ {0, 1}`, introduce latent liability `l_i` such that

`z_i = 1(l_i > 0)`.

Then the full conditional is

`l_i | z_i, μ_i ~ N(μ_i, s^2)` truncated to:
- `(-Inf, 0]` when `z_i = 0`
- `[0, Inf)` when `z_i = 1`
"""
function annotation_binary_bounds!(lower::AbstractVector, upper::AbstractVector, response::AbstractVector{<:Real})
    for i in eachindex(response, lower, upper)
        if response[i] == 0
            lower[i] = -Inf
            upper[i] = 0.0
        else
            lower[i] = 0.0
            upper[i] = Inf
        end
    end
    return nothing
end

"""
    sample_binary_annotation_liabilities!(liability, mu, lower, upper, response; latent_sd)

Sample the latent liabilities for one binary probit step after the truncation
bounds have been determined by [`annotation_binary_bounds!`](@ref).

This is the common latent-variable update used by:
- annotated BayesC for the single inclusion step
- annotated BayesR for each of the nested step-up indicators
"""
function sample_binary_annotation_liabilities!(liability::AbstractVector,
                                              mu::AbstractVector,
                                              lower::AbstractVector,
                                              upper::AbstractVector,
                                              response::AbstractVector{<:Real};
                                              latent_sd::Real)
    annotation_binary_bounds!(lower, upper, response)
    for i in eachindex(response, liability, mu, lower, upper)
        if lower[i] == upper[i]
            liability[i] = lower[i]
        else
            liability[i] = rand(truncated(Normal(mu[i], latent_sd), lower[i], upper[i]))
        end
    end
    return nothing
end

"""
    gibbs_update_annotation_coefficients!(ann)

BayesC one-step coefficient update.

This keeps the existing JWAS BayesC behavior: a single multivariate Gibbs update
on the latent regression

`l = Xα + ε`.
"""
function gibbs_update_annotation_coefficients!(ann)
    rhs = ann.design_matrix' * ann.liability
    Gibbs(ann.lhs, ann.coefficients, rhs, ann.variance)
    ann.mu .= ann.design_matrix * ann.coefficients
    return nothing
end

"""
    gibbs_update_annotation_coefficients!(coeffs, X, latent_residual, variance)

BayesR step-wise coefficient update for one binary probit submodel.

Write the latent regression as

`l = Xα + ε`, with `ε ~ N(0, I)`.

After sampling the latent liabilities, we work with the residual

`r = l - Xα`.

Then each coefficient is updated by a standard scalar Gibbs step:
- the intercept uses a flat prior
- slopes use `α_k ~ N(0, σ^2_α)`
"""
function gibbs_update_annotation_coefficients!(coeffs::AbstractVector,
                                               X::AbstractMatrix,
                                               latent_residual::AbstractVector,
                                               variance::Real)
    nobs = size(X, 1)

    old_sample = coeffs[1]
    rhs = sum(latent_residual) + nobs * old_sample
    inv_lhs = 1.0 / nobs
    ahat = inv_lhs * rhs
    coeffs[1] = randn() * sqrt(inv_lhs) + ahat
    latent_residual .+= old_sample - coeffs[1]

    if size(X, 2) > 1
        for k in 2:size(X, 2)
            old_sample = coeffs[k]
            xk = view(X, :, k)
            anno_diag = dot(xk, xk)
            inv_lhs = 1.0 / (anno_diag + 1.0 / variance)
            ahat = inv_lhs * (dot(xk, latent_residual) + anno_diag * old_sample)
            coeffs[k] = randn() * sqrt(inv_lhs) + ahat
            latent_residual .+= xk .* (old_sample - coeffs[k])
        end
    end
    return nothing
end

"""
    sample_annotation_effect_variance!(variance, coeffs)

Update the slope variance for one annotation step using the same scaled
inverse-chi-square form as Jian's `sbayesrc.R`:

`σ^2_α = (Σ_{k>1} α_k^2 + 2) / χ^2_{p+1}`

where `p` is the total number of annotation coefficients including the intercept.
"""
function sample_annotation_effect_variance!(variance::AbstractVector, step::Integer, coeffs::AbstractVector)
    variance[step] = (sum(abs2, view(coeffs, 2:length(coeffs))) + 2.0) / rand(Chisq(length(coeffs) + 1.0))
    return nothing
end

"""
    update_annotation_bounds_single_step!(Mi)

Annotated BayesC is the one-step binary special case. The thresholds remain the
current BayesC convention, and this refactor preserves that behavior.
"""
function update_annotation_bounds_single_step!(Mi)
    ann = Mi.annotations
    ann.lower_bound .= ann.thresholds[Int.(Mi.δ[1]) .+ 1]
    ann.upper_bound .= ann.thresholds[Int.(Mi.δ[1]) .+ 2]
    return nothing
end

function sample_annotation_liabilities_single_step!(Mi)
    ann = Mi.annotations
    # Annotated BayesC uses the same standard-probit convention as annotated
    # BayesR: the latent annotation error variance is fixed to 1 for
    # identifiability, while ann.variance is reserved for coefficient shrinkage.
    sample_binary_annotation_liabilities!(
        ann.liability,
        ann.mu,
        ann.lower_bound,
        ann.upper_bound,
        Mi.δ[1];
        latent_sd=1.0,
    )
    return nothing
end

function update_annotation_priors_single_step!(Mi)
    ann = Mi.annotations
    update_annotation_bounds_single_step!(Mi)
    sample_annotation_liabilities_single_step!(Mi)
    gibbs_update_annotation_coefficients!(ann)
    pi_values = clamp.(1 .- cdf.(Normal(), ann.mu), eps(Float64), 1 - eps(Float64))
    if Mi.π isa AbstractVector && length(Mi.π) == length(pi_values)
        Mi.π .= pi_values
    else
        Mi.π = copy(pi_values)
    end
    return nothing
end

"""
    bayesr_annotation_step_indicators(delta)

Build the three nested BayesR step-up indicators:

- `z1_j = 1(δ_j > 1)`
- `z2_j = 1(δ_j > 2)`
- `z3_j = 1(δ_j > 3)`

and the corresponding active subsets used by the conditional probit updates.
"""
function bayesr_annotation_step_indicators(delta::AbstractVector{<:Integer})
    z1 = Int.(delta .> 1)
    z2 = Int.(delta .> 2)
    z3 = Int.(delta .> 3)
    active = (
        collect(eachindex(delta)),
        findall(!iszero, z1),
        findall(!iszero, z2),
    )
    return (z1, z2, z3), active
end

"""
    sample_bayesr_annotation_step!(ann, step, response, active)

Sample one conditional BayesR annotation step.

For step-specific conditional probabilities `(p1, p2, p3)`, BayesR reconstructs
the 4-class per-marker prior as

- `π_{j1} = 1 - p1_j`
- `π_{j2} = p1_j (1 - p2_j)`
- `π_{j3} = p1_j p2_j (1 - p3_j)`
- `π_{j4} = p1_j p2_j p3_j`

This helper updates one binary probit submodel that contributes to `p_step`.
"""
function sample_bayesr_annotation_step!(ann, step::Integer, response::AbstractVector{<:Integer}, active::AbstractVector{<:Integer})
    coeffs = view(ann.coefficients, :, step)
    ann.mu[:, step] .= ann.design_matrix * coeffs
    ann.lower_bound[:, step] .= -Inf
    ann.upper_bound[:, step] .= Inf
    isempty(active) && return nothing

    X = ann.design_matrix[active, :]
    mu_active = view(ann.mu, active, step)
    liability_active = view(ann.liability, active, step)
    lower_active = view(ann.lower_bound, active, step)
    upper_active = view(ann.upper_bound, active, step)
    response_active = view(response, active)

    sample_binary_annotation_liabilities!(
        liability_active,
        mu_active,
        lower_active,
        upper_active,
        response_active;
        latent_sd=1.0,
    )

    latent_residual = liability_active .- mu_active
    gibbs_update_annotation_coefficients!(coeffs, X, latent_residual, ann.variance[step])

    if size(X, 2) > 1
        sample_annotation_effect_variance!(ann.variance, step, coeffs)
    end

    ann.mu[:, step] .= ann.design_matrix * coeffs
    return nothing
end

function rebuild_bayesr_annotation_priors!(ann)
    probs = clamp.(cdf.(Normal(), ann.mu), eps(Float64), 1 - eps(Float64))
    ann.snp_pi[:, 1] .= 1 .- probs[:, 1]
    ann.snp_pi[:, 2] .= probs[:, 1] .* (1 .- probs[:, 2])
    ann.snp_pi[:, 3] .= probs[:, 1] .* probs[:, 2] .* (1 .- probs[:, 3])
    ann.snp_pi[:, 4] .= probs[:, 1] .* probs[:, 2] .* probs[:, 3]
    return nothing
end

function update_annotation_priors!(Mi)
    if !has_marker_annotations(Mi)
        return nothing
    end
    ann = Mi.annotations
    if ann.nsteps == 1
        return update_annotation_priors_single_step!(Mi)
    end

    delta = Int.(Mi.δ[1])
    responses, active_sets = bayesr_annotation_step_indicators(delta)
    for step in 1:ann.nsteps
        sample_bayesr_annotation_step!(ann, step, responses[step], active_sets[step])
    end
    rebuild_bayesr_annotation_priors!(ann)
    mean_pi = vec(mean(ann.snp_pi, dims=1))
    if Mi.π isa AbstractVector
        Mi.π .= mean_pi
    else
        Mi.π = mean_pi
    end
    return nothing
end
