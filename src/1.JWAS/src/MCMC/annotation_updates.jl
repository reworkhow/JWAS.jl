function has_marker_annotations(Mi)
    return Mi.annotations !== false
end

function initialize_annotation_indicators!(Mi)
    if !has_marker_annotations(Mi) || Mi.nMarkers <= 1
        return nothing
    end
    ann = Mi.annotations
    if ann.nsteps != 1
        delta = Mi.δ[1]
        priors = ann.snp_pi === false ? repeat(reshape(Float64.(Mi.π), 1, :), Mi.nMarkers, 1) : ann.snp_pi
        for j in eachindex(delta)
            delta[j] = rand(Categorical(Vector(view(priors, j, :))))
        end
        return nothing
    end
    exclude_prob = Mi.π isa AbstractVector ? Float64(Mi.π[1]) : Float64(Mi.π)
    exclude_prob = clamp(exclude_prob, 0.0, 1.0)
    delta = Mi.δ[1]
    included = one(eltype(delta))
    excluded = zero(eltype(delta))
    minority_count = clamp(round(Int, 0.1 * Mi.nMarkers), 1, Mi.nMarkers - 1)

    if exclude_prob == 0.0
        @info "Annotated BayesC initialization: starting pi=0.0 is degenerate; using 10% excluded markers."
        fill!(delta, included)
        delta[sample(1:Mi.nMarkers, minority_count; replace=false)] .= excluded
        return nothing
    elseif exclude_prob == 1.0
        @info "Annotated BayesC initialization: starting pi=1.0 is degenerate; using 10% included markers."
        fill!(delta, excluded)
        delta[sample(1:Mi.nMarkers, minority_count; replace=false)] .= included
        return nothing
    end

    for j in eachindex(delta)
        delta[j] = rand() < exclude_prob ? excluded : included
    end

    if all(iszero, delta)
        @info "Annotated BayesC initialization: sampled all markers as excluded; reinitializing with 10% included markers."
        delta[sample(1:Mi.nMarkers, minority_count; replace=false)] .= included
    elseif all(isone, delta)
        @info "Annotated BayesC initialization: sampled all markers as included; reinitializing with 10% excluded markers."
        delta[sample(1:Mi.nMarkers, minority_count; replace=false)] .= excluded
    end
    return nothing
end

function update_annotation_bounds_single_step!(Mi)
    ann = Mi.annotations
    ann.lower_bound .= ann.thresholds[Int.(Mi.δ[1]) .+ 1]
    ann.upper_bound .= ann.thresholds[Int.(Mi.δ[1]) .+ 2]
    return nothing
end

function sample_annotation_liabilities_single_step!(Mi)
    ann = Mi.annotations
    update_annotation_bounds_single_step!(Mi)
    for i in 1:Mi.nMarkers
        lower = ann.lower_bound[i]
        upper = ann.upper_bound[i]
        if lower == upper
            ann.liability[i] = lower
        else
            ann.liability[i] = rand(truncated(Normal(ann.mu[i], sqrt(ann.variance)), lower, upper))
        end
    end
    return nothing
end

function update_annotation_priors_single_step!(Mi)
    ann = Mi.annotations
    sample_annotation_liabilities_single_step!(Mi)
    rhs = ann.design_matrix' * ann.liability
    Gibbs(ann.lhs, ann.coefficients, rhs, ann.variance)
    ann.mu .= ann.design_matrix * ann.coefficients
    dist = Normal(0, sqrt(ann.variance))
    Mi.π .= clamp.(1 .- cdf.(dist, ann.mu), eps(Float64), 1 - eps(Float64))
    return nothing
end

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

function sample_bayesr_annotation_step!(ann, step::Integer, response::AbstractVector{<:Integer}, active::AbstractVector{<:Integer})
    coeffs = view(ann.coefficients, :, step)
    ann.mu[:, step] .= ann.design_matrix * coeffs
    ann.lower_bound[:, step] .= -Inf
    ann.upper_bound[:, step] .= Inf
    isempty(active) && return nothing

    X = ann.design_matrix[active, :]
    mu_active = ann.mu[active, step]
    variance = ann.variance[step]
    liability_active = view(ann.liability, active, step)

    for (local_i, global_i) in pairs(active)
        if response[global_i] == 0
            ann.lower_bound[global_i, step] = -Inf
            ann.upper_bound[global_i, step] = 0.0
        else
            ann.lower_bound[global_i, step] = 0.0
            ann.upper_bound[global_i, step] = Inf
        end
        liability_active[local_i] = rand(truncated(Normal(mu_active[local_i], 1.0), ann.lower_bound[global_i, step], ann.upper_bound[global_i, step]))
    end

    lj = liability_active .- mu_active

    old_sample = coeffs[1]
    rhs = sum(lj) + length(active) * old_sample
    inv_lhs = 1.0 / length(active)
    ahat = inv_lhs * rhs
    coeffs[1] = randn() * sqrt(inv_lhs) + ahat
    lj .+= old_sample - coeffs[1]

    if size(X, 2) > 1
        for k in 2:size(X, 2)
            old_sample = coeffs[k]
            anno_diag = dot(view(X, :, k), view(X, :, k))
            inv_lhs = 1.0 / (anno_diag + 1.0 / variance)
            ahat = inv_lhs * (dot(view(X, :, k), lj) + anno_diag * old_sample)
            coeffs[k] = randn() * sqrt(inv_lhs) + ahat
            lj .+= view(X, :, k) .* (old_sample - coeffs[k])
        end
        ann.variance[step] = (sum(abs2, view(coeffs, 2:length(coeffs))) + 2.0) / rand(Chisq(size(X, 2) + 1.0))
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
