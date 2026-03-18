using CSV
using DataFrames
using Random
using Statistics

function normalize_pi_class_labels(labels)
    normalized = String[]
    for label in labels
        label_str = string(label)
        parsed = tryparse(Int, label_str)
        if parsed === nothing
            parsed_float = tryparse(Float64, label_str)
            if parsed_float !== nothing && isfinite(parsed_float) && isinteger(parsed_float)
                push!(normalized, "class$(Int(round(parsed_float)))")
            elseif startswith(label_str, "class")
                push!(normalized, label_str)
            else
                push!(normalized, label_str)
            end
        else
            push!(normalized, "class$(parsed)")
        end
    end
    return normalized
end

function build_bayesr_parity_dataset(; seed=2026, n_obs=40, n_markers=12)
    Random.seed!(seed)

    ids = ["id_$(i)" for i in 1:n_obs]
    marker_ids = ["m$(j)" for j in 1:n_markers]

    X = Matrix{Float64}(undef, n_obs, n_markers)
    allele_freq = rand(n_markers) .* 0.3 .+ 0.1
    for j in 1:n_markers
        p = allele_freq[j]
        for i in 1:n_obs
            X[i, j] = (rand() < p) + (rand() < p)
        end
    end

    beta_true = zeros(Float64, n_markers)
    for (idx, effect) in zip(1:min(3, n_markers), (0.8, -0.5, 0.3))
        beta_true[idx] = effect
    end
    mu_true = 1.0
    y = mu_true .+ X * beta_true .+ randn(n_obs)

    return (
        ids=ids,
        marker_ids=marker_ids,
        X=X,
        y=y,
        allele_freq=allele_freq,
        beta_true=beta_true,
        mu_true=mu_true,
    )
end

function build_bayesr_parity_initial_state(y, n_markers)
    mu0 = mean(y)
    return (
        beta0=zeros(Float64, n_markers),
        delta0=ones(Int, n_markers),
        mu0=mu0,
        ycorr0=Float64.(y .- mu0),
    )
end

function write_parity_dataset(outdir;
                              ids,
                              marker_ids,
                              X,
                              y,
                              gamma,
                              start_pi,
                              estimate_pi,
                              chain_length,
                              burnin,
                              start_h2,
                              start_sigma_sq,
                              start_vare,
                              seed,
                              initial_state=build_bayesr_parity_initial_state(y, length(marker_ids)))
    mkpath(outdir)

    geno_df = DataFrame(ID=ids)
    for (j, marker_id) in enumerate(marker_ids)
        geno_df[!, marker_id] = X[:, j]
    end
    CSV.write(joinpath(outdir, "genotypes.csv"), geno_df)

    pheno_df = DataFrame(ID=ids, y1=y)
    CSV.write(joinpath(outdir, "phenotypes.csv"), pheno_df)

    config_df = DataFrame(
        key=[
            "seed",
            "chain_length",
            "burnin",
            "estimate_pi",
            "start_h2",
            "start_sigma_sq",
            "start_vare",
            "gamma",
            "start_pi",
        ],
        value=[
            string(seed),
            string(chain_length),
            string(burnin),
            string(estimate_pi),
            string(start_h2),
            string(start_sigma_sq),
            string(start_vare),
            join(string.(gamma), ","),
            join(string.(start_pi), ","),
        ],
    )
    CSV.write(joinpath(outdir, "config.csv"), config_df)

    initial_state_df = DataFrame(
        marker_id=marker_ids,
        beta0=Float64.(initial_state.beta0),
        delta0=Int.(initial_state.delta0),
    )
    CSV.write(joinpath(outdir, "initial_state.csv"), initial_state_df)

    initial_scalars_df = DataFrame(
        key=["mu0", "sigmaSq0", "vare0"],
        value=[string(initial_state.mu0), string(start_sigma_sq), string(start_vare)],
    )
    CSV.write(joinpath(outdir, "initial_scalars.csv"), initial_scalars_df)
end

function write_jwas_parity_summary(output, outdir; sigma_sq, pi_values=nothing)
    mkpath(outdir)

    marker_effects = output["marker effects geno"]
    residual_variance = output["residual variance"]

    scalar_metrics = DataFrame(
        metric=["sigmaSq", "residual_variance", "mean_nonzero_frequency"],
        value=[
            sigma_sq,
            residual_variance[1, :Estimate],
            mean(marker_effects[!, :Model_Frequency]),
        ],
    )
    CSV.write(joinpath(outdir, "scalar_metrics.csv"), scalar_metrics)

    if haskey(output, "pi_geno")
        pi_summary = output["pi_geno"]
        pi_out = DataFrame(
            class=normalize_pi_class_labels(pi_summary[!, :π]),
            estimate=pi_summary[!, :Estimate],
        )
    elseif pi_values !== nothing
        pi_out = DataFrame(
            class=["class$(i)" for i in eachindex(pi_values)],
            estimate=Float64.(pi_values),
        )
    else
        error("write_jwas_parity_summary requires either output[\"pi_geno\"] or pi_values.")
    end
    CSV.write(joinpath(outdir, "pi.csv"), pi_out)

    marker_out = DataFrame(
        marker_id=marker_effects[!, :Marker_ID],
        estimate=marker_effects[!, :Estimate],
        model_frequency=marker_effects[!, :Model_Frequency],
    )
    CSV.write(joinpath(outdir, "marker_effects.csv"), marker_out)
end

function read_parity_summary(outdir)
    return (
        scalar_metrics = CSV.read(joinpath(outdir, "scalar_metrics.csv"), DataFrame),
        pi = CSV.read(joinpath(outdir, "pi.csv"), DataFrame),
        marker_effects = CSV.read(joinpath(outdir, "marker_effects.csv"), DataFrame),
    )
end

function write_parity_trace(path, trace_df)
    mkpath(dirname(path))
    CSV.write(path, trace_df)
end

function read_parity_trace(path)
    CSV.read(path, DataFrame)
end

function read_parity_initial_state(datadir)
    initial_state_df = CSV.read(joinpath(datadir, "initial_state.csv"), DataFrame)
    initial_scalars_df = CSV.read(joinpath(datadir, "initial_scalars.csv"), DataFrame)
    scalar_map = Dict(String(row.key) => parse(Float64, string(row.value)) for row in eachrow(initial_scalars_df))
    return (
        marker_id=String.(initial_state_df.marker_id),
        beta0=Float64.(initial_state_df.beta0),
        delta0=Int.(initial_state_df.delta0),
        mu0=scalar_map["mu0"],
        sigmaSq0=scalar_map["sigmaSq0"],
        vare0=scalar_map["vare0"],
    )
end

function compare_scalar_metrics(jwas_scalars, ref_scalars; atol_map=Dict{String, Float64}(), rtol_map=Dict{String, Float64}())
    merged = innerjoin(
        rename(jwas_scalars, :value => :jwas_value),
        rename(ref_scalars, :value => :ref_value),
        on=:metric,
    )
    merged.abs_diff = abs.(merged.jwas_value .- merged.ref_value)
    merged.rel_diff = merged.abs_diff ./ max.(abs.(merged.ref_value), eps(Float64))
    merged.atol = [get(atol_map, String(metric), NaN) for metric in merged.metric]
    merged.rtol = [get(rtol_map, String(metric), NaN) for metric in merged.metric]
    return merged
end

function compare_pi(jwas_pi, ref_pi; atol=NaN)
    merged = innerjoin(
        rename(jwas_pi, :estimate => :jwas_estimate),
        rename(ref_pi, :estimate => :ref_estimate),
        on=:class,
    )
    merged.abs_diff = abs.(merged.jwas_estimate .- merged.ref_estimate)
    merged.atol = fill(atol, nrow(merged))
    return merged
end

function compare_marker_effects(jwas_markers, ref_markers; top_k=10)
    merged = innerjoin(
        rename(jwas_markers, :estimate => :jwas_estimate, :model_frequency => :jwas_model_frequency),
        rename(ref_markers, :estimate => :ref_estimate, :model_frequency => :ref_model_frequency),
        on=:marker_id,
    )
    merged.abs_diff = abs.(merged.jwas_estimate .- merged.ref_estimate)
    merged.model_frequency_abs_diff = abs.(merged.jwas_model_frequency .- merged.ref_model_frequency)
    sort!(merged, :abs_diff, rev=true)
    return first(merged, min(top_k, nrow(merged)))
end

function compare_parity_summaries(jwas_scalars, ref_scalars, jwas_pi, ref_pi, jwas_markers, ref_markers)
    scalar_report = compare_scalar_metrics(jwas_scalars, ref_scalars)
    pi_report = compare_pi(jwas_pi, ref_pi)
    marker_report = innerjoin(
        rename(jwas_markers, :estimate => :jwas_estimate, :model_frequency => :jwas_model_frequency),
        rename(ref_markers, :estimate => :ref_estimate, :model_frequency => :ref_model_frequency),
        on=:marker_id,
    )
    marker_report.abs_diff = abs.(marker_report.jwas_estimate .- marker_report.ref_estimate)
    marker_report.model_frequency_abs_diff = abs.(marker_report.jwas_model_frequency .- marker_report.ref_model_frequency)
    marker_correlation = cor(marker_report.jwas_estimate, marker_report.ref_estimate)

    return (
        scalar_report=scalar_report,
        pi_report=pi_report,
        marker_report=marker_report,
        marker_correlation=marker_correlation,
    )
end
