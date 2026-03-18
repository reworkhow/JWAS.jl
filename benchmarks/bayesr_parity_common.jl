using CSV
using DataFrames
using Random
using Statistics

function write_jwas_parity_summary(output, outdir; sigma_sq)
    mkpath(outdir)

    marker_effects = output["marker effects geno"]
    residual_variance = output["residual variance"]
    pi_summary = output["pi_geno"]

    scalar_metrics = DataFrame(
        metric=["sigmaSq", "residual_variance", "mean_nonzero_frequency"],
        value=[
            sigma_sq,
            residual_variance[1, :Estimate],
            mean(marker_effects[!, :Model_Frequency]),
        ],
    )
    CSV.write(joinpath(outdir, "scalar_metrics.csv"), scalar_metrics)

    pi_out = DataFrame(
        class=string.(pi_summary[!, :π]),
        estimate=pi_summary[!, :Estimate],
    )
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
