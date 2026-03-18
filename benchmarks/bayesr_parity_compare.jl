using CSV
using DataFrames
using Printf
using Statistics

include(joinpath(@__DIR__, "bayesr_parity_common.jl"))

if !(length(ARGS) in (2, 3))
    error("Usage: julia --project=. --startup-file=no benchmarks/bayesr_parity_compare.jl <outdir> <fixed_pi|estimate_pi> [summary|trace]")
end

outdir = ARGS[1]
mode = ARGS[2]
mode in ("fixed_pi", "estimate_pi") || error("Mode must be fixed_pi or estimate_pi.")
run_kind = length(ARGS) == 3 ? ARGS[3] : "summary"
run_kind in ("summary", "trace") || error("Run kind must be summary or trace.")

jwas_dir = joinpath(outdir, mode == "estimate_pi" ? "jwas_estimate_pi" : "jwas_fixed_pi")
ref_dir = joinpath(outdir, mode == "estimate_pi" ? "ref_estimate_pi" : "ref_fixed_pi")

if run_kind == "trace"
    mode == "fixed_pi" || error("Trace comparison currently supports fixed_pi only.")
    jwas_trace = read_parity_trace(joinpath(jwas_dir, "trace_fixed_pi.csv"))
    ref_trace = read_parity_trace(joinpath(ref_dir, "trace_fixed_pi.csv"))
    report = compare_trace_metrics(jwas_trace, ref_trace)
    CSV.write(joinpath(outdir, "comparison_trace_fixed_pi.csv"), report)

    first_sigma_gap = findfirst(>(0), report.sigmaSq_abs_diff)
    @printf("first sigmaSq gap iteration: %s\n", isnothing(first_sigma_gap) ? "none" : string(report.iter[first_sigma_gap]))
    @printf("max sigmaSq abs diff: %.6f\n", maximum(report.sigmaSq_abs_diff))
    @printf("max ssq abs diff: %.6f\n", maximum(report.ssq_abs_diff))
    @printf("max nnz abs diff: %.0f\n", maximum(report.nnz_abs_diff))
    @printf("max vare abs diff: %.6f\n", maximum(report.vare_abs_diff))
else
    jwas = read_parity_summary(jwas_dir)
    ref = read_parity_summary(ref_dir)
    report = compare_parity_summaries(jwas.scalar_metrics, ref.scalar_metrics, jwas.pi, ref.pi, jwas.marker_effects, ref.marker_effects)

    CSV.write(joinpath(outdir, "comparison_scalar_metrics.csv"), report.scalar_report)
    CSV.write(joinpath(outdir, "comparison_pi.csv"), report.pi_report)
    CSV.write(joinpath(outdir, "comparison_marker_effects.csv"), report.marker_report)

    @printf("max scalar relative diff: %.6f\n", maximum(report.scalar_report.rel_diff))
    @printf("max pi abs diff: %.6f\n", maximum(report.pi_report.abs_diff))
    @printf("marker-effect correlation: %.6f\n", report.marker_correlation)
end
