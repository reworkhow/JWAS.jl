using CSV
using DataFrames
using Printf
using Statistics

include(joinpath(@__DIR__, "bayesr_parity_common.jl"))

if length(ARGS) != 2
    error("Usage: julia --project=. --startup-file=no benchmarks/bayesr_parity_compare.jl <outdir> <fixed_pi|estimate_pi>")
end

outdir = ARGS[1]
mode = ARGS[2]
mode in ("fixed_pi", "estimate_pi") || error("Mode must be fixed_pi or estimate_pi.")

jwas_dir = joinpath(outdir, mode == "estimate_pi" ? "jwas_estimate_pi" : "jwas_fixed_pi")
ref_dir = joinpath(outdir, mode == "estimate_pi" ? "ref_estimate_pi" : "ref_fixed_pi")

jwas = read_parity_summary(jwas_dir)
ref = read_parity_summary(ref_dir)
report = compare_parity_summaries(jwas.scalar_metrics, ref.scalar_metrics, jwas.pi, ref.pi, jwas.marker_effects, ref.marker_effects)

CSV.write(joinpath(outdir, "comparison_scalar_metrics.csv"), report.scalar_report)
CSV.write(joinpath(outdir, "comparison_pi.csv"), report.pi_report)
CSV.write(joinpath(outdir, "comparison_marker_effects.csv"), report.marker_report)

@printf("max scalar relative diff: %.6f\n", maximum(report.scalar_report.rel_diff))
@printf("max pi abs diff: %.6f\n", maximum(report.pi_report.abs_diff))
@printf("marker-effect correlation: %.6f\n", report.marker_correlation)
