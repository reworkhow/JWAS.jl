using CSV
using DataFrames
using Printf

include(joinpath(@__DIR__, "bayesr_parity_common.jl"))

length(ARGS) == 1 || error("Usage: julia --project=. --startup-file=no benchmarks/bayesr_parity_replay_compare.jl <outdir>")

outdir = ARGS[1]
jwas_dir = joinpath(outdir, "jwas_fixed_pi")
ref_dir = joinpath(outdir, "ref_fixed_pi")

jwas_markers = CSV.read(joinpath(jwas_dir, "replay_marker_iteration1.csv"), DataFrame)
ref_markers = CSV.read(joinpath(ref_dir, "replay_marker_iteration1.csv"), DataFrame)
jwas_scalars = CSV.read(joinpath(jwas_dir, "replay_scalar_iteration1.csv"), DataFrame)
ref_scalars = CSV.read(joinpath(ref_dir, "replay_scalar_iteration1.csv"), DataFrame)

marker_report = compare_replay_marker_tables(jwas_markers, ref_markers)
scalar_report = compare_replay_scalar_tables(jwas_scalars, ref_scalars)

CSV.write(joinpath(outdir, "comparison_replay_marker_iteration1.csv"), marker_report)
CSV.write(joinpath(outdir, "comparison_replay_scalar_iteration1.csv"), scalar_report)

marker_diff_cols = filter(name -> endswith(name, "_abs_diff"), names(marker_report))
scalar_max = maximum(scalar_report.abs_diff)
marker_max = isempty(marker_diff_cols) ? 0.0 : maximum(maximum(marker_report[!, col]) for col in marker_diff_cols)

@printf("max scalar abs diff: %.12f\n", scalar_max)
@printf("max marker abs diff: %.12f\n", marker_max)

if !isempty(marker_diff_cols)
    for col in marker_diff_cols
        first_gap = findfirst(>(0), marker_report[!, col])
        if first_gap !== nothing
            @printf("first %s mismatch: marker %s\n", col, marker_report.marker_id[first_gap])
            break
        end
    end
end

first_scalar_gap = findfirst(>(0), scalar_report.abs_diff)
if first_scalar_gap !== nothing
    @printf("first scalar mismatch: %s\n", scalar_report.field[first_scalar_gap])
end
