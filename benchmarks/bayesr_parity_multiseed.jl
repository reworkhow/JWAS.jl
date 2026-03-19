using CSV
using DataFrames
using Statistics

include(joinpath(@__DIR__, "bayesr_parity_common.jl"))

function parse_seeds(raw::AbstractString)
    parse.(Int, split(raw, ","))
end

function parity_seed_row(seed_outdir, mode, seed)
    jwas_dir = joinpath(seed_outdir, mode == "estimate_pi" ? "jwas_estimate_pi" : "jwas_fixed_pi")
    ref_dir = joinpath(seed_outdir, mode == "estimate_pi" ? "ref_estimate_pi" : "ref_fixed_pi")

    jwas = read_parity_summary(jwas_dir)
    ref = read_parity_summary(ref_dir)

    jwas_scalars = Dict(String(row.metric) => Float64(row.value) for row in eachrow(jwas.scalar_metrics))
    ref_scalars = Dict(String(row.metric) => Float64(row.value) for row in eachrow(ref.scalar_metrics))

    row = (
        seed=seed,
        jwas_sigmaSq=jwas_scalars["sigmaSq"],
        ref_sigmaSq=ref_scalars["sigmaSq"],
        sigma_rel_diff=abs(jwas_scalars["sigmaSq"] - ref_scalars["sigmaSq"]) / abs(ref_scalars["sigmaSq"]),
        jwas_vare=jwas_scalars["residual_variance"],
        ref_vare=ref_scalars["residual_variance"],
        vare_rel_diff=abs(jwas_scalars["residual_variance"] - ref_scalars["residual_variance"]) / abs(ref_scalars["residual_variance"]),
        jwas_nonzero=jwas_scalars["mean_nonzero_frequency"],
        ref_nonzero=ref_scalars["mean_nonzero_frequency"],
        nonzero_abs_diff=abs(jwas_scalars["mean_nonzero_frequency"] - ref_scalars["mean_nonzero_frequency"]),
    )

    if mode == "estimate_pi"
        pi_cmp = compare_pi(jwas.pi, ref.pi)
        return merge(row, (max_pi_abs_diff=maximum(pi_cmp.abs_diff),))
    end

    return row
end

if !(length(ARGS) in (2, 3))
    error("Usage: julia --project=. --startup-file=no benchmarks/bayesr_parity_multiseed.jl <outdir> <fixed_pi|estimate_pi> [seed1,seed2,...]")
end

outdir = abspath(ARGS[1])
mode = ARGS[2]
mode in ("fixed_pi", "estimate_pi") || error("Mode must be fixed_pi or estimate_pi.")
seeds = length(ARGS) == 3 ? parse_seeds(ARGS[3]) : collect(2021:2030)
chain_length = parse(Int, get(ENV, "JWAS_PARITY_CHAIN_LENGTH", "10000"))
burnin = parse(Int, get(ENV, "JWAS_PARITY_BURNIN", "2000"))
repo_root = dirname(@__DIR__)
jwas_script = joinpath(@__DIR__, "bayesr_parity_jwas.jl")
r_script = joinpath(@__DIR__, "bayesr_parity_reference.R")

mkpath(outdir)
runs = NamedTuple[]

for seed in seeds
    seed_outdir = joinpath(outdir, "seed_$(seed)")
    isdir(seed_outdir) && rm(seed_outdir; recursive=true, force=true)

    env = copy(ENV)
    env["JWAS_PARITY_OUTDIR"] = seed_outdir
    env["JWAS_PARITY_MODE"] = mode
    env["JWAS_PARITY_CHAIN_LENGTH"] = string(chain_length)
    env["JWAS_PARITY_BURNIN"] = string(burnin)
    env["JWAS_PARITY_SEED"] = string(seed)

    run(Cmd(setenv(`$(Base.julia_cmd()) --project=. --startup-file=no $jwas_script`, env); dir=repo_root))
    run(Cmd(`Rscript $r_script $seed_outdir $mode`; dir=repo_root))

    push!(runs, parity_seed_row(seed_outdir, mode, seed))
end

runs_df = DataFrame(runs)
summary_df = summarize_multiseed_parity(runs_df)

sigma_scale_df = DataFrame(
    source=["JWAS", "R"],
    mean_sigmaSq=[mean(runs_df.jwas_sigmaSq), mean(runs_df.ref_sigmaSq)],
    sd_sigmaSq=[std(runs_df.jwas_sigmaSq), std(runs_df.ref_sigmaSq)],
)

CSV.write(joinpath(outdir, "multiseed_runs.csv"), runs_df)
CSV.write(joinpath(outdir, "multiseed_summary.csv"), summary_df)
CSV.write(joinpath(outdir, "multiseed_sigma_scale.csv"), sigma_scale_df)

println("WROTE ", joinpath(outdir, "multiseed_runs.csv"))
println("WROTE ", joinpath(outdir, "multiseed_summary.csv"))
println("WROTE ", joinpath(outdir, "multiseed_sigma_scale.csv"))
