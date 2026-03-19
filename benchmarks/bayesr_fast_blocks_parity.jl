using CSV
using DataFrames
using Statistics
using JWAS

include(joinpath(@__DIR__, "bayesr_parity_common.jl"))

const DEFAULT_GAMMA = Float64[0.0, 0.01, 0.1, 1.0]
const DEFAULT_START_PI = Float64[0.95, 0.03, 0.015, 0.005]

function env_int(name, default)
    return parse(Int, get(ENV, name, string(default)))
end

function env_float(name, default)
    return parse(Float64, get(ENV, name, string(default)))
end

function env_string(name, default)
    return String(get(ENV, name, default))
end

function env_fast_blocks(name, default)
    raw = get(ENV, name, default)
    lower = lowercase(String(raw))
    if lower == "true"
        return true
    elseif lower == "false"
        return false
    else
        return parse(Int, lower)
    end
end

function parse_seed_list(raw::AbstractString)
    isempty(strip(raw)) && error("Seed list must not be empty.")
    return parse.(Int, strip.(split(raw, ",")))
end

function resolve_fast_block_size(fast_blocks, n_obs)
    if fast_blocks == false
        return 1
    elseif fast_blocks == true
        return Int(floor(sqrt(n_obs)))
    elseif fast_blocks isa Number
        return Int(floor(fast_blocks))
    else
        error("Unsupported fast_blocks setting for benchmark: $fast_blocks")
    end
end

function block_aware_burnin(chain_length, burnin, fast_blocks, n_obs)
    block_size = resolve_fast_block_size(fast_blocks, n_obs)
    adjusted_chain_length = max(1, Int(floor(chain_length / block_size)))
    adjusted_burnin = fast_blocks == false ? burnin : Int(floor(burnin / block_size))
    return min(adjusted_burnin, adjusted_chain_length - 1)
end

function run_bayesr_case(case_outdir;
                         data,
                         mcmc_seed,
                         chain_length,
                         burnin,
                         start_h2,
                         start_sigma_sq,
                         start_vare,
                         fast_blocks,
                         estimate_pi,
                         estimate_marker_variance,
                         fixed_hyperparameters)
    n_obs = length(data.ids)
    run_burnin = block_aware_burnin(chain_length, burnin, fast_blocks, n_obs)
    block_size = resolve_fast_block_size(fast_blocks, n_obs)

    datadir = joinpath(case_outdir, "data")
    write_parity_dataset(datadir;
                         ids=data.ids,
                         marker_ids=data.marker_ids,
                         X=data.X,
                         y=data.y,
                         gamma=DEFAULT_GAMMA,
                         start_pi=DEFAULT_START_PI,
                         estimate_pi=estimate_pi,
                         chain_length=chain_length,
                         burnin=burnin,
                         start_h2=start_h2,
                         start_sigma_sq=start_sigma_sq,
                         start_vare=start_vare,
                         seed=mcmc_seed)

    geno_path = joinpath(datadir, "genotypes.csv")
    pheno_path = joinpath(datadir, "phenotypes.csv")
    phenotypes = CSV.read(pheno_path, DataFrame)

    global geno = get_genotypes(geno_path, start_sigma_sq,
                                separator=',',
                                method="BayesR",
                                Pi=DEFAULT_START_PI,
                                estimatePi=estimate_pi,
                                G_is_marker_variance=true,
                                estimate_variance=estimate_marker_variance,
                                estimate_scale=false,
                                quality_control=false,
                                center=false)
    model = build_model("y1 = intercept + geno", start_vare)

    output_dir = joinpath(case_outdir, fast_blocks == false ? "dense_run" : "fast_blocks_run")
    isdir(output_dir) && rm(output_dir; recursive=true, force=true)
    elapsed = @elapsed output = runMCMC(model, phenotypes,
                                        chain_length=chain_length,
                                        burnin=run_burnin,
                                        output_samples_frequency=1,
                                        output_folder=output_dir,
                                        seed=mcmc_seed,
                                        printout_model_info=false,
                                        outputEBV=false,
                                        output_heritability=false,
                                        fast_blocks=fast_blocks)

    summary_dir = joinpath(case_outdir, fast_blocks == false ? "dense_summary" : "fast_blocks_summary")
    write_jwas_parity_summary(output, summary_dir;
                              sigma_sq=model.M[1].meanVara,
                              pi_values=model.M[1].π,
                              fixed_hyperparameters=fixed_hyperparameters)

    return (summary_dir=summary_dir, elapsed=elapsed, run_burnin=run_burnin, block_size=block_size)
end

function metric_value(df::DataFrame, metric_name::AbstractString)
    idx = findfirst(==(metric_name), String.(df.metric))
    idx === nothing && error("Metric $metric_name not found.")
    return Float64(df.value[idx])
end

function run_fixed_hyperparameter_mode(outdir;
                                       data_seed,
                                       mcmc_seed,
                                       n_obs,
                                       n_markers,
                                       chain_length,
                                       burnin,
                                       start_h2,
                                       block_setting)
    mkpath(outdir)
    data = build_bayesr_parity_dataset(seed=data_seed, n_obs=n_obs, n_markers=n_markers)
    vary = var(data.y)
    start_vare = vary * (1 - start_h2)
    start_sigma_sq = vary * start_h2 / (n_markers * sum(DEFAULT_GAMMA .* DEFAULT_START_PI))

    dense_run = run_bayesr_case(joinpath(outdir, "dense_case");
                                data=data,
                                mcmc_seed=mcmc_seed,
                                chain_length=chain_length,
                                burnin=burnin,
                                start_h2=start_h2,
                                start_sigma_sq=start_sigma_sq,
                                start_vare=start_vare,
                                fast_blocks=false,
                                estimate_pi=false,
                                estimate_marker_variance=false,
                                fixed_hyperparameters=true)

    block_run = run_bayesr_case(joinpath(outdir, "fast_blocks_case");
                                data=data,
                                mcmc_seed=mcmc_seed,
                                chain_length=chain_length,
                                burnin=burnin,
                                start_h2=start_h2,
                                start_sigma_sq=start_sigma_sq,
                                start_vare=start_vare,
                                fast_blocks=block_setting,
                                estimate_pi=false,
                                estimate_marker_variance=false,
                                fixed_hyperparameters=true)

    dense = read_parity_summary(dense_run.summary_dir)
    block = read_parity_summary(block_run.summary_dir)
    comparison = compare_parity_summaries(block.scalar_metrics, dense.scalar_metrics,
                                          block.pi, dense.pi,
                                          block.marker_effects, dense.marker_effects)
    top_markers = compare_marker_effects(block.marker_effects, dense.marker_effects)

    runtime_df = bayesr_runtime_report(
        dense_seconds=dense_run.elapsed,
        block_seconds=block_run.elapsed,
        requested_burnin=burnin,
        dense_burnin=dense_run.run_burnin,
        block_burnin=block_run.run_burnin,
        dense_block_size=dense_run.block_size,
        block_block_size=block_run.block_size,
    )

    CSV.write(joinpath(outdir, "comparison_scalar_metrics.csv"), comparison.scalar_report)
    CSV.write(joinpath(outdir, "comparison_pi.csv"), comparison.pi_report)
    CSV.write(joinpath(outdir, "comparison_marker_effects_top.csv"), top_markers)
    CSV.write(joinpath(outdir, "runtime.csv"), runtime_df)

    println("WROTE ", joinpath(outdir, "comparison_scalar_metrics.csv"))
    println("WROTE ", joinpath(outdir, "comparison_pi.csv"))
    println("WROTE ", joinpath(outdir, "comparison_marker_effects_top.csv"))
    println("WROTE ", joinpath(outdir, "runtime.csv"))
end

function run_dense_multiseed_mode(outdir;
                                  data_seed,
                                  seeds,
                                  n_obs,
                                  n_markers,
                                  chain_length,
                                  burnin,
                                  start_h2)
    mkpath(outdir)
    data = build_bayesr_parity_dataset(seed=data_seed, n_obs=n_obs, n_markers=n_markers)
    vary = var(data.y)
    start_vare = vary * (1 - start_h2)
    start_sigma_sq = vary * start_h2 / (n_markers * sum(DEFAULT_GAMMA .* DEFAULT_START_PI))

    rows = DataFrame(
        seed=Int[],
        residual_variance=Float64[],
        mean_nonzero_frequency=Float64[],
        marker_abs_mean=Float64[],
        seconds=Float64[],
    )

    for seed in seeds
        run_result = run_bayesr_case(joinpath(outdir, "dense_seed_$(seed)");
                                     data=data,
                                     mcmc_seed=seed,
                                     chain_length=chain_length,
                                     burnin=burnin,
                                     start_h2=start_h2,
                                     start_sigma_sq=start_sigma_sq,
                                     start_vare=start_vare,
                                     fast_blocks=false,
                                     estimate_pi=false,
                                     estimate_marker_variance=false,
                                     fixed_hyperparameters=true)
        summary = read_parity_summary(run_result.summary_dir)
        push!(rows, (
            seed,
            metric_value(summary.scalar_metrics, "residual_variance"),
            metric_value(summary.scalar_metrics, "mean_nonzero_frequency"),
            mean(abs.(summary.marker_effects.estimate)),
            run_result.elapsed,
        ))
    end

    summary_df = summarize_within_method_multiseed(rows)
    CSV.write(joinpath(outdir, "multiseed_runs.csv"), rows)
    CSV.write(joinpath(outdir, "multiseed_summary.csv"), summary_df)

    println("WROTE ", joinpath(outdir, "multiseed_runs.csv"))
    println("WROTE ", joinpath(outdir, "multiseed_summary.csv"))
end

if length(ARGS) > 1
    error("Usage: julia --project=. --startup-file=no benchmarks/bayesr_fast_blocks_parity.jl [outdir]")
end

outdir = length(ARGS) == 1 ? abspath(ARGS[1]) : joinpath(@__DIR__, "out", "bayesr_fast_blocks_parity")
mode = env_string("JWAS_BAYESR_BLOCK_MODE", "fixed_hyperparameters")
data_seed = env_int("JWAS_BAYESR_BLOCK_DATA_SEED", 2026)
seed = env_int("JWAS_BAYESR_BLOCK_SEED", 2026)
chain_length = env_int("JWAS_BAYESR_BLOCK_CHAIN_LENGTH", 300)
burnin = env_int("JWAS_BAYESR_BLOCK_BURNIN", 100)
n_obs = env_int("JWAS_BAYESR_BLOCK_N_OBS", 60)
n_markers = env_int("JWAS_BAYESR_BLOCK_N_MARKERS", 40)
start_h2 = env_float("JWAS_BAYESR_BLOCK_START_H2", 0.5)
block_setting = env_fast_blocks("JWAS_BAYESR_BLOCK_SETTING", "true")
seed_list = parse_seed_list(env_string("JWAS_BAYESR_BLOCK_SEEDS", "2026,2027,2028,2029,2030"))

if mode == "fixed_hyperparameters"
    run_fixed_hyperparameter_mode(outdir;
                                  data_seed=data_seed,
                                  mcmc_seed=seed,
                                  n_obs=n_obs,
                                  n_markers=n_markers,
                                  chain_length=chain_length,
                                  burnin=burnin,
                                  start_h2=start_h2,
                                  block_setting=block_setting)
elseif mode == "dense_multiseed"
    run_dense_multiseed_mode(outdir;
                             data_seed=data_seed,
                             seeds=seed_list,
                             n_obs=n_obs,
                             n_markers=n_markers,
                             chain_length=chain_length,
                             burnin=burnin,
                             start_h2=start_h2)
else
    error("Unsupported benchmark mode: $mode")
end
