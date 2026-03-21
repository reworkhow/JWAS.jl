using LinearAlgebra
using Random

include(joinpath(@__DIR__, "..", "bayesr_fast_blocks_parity.jl"))

function initialize_local_bayesr_state(data;
                                       start_sigma_sq,
                                       fast_blocks)
    tmpdir = mktempdir()
    geno_path = joinpath(tmpdir, "genotypes.csv")
    geno_df = DataFrame(ID=data.ids)
    for (j, marker_id) in enumerate(data.marker_ids)
        geno_df[!, marker_id] = data.X[:, j]
    end
    CSV.write(geno_path, geno_df)

    geno = get_genotypes(geno_path, start_sigma_sq,
                         separator=',',
                         method="BayesR",
                         Pi=DEFAULT_START_PI,
                         estimatePi=false,
                         G_is_marker_variance=true,
                         estimate_variance=false,
                         estimate_scale=false,
                         quality_control=false,
                         center=false)

    T = eltype(geno.genotypes)
    block_starts = trace_fast_block_starts(fast_blocks, geno.nObs, geno.nMarkers)
    mGibbs = JWAS.GibbsMats(geno.genotypes, ones(T, geno.nObs), fast_blocks=block_starts)
    geno.mArray = mGibbs.xArray
    geno.mRinvArray = mGibbs.xRinvArray
    geno.mpRinvm = mGibbs.xpRinvx
    geno.MArray = mGibbs.XArray
    geno.MRinvArray = mGibbs.XRinvArray
    geno.MpRinvM = mGibbs.XpRinvX
    geno.α = [zeros(T, geno.nMarkers)]
    geno.δ = [ones(Int, geno.nMarkers)]
    geno.G.val = T(start_sigma_sq)

    return geno, T, block_starts
end

function run_local_bayesr_trace(data;
                                mcmc_seed,
                                chain_length,
                                start_sigma_sq,
                                start_vare,
                                fast_blocks)
    geno, T, block_starts = initialize_local_bayesr_state(data;
                                                           start_sigma_sq=start_sigma_sq,
                                                           fast_blocks=fast_blocks)
    Random.seed!(mcmc_seed)

    y = T.(data.y)
    mu = T(mean(data.y))
    vare = T(start_vare)
    ycorr = T.(y .- mu)
    invweights = ones(T, length(ycorr))
    nue = 4.0
    scalee = start_vare * (nue - 2.0) / nue

    rows = NamedTuple[]
    for iter in 1:chain_length
        ycorr .+= mu
        rhs_mu = sum(ycorr)
        inv_lhs_mu = one(T) / length(y)
        mu_hat = inv_lhs_mu * rhs_mu
        mu = mu_hat + randn() * sqrt(inv_lhs_mu * vare)
        ycorr .-= mu

        if fast_blocks === false
            JWAS.BayesR!(geno, ycorr, vare)
        else
            JWAS.BayesR_block!(geno, ycorr, vare, invweights)
        end

        vare = T(JWAS.sample_variance(ycorr, length(ycorr), nue, scalee, false))
        alpha = geno.α[1]
        push!(rows, (
            iter=iter,
            residual_variance=Float64(vare),
            ycorr_norm=Float64(norm(ycorr)),
            alpha_norm=Float64(norm(alpha)),
            alpha_abs_mean=Float64(mean(abs.(alpha))),
            nnz=count(>(1), geno.δ[1]),
            max_abs_alpha=Float64(maximum(abs.(alpha))),
        ))
    end

    return DataFrame(rows), block_starts
end

function local_bayesr_output(marker_ids;
                             mean_sigma_sq,
                             mean_vare,
                             mean_pi,
                             mean_alpha,
                             mean_model_frequency)
    return Dict(
        "pi_geno" => DataFrame(
            π=["class$(i)" for i in eachindex(mean_pi)],
            Estimate=mean_pi,
            SD=zeros(length(mean_pi)),
        ),
        "marker effects geno" => DataFrame(
            Trait=fill("y1", length(marker_ids)),
            Marker_ID=marker_ids,
            Estimate=mean_alpha,
            SD=zeros(length(marker_ids)),
            Model_Frequency=mean_model_frequency,
        ),
        "residual variance" => DataFrame(
            Covariance=["y1_y1"],
            Estimate=[mean_vare],
            SD=[0.0],
        ),
        "sigmaSq" => mean_sigma_sq,
    )
end

function run_local_bayesr_summary(data,
                                  mcmc_seed,
                                  chain_length,
                                  burnin,
                                  start_sigma_sq,
                                  start_vare,
                                  fast_blocks;
                                  force_single_rep=false)
    chain_length > burnin || error("chain_length must exceed burnin for posterior summaries.")
    geno, T, block_starts = initialize_local_bayesr_state(data;
                                                           start_sigma_sq=start_sigma_sq,
                                                           fast_blocks=fast_blocks)
    Random.seed!(mcmc_seed)

    y = T.(data.y)
    mu = T(mean(data.y))
    vare = T(start_vare)
    ycorr = T.(y .- mu)
    invweights = ones(T, length(ycorr))
    nue = 4.0
    scalee = start_vare * (nue - 2.0) / nue

    posterior_count = 0
    sigma_sq_sum = 0.0
    vare_sum = 0.0
    pi_sum = zeros(Float64, length(geno.π))
    alpha_sum = zeros(Float64, geno.nMarkers)
    model_frequency_sum = zeros(Float64, geno.nMarkers)

    elapsed = @elapsed begin
        for iter in 1:chain_length
            ycorr .+= mu
            rhs_mu = sum(ycorr)
            inv_lhs_mu = one(T) / length(y)
            mu_hat = inv_lhs_mu * rhs_mu
            mu = mu_hat + randn() * sqrt(inv_lhs_mu * vare)
            ycorr .-= mu

            if fast_blocks === false
                JWAS.BayesR!(geno, ycorr, vare)
            elseif force_single_rep
                JWAS.BayesR_block!(geno, ycorr, vare, invweights, iter, typemax(Int))
            else
                JWAS.BayesR_block!(geno, ycorr, vare, invweights)
            end

            geno.π .= JWAS.samplePi(geno.δ[1], length(geno.π))
            JWAS.sample_marker_effect_variance(geno)
            vare = T(JWAS.sample_variance(ycorr, length(ycorr), nue, scalee, invweights))

            if iter > burnin
                posterior_count += 1
                sigma_sq_sum += Float64(geno.G.val)
                vare_sum += Float64(vare)
                pi_sum .+= Float64.(geno.π)
                alpha_sum .+= Float64.(geno.α[1])
                model_frequency_sum .+= Float64.(geno.δ[1] .> 1)
            end
        end
    end

    output = local_bayesr_output(
        data.marker_ids;
        mean_sigma_sq=sigma_sq_sum / posterior_count,
        mean_vare=vare_sum / posterior_count,
        mean_pi=pi_sum ./ posterior_count,
        mean_alpha=alpha_sum ./ posterior_count,
        mean_model_frequency=model_frequency_sum ./ posterior_count,
    )

    return (output=output, elapsed=elapsed, block_starts=block_starts)
end

function run_fixed_hyperparameter_trace_mode(outdir;
                                             data_seed,
                                             mcmc_seed,
                                             n_obs,
                                             n_markers,
                                             chain_length,
                                             burnin,
                                             start_h2,
                                             block_setting)
    burnin == 0 || error("fixed_hyperparameters_trace mode requires burnin=0.")
    mkpath(outdir)
    data = build_bayesr_parity_dataset(seed=data_seed, n_obs=n_obs, n_markers=n_markers)
    vary = var(data.y)
    start_vare = vary * (1 - start_h2)
    start_sigma_sq = vary * start_h2 / (n_markers * sum(DEFAULT_GAMMA .* DEFAULT_START_PI))

    datadir = joinpath(outdir, "data")
    write_parity_dataset(datadir;
                         ids=data.ids,
                         marker_ids=data.marker_ids,
                         X=data.X,
                         y=data.y,
                         gamma=DEFAULT_GAMMA,
                         start_pi=DEFAULT_START_PI,
                         estimate_pi=false,
                         chain_length=chain_length,
                         burnin=burnin,
                         start_h2=start_h2,
                         start_sigma_sq=start_sigma_sq,
                         start_vare=start_vare,
                         seed=mcmc_seed)

    dense_trace, _ = run_local_bayesr_trace(data;
                                            mcmc_seed=mcmc_seed,
                                            chain_length=chain_length,
                                            start_sigma_sq=start_sigma_sq,
                                            start_vare=start_vare,
                                            fast_blocks=false)
    block_trace, block_starts = run_local_bayesr_trace(data;
                                                       mcmc_seed=mcmc_seed,
                                                       chain_length=chain_length,
                                                       start_sigma_sq=start_sigma_sq,
                                                       start_vare=start_vare,
                                                       fast_blocks=block_setting)
    comparison = compare_trace_metrics(block_trace, dense_trace)
    if block_starts != false
        comparison.block_size = fill(resolve_fast_block_size(block_setting, n_obs), nrow(comparison))
    end

    CSV.write(joinpath(outdir, "dense_trace.csv"), dense_trace)
    CSV.write(joinpath(outdir, "fast_blocks_trace.csv"), block_trace)
    CSV.write(joinpath(outdir, "trace_comparison.csv"), comparison)

    println("WROTE ", joinpath(outdir, "dense_trace.csv"))
    println("WROTE ", joinpath(outdir, "fast_blocks_trace.csv"))
    println("WROTE ", joinpath(outdir, "trace_comparison.csv"))
end

function run_default_blocks_single_rep_mode(outdir;
                                            data_seed,
                                            mcmc_seed,
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

    dense_run = run_local_bayesr_summary(data,
                                         mcmc_seed,
                                         chain_length,
                                         burnin,
                                         start_sigma_sq,
                                         start_vare,
                                         false)
    block_run = run_local_bayesr_summary(data,
                                         mcmc_seed,
                                         chain_length,
                                         burnin,
                                         start_sigma_sq,
                                         start_vare,
                                         true;
                                         force_single_rep=true)

    dense_summary_dir = joinpath(outdir, "dense_summary")
    block_summary_dir = joinpath(outdir, "default_blocks_single_rep_summary")
    write_jwas_parity_summary(dense_run.output, dense_summary_dir;
                              sigma_sq=dense_run.output["sigmaSq"])
    write_jwas_parity_summary(block_run.output, block_summary_dir;
                              sigma_sq=block_run.output["sigmaSq"])

    dense = read_parity_summary(dense_summary_dir)
    block = read_parity_summary(block_summary_dir)
    comparison = compare_parity_summaries(block.scalar_metrics, dense.scalar_metrics,
                                          block.pi, dense.pi,
                                          block.marker_effects, dense.marker_effects)
    top_markers = compare_marker_effects(block.marker_effects, dense.marker_effects)

    runtime_df = bayesr_runtime_report(
        dense_seconds=dense_run.elapsed,
        block_seconds=block_run.elapsed,
        requested_burnin=burnin,
        dense_burnin=burnin,
        block_burnin=burnin,
        dense_block_size=1,
        block_block_size=resolve_fast_block_size(true, n_obs),
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

function main_debug(argv=ARGS)
    length(argv) > 1 && error("Usage: julia --project=. --startup-file=no benchmarks/debug/bayesr_fast_blocks_debug.jl [outdir]")

    outdir = length(argv) == 1 ? abspath(argv[1]) : joinpath(@__DIR__, "..", "out", "bayesr_fast_blocks_debug")
    mode = env_string("JWAS_BAYESR_BLOCK_MODE", "fixed_hyperparameters_trace")
    data_seed = env_int("JWAS_BAYESR_BLOCK_DATA_SEED", 2026)
    seed = env_int("JWAS_BAYESR_BLOCK_SEED", 2026)
    default_chain_length = mode == "fixed_hyperparameters_trace" ? 100 : 300
    default_burnin = mode == "fixed_hyperparameters_trace" ? 0 : 100
    chain_length = env_int("JWAS_BAYESR_BLOCK_CHAIN_LENGTH", default_chain_length)
    burnin = env_int("JWAS_BAYESR_BLOCK_BURNIN", default_burnin)
    n_obs = env_int("JWAS_BAYESR_BLOCK_N_OBS", 60)
    n_markers = env_int("JWAS_BAYESR_BLOCK_N_MARKERS", 40)
    start_h2 = env_float("JWAS_BAYESR_BLOCK_START_H2", 0.5)
    block_setting = env_fast_blocks("JWAS_BAYESR_BLOCK_SETTING", "true")

    if mode == "fixed_hyperparameters_trace"
        run_fixed_hyperparameter_trace_mode(outdir;
                                            data_seed=data_seed,
                                            mcmc_seed=seed,
                                            n_obs=n_obs,
                                            n_markers=n_markers,
                                            chain_length=chain_length,
                                            burnin=burnin,
                                            start_h2=start_h2,
                                            block_setting=block_setting)
    elseif mode == "default_blocks_single_rep"
        run_default_blocks_single_rep_mode(outdir;
                                           data_seed=data_seed,
                                           mcmc_seed=seed,
                                           n_obs=n_obs,
                                           n_markers=n_markers,
                                           chain_length=chain_length,
                                           burnin=burnin,
                                           start_h2=start_h2)
    else
        error("Unsupported debug benchmark mode: $mode")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_debug()
end
