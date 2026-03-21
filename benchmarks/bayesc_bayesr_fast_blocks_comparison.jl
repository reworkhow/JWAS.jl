using CSV
using DataFrames
using Statistics
using JWAS

include(joinpath(@__DIR__, "bayesr_parity_common.jl"))

const BENCH_BAYESR_GAMMA = Float64[0.0, 0.01, 0.1, 1.0]
const BENCH_BAYESR_START_PI = Float64[0.95, 0.03, 0.015, 0.005]

function env_int(name, default)
    return parse(Int, get(ENV, name, string(default)))
end

function env_float(name, default)
    return parse(Float64, get(ENV, name, string(default)))
end

function env_string(name, default)
    return String(get(ENV, name, default))
end

function parse_seed_list(raw::AbstractString)
    isempty(strip(raw)) && error("Seed list must not be empty.")
    return parse.(Int, strip.(split(raw, ",")))
end

function write_method_block_dataset(outdir; ids, marker_ids, X, y)
    mkpath(outdir)
    geno_df = DataFrame(ID=ids)
    for (j, marker_id) in enumerate(marker_ids)
        geno_df[!, marker_id] = X[:, j]
    end
    CSV.write(joinpath(outdir, "genotypes.csv"), geno_df)
    CSV.write(joinpath(outdir, "phenotypes.csv"), DataFrame(ID=ids, y1=y))
end

function resolve_fast_block_size(fast_blocks, n_obs)
    if fast_blocks === false
        return 1
    elseif fast_blocks === true
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
    adjusted_burnin = fast_blocks === false ? burnin : Int(floor(burnin / block_size))
    return min(adjusted_burnin, adjusted_chain_length - 1)
end

function method_case_row(method::AbstractString,
                         variant::AbstractString,
                         seed::Integer,
                         output,
                         phenotypes::DataFrame,
                         model,
                         seconds::Real,
                         block_size::Integer,
                         effective_burnin::Integer)
    marker_effects = output["marker effects geno"]
    residual_variance = Float64(output["residual variance"][1, :Estimate])
    mean_model_frequency = Float64(mean(marker_effects[!, :Model_Frequency]))
    marker_abs_mean = Float64(mean(abs.(marker_effects[!, :Estimate])))
    phenotype_ebv_correlation = NaN

    if haskey(output, "EBV_y1")
        ebv_df = output["EBV_y1"]
        joined = innerjoin(phenotypes[:, [:ID, :y1]], ebv_df[:, [:ID, :EBV]], on=:ID)
        phenotype_ebv_correlation = cor(Float64.(joined[!, :y1]), Float64.(joined[!, :EBV]))
    end

    pi_scalar = NaN
    pi_class1 = NaN
    pi_class2 = NaN
    pi_class3 = NaN
    pi_class4 = NaN
    marker_variance = NaN
    sigmaSq = NaN

    if haskey(output, "pi_geno")
        pi_df = output["pi_geno"]
        if method == "BayesC"
            pi_scalar = Float64(pi_df[1, :Estimate])
            marker_variance = Float64(model.M[1].meanVara)
        elseif method == "BayesR"
            pi_lookup = Dict(
                normalize_pi_class_labels(pi_df[!, :π]) .=> Float64.(pi_df[!, :Estimate])
            )
            pi_class1 = get(pi_lookup, "class1", NaN)
            pi_class2 = get(pi_lookup, "class2", NaN)
            pi_class3 = get(pi_lookup, "class3", NaN)
            pi_class4 = get(pi_lookup, "class4", NaN)
            sigmaSq = Float64(model.M[1].meanVara)
        end
    end

    return (
        method=String(method),
        method_variant=String(variant),
        seed=Int(seed),
        residual_variance=residual_variance,
        mean_model_frequency=mean_model_frequency,
        marker_abs_mean=marker_abs_mean,
        phenotype_ebv_correlation=phenotype_ebv_correlation,
        pi_scalar=pi_scalar,
        marker_variance=marker_variance,
        sigmaSq=sigmaSq,
        pi_class1=pi_class1,
        pi_class2=pi_class2,
        pi_class3=pi_class3,
        pi_class4=pi_class4,
        seconds=Float64(seconds),
        block_size=Int(block_size),
        effective_burnin=Int(effective_burnin),
    )
end

function run_method_case(case_outdir;
                         data,
                         method,
                         variant,
                         mcmc_seed,
                         chain_length,
                         burnin,
                         start_h2,
                         fast_blocks)
    mkpath(case_outdir)
    write_method_block_dataset(case_outdir;
                               ids=data.ids,
                               marker_ids=data.marker_ids,
                               X=data.X,
                               y=data.y)

    geno_path = joinpath(case_outdir, "genotypes.csv")
    pheno_path = joinpath(case_outdir, "phenotypes.csv")
    phenotypes = CSV.read(pheno_path, DataFrame)
    vary = var(data.y)
    start_vare = vary * (1 - start_h2)
    run_burnin = block_aware_burnin(chain_length, burnin, fast_blocks, length(data.ids))
    block_size = resolve_fast_block_size(fast_blocks, length(data.ids))

    if method == "BayesR"
        start_pi = copy(BENCH_BAYESR_START_PI)
        start_sigma_sq = vary * start_h2 / (length(data.marker_ids) * sum(BENCH_BAYESR_GAMMA .* start_pi))
        global geno = get_genotypes(geno_path, start_sigma_sq,
                                    separator=',',
                                    method=method,
                                    Pi=copy(start_pi),
                                    estimatePi=true,
                                    G_is_marker_variance=true,
                                    estimate_variance=true,
                                    estimate_scale=false,
                                    quality_control=false,
                                    center=false)
        model = build_model("y1 = intercept + geno", start_vare)
        outputEBV(model, geno.obsID)
    elseif method == "BayesC"
        start_genetic_variance = vary * start_h2
        global geno = get_genotypes(geno_path, start_genetic_variance,
                                    separator=',',
                                    method=method,
                                    Pi=0.0,
                                    estimatePi=true,
                                    estimate_variance=true,
                                    estimate_scale=false,
                                    quality_control=false,
                                    center=false)
        model = build_model("y1 = intercept + geno", start_vare)
        outputEBV(model, geno.obsID)
    else
        error("Unsupported method: $method")
    end

    output_dir = joinpath(case_outdir, "run")
    isdir(output_dir) && rm(output_dir; recursive=true, force=true)
    elapsed = @elapsed output = runMCMC(model, phenotypes,
                                        chain_length=chain_length,
                                        burnin=run_burnin,
                                        output_samples_frequency=1,
                                        output_folder=output_dir,
                                        seed=mcmc_seed,
                                        printout_model_info=false,
                                        outputEBV=true,
                                        output_heritability=false,
                                        fast_blocks=fast_blocks)

    return method_case_row(method, variant, mcmc_seed, output, phenotypes, model, elapsed, block_size, run_burnin)
end

function pairwise_summary_rows(runs::DataFrame)
    rows = NamedTuple[]
    summary_metrics = (
        :residual_variance,
        :mean_model_frequency,
        :marker_abs_mean,
        :phenotype_ebv_correlation,
        :pi_scalar,
        :marker_variance,
        :sigmaSq,
        :pi_class1,
        :pi_class2,
        :pi_class3,
        :pi_class4,
    )

    for method in sort(unique(String.(runs.method)))
        dense = runs[(runs.method .== method) .& endswith.(String.(runs.method_variant), "_dense"), :]
        for variant in sort(unique(String.(runs[runs.method .== method, :method_variant])))
            endswith(variant, "_dense") && continue
            variant_runs = runs[runs.method_variant .== variant, :]
            merged = innerjoin(
                dense[:, [:seed, summary_metrics...]],
                variant_runs[:, [:seed, summary_metrics..., :seconds]],
                on=:seed,
                makeunique=true,
            )
            for metric in summary_metrics
                dense_col = metric
                variant_col = Symbol(string(metric) * "_1")
                values = Float64.(merged[!, variant_col])
                dense_values = Float64.(merged[!, dense_col])
                if all(isnan, values) || all(isnan, dense_values)
                    continue
                end
                diffs = abs.(values .- dense_values)
                push!(rows, (
                    method=method,
                    comparison=String(variant * "_vs_" * method * "_dense"),
                    metric=String(metric),
                    summary_kind="abs_diff",
                    mean_value=mean(diffs),
                    max_value=maximum(diffs),
                ))
            end

            speedups = Float64.(dense.seconds) ./ max.(Float64.(variant_runs.seconds), eps(Float64))
            push!(rows, (
                method=method,
                comparison=String(variant * "_vs_" * method * "_dense"),
                metric="speedup_vs_dense",
                summary_kind="ratio",
                mean_value=mean(speedups),
                max_value=maximum(speedups),
            ))
        end
    end

    return DataFrame(rows)
end

function main(argv=ARGS)
    length(argv) > 1 && error("Usage: julia --project=. --startup-file=no benchmarks/bayesc_bayesr_fast_blocks_comparison.jl [outdir]")

    outdir = length(argv) == 1 ? abspath(argv[1]) : joinpath(@__DIR__, "out", "bayesc_bayesr_fast_blocks_comparison")
    data_seed = env_int("JWAS_METHOD_BLOCK_DATA_SEED", 2026)
    chain_length = env_int("JWAS_METHOD_BLOCK_CHAIN_LENGTH", 10000)
    burnin = env_int("JWAS_METHOD_BLOCK_BURNIN", 2000)
    n_obs = env_int("JWAS_METHOD_BLOCK_N_OBS", 60)
    n_markers = env_int("JWAS_METHOD_BLOCK_N_MARKERS", 40)
    start_h2 = env_float("JWAS_METHOD_BLOCK_START_H2", 0.5)
    seeds = parse_seed_list(env_string("JWAS_METHOD_BLOCK_SEEDS", "2026,2027,2028,2029,2030"))

    mkpath(outdir)
    data = build_bayesr_parity_dataset(seed=data_seed, n_obs=n_obs, n_markers=n_markers)
    rows = NamedTuple[]

    cases = (
        ("BayesC", "BayesC_dense", false),
        ("BayesC", "BayesC_fast_blocks_default", true),
        ("BayesC", "BayesC_fast_blocks_1", 1),
        ("BayesR", "BayesR_dense", false),
        ("BayesR", "BayesR_fast_blocks_default", true),
        ("BayesR", "BayesR_fast_blocks_1", 1),
    )

    for seed in seeds
        for (method, variant, fast_blocks) in cases
            case_outdir = joinpath(outdir, "$(variant)_seed_$(seed)")
            push!(rows, run_method_case(case_outdir;
                                        data=data,
                                        method=method,
                                        variant=variant,
                                        mcmc_seed=seed,
                                        chain_length=chain_length,
                                        burnin=burnin,
                                        start_h2=start_h2,
                                        fast_blocks=fast_blocks))
        end
    end

    runs_df = DataFrame(rows)
    summary_df = pairwise_summary_rows(runs_df)
    CSV.write(joinpath(outdir, "comparison_runs.csv"), runs_df)
    CSV.write(joinpath(outdir, "comparison_pairwise_summary.csv"), summary_df)

    println("WROTE ", joinpath(outdir, "comparison_runs.csv"))
    println("WROTE ", joinpath(outdir, "comparison_pairwise_summary.csv"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
