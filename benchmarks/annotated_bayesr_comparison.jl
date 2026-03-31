using CSV
using DataFrames
using Distributions
using Random
using Statistics
using JWAS

include(joinpath(@__DIR__, "bayesr_parity_common.jl"))

const ANNOT_BENCH_BAYESR_GAMMA = Float64[0.0, 0.01, 0.1, 1.0]
const ANNOT_BENCH_BAYESR_START_PI = Float64[0.95, 0.03, 0.015, 0.005]

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

function safe_mean(values)
    isempty(values) && return NaN
    return mean(Float64.(values))
end

function write_benchmark_dataset(outdir; ids, marker_ids, X, y)
    mkpath(outdir)
    geno_df = DataFrame(ID=ids)
    for (j, marker_id) in enumerate(marker_ids)
        geno_df[!, marker_id] = X[:, j]
    end
    CSV.write(joinpath(outdir, "genotypes.csv"), geno_df)
    CSV.write(joinpath(outdir, "phenotypes.csv"), DataFrame(ID=ids, y1=y))
end

function joint_probs_from_conditionals(p1::Real, p2::Real, p3::Real)
    0.0 <= p1 <= 1.0 || error("p1 must lie in [0, 1].")
    0.0 <= p2 <= 1.0 || error("p2 must lie in [0, 1].")
    0.0 <= p3 <= 1.0 || error("p3 must lie in [0, 1].")
    probs = Float64[
        1 - p1,
        p1 * (1 - p2),
        p1 * p2 * (1 - p3),
        p1 * p2 * p3,
    ]
    isapprox(sum(probs), 1.0; atol=1e-12) || error("Conditional probabilities must map to a valid class distribution.")
    return probs
end

function scenario_class_probabilities(name::AbstractString)
    if name == "sparse_upper_classes"
        return (
            baseline=Float64[0.97, 0.02, 0.008, 0.002],
            enriched=Float64[0.80, 0.08, 0.07, 0.05],
        )
    elseif name == "less_sparse_upper_classes"
        return (
            baseline=Float64[0.90, 0.05, 0.03, 0.02],
            enriched=Float64[0.65, 0.15, 0.12, 0.08],
        )
    elseif name == "stepwise_annotation_signal"
        return (
            baseline=joint_probs_from_conditionals(0.10, 0.20, 0.20),
            enriched=joint_probs_from_conditionals(0.30, 0.60, 0.60),
        )
    end
    error("Unsupported annotated benchmark scenario: $name")
end

function sample_truth_classes!(rng, annotation_1, class_probs)
    n_markers = length(annotation_1)
    true_class = Vector{Int}(undef, n_markers)
    baseline_probs = class_probs.baseline
    enriched_probs = class_probs.enriched

    for j in 1:n_markers
        probs = annotation_1[j] == 1.0 ? enriched_probs : baseline_probs
        true_class[j] = rand(rng, Categorical(probs))
    end

    if !any(true_class .> 1)
        forced_idx = findfirst(==(1.0), annotation_1)
        forced_idx = forced_idx === nothing ? 1 : forced_idx
        true_class[forced_idx] = 4
    end
    if !any(true_class .== 1)
        forced_idx = findfirst(==(0.0), annotation_1)
        forced_idx = forced_idx === nothing ? 1 : forced_idx
        true_class[forced_idx] = 1
    end

    return true_class
end

function build_annotated_benchmark_dataset(; seed=20260328, n_obs=200, n_markers=1000, h2=0.45, scenario="sparse_upper_classes")
    rng = MersenneTwister(seed)
    class_probs = scenario_class_probabilities(scenario)

    ids = ["id_$(i)" for i in 1:n_obs]
    marker_ids = ["m$(j)" for j in 1:n_markers]

    X = Matrix{Float64}(undef, n_obs, n_markers)
    allele_freq = rand(rng, n_markers) .* 0.40 .+ 0.05
    for j in 1:n_markers
        p = allele_freq[j]
        for i in 1:n_obs
            X[i, j] = (rand(rng) < p) + (rand(rng) < p)
        end
    end

    n_enriched = max(1, round(Int, 0.25 * n_markers))
    annotation_1 = zeros(Float64, n_markers)
    annotation_1[randperm(rng, n_markers)[1:n_enriched]] .= 1.0

    n_nuisance = max(1, round(Int, 0.50 * n_markers))
    annotation_2 = zeros(Float64, n_markers)
    annotation_2[randperm(rng, n_markers)[1:n_nuisance]] .= 1.0

    true_class = sample_truth_classes!(rng, annotation_1, class_probs)

    raw_beta = zeros(Float64, n_markers)
    for j in 1:n_markers
        class_j = true_class[j]
        if class_j > 1
            raw_beta[j] = randn(rng) * sqrt(ANNOT_BENCH_BAYESR_GAMMA[class_j])
        end
    end

    genetic_raw = X * raw_beta
    genetic_var = var(genetic_raw)
    if !isfinite(genetic_var) || genetic_var <= 0
        error("Synthetic annotated benchmark generated zero genetic variance.")
    end

    beta_true = raw_beta / sqrt(genetic_var)
    genetic_value = X * beta_true
    residual_sd = sqrt((1 - h2) / h2)
    y = 1.0 .+ genetic_value .+ randn(rng, n_obs) .* residual_sd

    truth_metadata = DataFrame(
        scenario=fill(String(scenario), n_markers),
        marker_id=marker_ids,
        annotation_1=annotation_1,
        annotation_2=annotation_2,
        true_class=true_class,
        true_class_label=["class$(k)" for k in true_class],
        true_effect=beta_true,
        is_causal=true_class .> 1,
        is_enriched=annotation_1 .== 1.0,
    )

    return (
        scenario=String(scenario),
        ids=ids,
        marker_ids=marker_ids,
        X=X,
        y=y,
        annotations=hcat(annotation_1, annotation_2),
        truth_metadata=truth_metadata,
        target_h2=h2,
    )
end

function ebv_correlation(output, phenotypes)
    haskey(output, "EBV_y1") || return NaN
    joined = innerjoin(phenotypes[:, [:ID, :y1]], output["EBV_y1"][:, [:ID, :EBV]], on=:ID)
    return cor(Float64.(joined[!, :y1]), Float64.(joined[!, :EBV]))
end

function marker_summary_rows(method_variant::AbstractString, seed::Integer, scenario::AbstractString, merged::DataFrame)
    rows = NamedTuple[]

    group_specs = [
        ("causality", "causal", merged.is_causal .== true),
        ("causality", "noncausal", merged.is_causal .== false),
        ("annotation_1", "enriched", merged.is_enriched .== true),
        ("annotation_1", "baseline", merged.is_enriched .== false),
    ]

    for class_label in sort(unique(String.(merged.true_class_label)))
        push!(group_specs, ("true_class", class_label, String.(merged.true_class_label) .== class_label))
    end

    for (group_type, group_label, mask) in group_specs
        subset = merged[mask, :]
        isempty(subset) && continue
        push!(rows, (
            scenario=String(scenario),
            method_variant=String(method_variant),
            seed=Int(seed),
            group_type=String(group_type),
            group_label=String(group_label),
            n_markers=nrow(subset),
            mean_pip=safe_mean(subset.Model_Frequency),
            mean_abs_effect=safe_mean(abs.(subset.Estimate)),
        ))
    end

    return rows
end

function annotation_coefficient_rows(method_variant::AbstractString, seed::Integer, scenario::AbstractString, output)
    key = "annotation coefficients geno"
    haskey(output, key) || return NamedTuple[]

    coeff_df = output[key]
    steps = if "Step" in names(coeff_df)
        String.(coeff_df[!, :Step])
    else
        fill("step1_zero_vs_nonzero", nrow(coeff_df))
    end

    rows = NamedTuple[]
    for rowi in 1:nrow(coeff_df)
        push!(rows, (
            scenario=String(scenario),
            method_variant=String(method_variant),
            seed=Int(seed),
            Annotation=String(coeff_df[rowi, :Annotation]),
            Step=String(steps[rowi]),
            Estimate=Float64(coeff_df[rowi, :Estimate]),
            SD=Float64(coeff_df[rowi, :SD]),
        ))
    end
    return rows
end

function method_case_row(method_variant::AbstractString,
                         seed::Integer,
                         scenario::AbstractString,
                         output,
                         phenotypes::DataFrame,
                         truth_metadata::DataFrame,
                         model,
                         seconds::Real)
    marker_effects = output["marker effects geno"][:, [:Marker_ID, :Estimate, :Model_Frequency]]
    merged = innerjoin(
        rename(marker_effects, :Marker_ID => :marker_id),
        truth_metadata,
        on=:marker_id,
    )

    causal_mask = merged.is_causal .== true
    enriched_mask = merged.is_enriched .== true
    n_causal = sum(causal_mask)
    top_k = max(1, min(nrow(merged), n_causal))
    top_idx = sortperm(Float64.(merged.Model_Frequency), rev=true)[1:top_k]
    top_causal_recall = mean(Float64.(merged.is_causal[top_idx]))

    residual_variance = Float64(output["residual variance"][1, :Estimate])
    phenotype_corr = ebv_correlation(output, phenotypes)
    mean_pip_all = safe_mean(merged.Model_Frequency)
    mean_pip_causal = safe_mean(merged.Model_Frequency[causal_mask])
    mean_pip_null = safe_mean(merged.Model_Frequency[.!causal_mask])
    mean_pip_enriched = safe_mean(merged.Model_Frequency[enriched_mask])
    mean_pip_non_enriched = safe_mean(merged.Model_Frequency[.!enriched_mask])

    sigmaSq = NaN
    pi_class1 = NaN
    pi_class2 = NaN
    pi_class3 = NaN
    pi_class4 = NaN
    if occursin("BayesR", method_variant)
        sigmaSq = Float64(model.M[1].meanVara)
        if haskey(output, "pi_geno")
            pi_lookup = Dict(
                normalize_pi_class_labels(output["pi_geno"][!, :π]) .=> Float64.(output["pi_geno"][!, :Estimate])
            )
            pi_class1 = get(pi_lookup, "class1", NaN)
            pi_class2 = get(pi_lookup, "class2", NaN)
            pi_class3 = get(pi_lookup, "class3", NaN)
            pi_class4 = get(pi_lookup, "class4", NaN)
        end
    end

    return (
        scenario=String(scenario),
        method=occursin("BayesC", method_variant) ? "BayesC" : "BayesR",
        method_variant=String(method_variant),
        seed=Int(seed),
        phenotype_ebv_correlation=Float64(phenotype_corr),
        residual_variance=Float64(residual_variance),
        mean_pip_all=Float64(mean_pip_all),
        mean_pip_causal=Float64(mean_pip_causal),
        mean_pip_null=Float64(mean_pip_null),
        mean_pip_enriched=Float64(mean_pip_enriched),
        mean_pip_non_enriched=Float64(mean_pip_non_enriched),
        top_causal_recall=Float64(top_causal_recall),
        sigmaSq=Float64(sigmaSq),
        pi_class1=Float64(pi_class1),
        pi_class2=Float64(pi_class2),
        pi_class3=Float64(pi_class3),
        pi_class4=Float64(pi_class4),
        seconds=Float64(seconds),
    ), merged
end

function summarize_runs(runs::DataFrame)
    summary_rows = NamedTuple[]
    metric_names = [
        :phenotype_ebv_correlation,
        :residual_variance,
        :mean_pip_all,
        :mean_pip_causal,
        :mean_pip_null,
        :mean_pip_enriched,
        :mean_pip_non_enriched,
        :top_causal_recall,
        :sigmaSq,
        :pi_class1,
        :pi_class2,
        :pi_class3,
        :pi_class4,
        :seconds,
    ]

    for subdf in groupby(runs, [:scenario, :method_variant])
        values = Dict{Symbol, Float64}()
        sds = Dict{Symbol, Float64}()
        for metric in metric_names
            observed = collect(skipmissing(subdf[!, metric]))
            numeric_values = Float64.(filter(isfinite, observed))
            values[metric] = isempty(numeric_values) ? NaN : mean(numeric_values)
            sds[metric] = isempty(numeric_values) ? NaN : std(numeric_values)
        end
        push!(summary_rows, (
            scenario=String(subdf.scenario[1]),
            method_variant=String(subdf.method_variant[1]),
            method=String(subdf.method[1]),
            n_seeds=nrow(subdf),
            phenotype_ebv_correlation=values[:phenotype_ebv_correlation],
            sd_phenotype_ebv_correlation=sds[:phenotype_ebv_correlation],
            residual_variance=values[:residual_variance],
            sd_residual_variance=sds[:residual_variance],
            mean_pip_all=values[:mean_pip_all],
            sd_mean_pip_all=sds[:mean_pip_all],
            mean_pip_causal=values[:mean_pip_causal],
            sd_mean_pip_causal=sds[:mean_pip_causal],
            mean_pip_null=values[:mean_pip_null],
            sd_mean_pip_null=sds[:mean_pip_null],
            mean_pip_enriched=values[:mean_pip_enriched],
            sd_mean_pip_enriched=sds[:mean_pip_enriched],
            mean_pip_non_enriched=values[:mean_pip_non_enriched],
            sd_mean_pip_non_enriched=sds[:mean_pip_non_enriched],
            top_causal_recall=values[:top_causal_recall],
            sd_top_causal_recall=sds[:top_causal_recall],
            sigmaSq=values[:sigmaSq],
            sd_sigmaSq=sds[:sigmaSq],
            pi_class1=values[:pi_class1],
            sd_pi_class1=sds[:pi_class1],
            pi_class2=values[:pi_class2],
            sd_pi_class2=sds[:pi_class2],
            pi_class3=values[:pi_class3],
            sd_pi_class3=sds[:pi_class3],
            pi_class4=values[:pi_class4],
            sd_pi_class4=sds[:pi_class4],
            seconds=values[:seconds],
            sd_seconds=sds[:seconds],
        ))
    end

    return DataFrame(summary_rows)
end

function run_method_case(case_outdir;
                         dataset,
                         geno_path,
                         pheno_path,
                         method_variant,
                         mcmc_seed,
                         chain_length,
                         burnin,
                         output_samples_frequency,
                         start_h2)
    mkpath(case_outdir)
    phenotypes = CSV.read(pheno_path, DataFrame)
    vary = var(dataset.y)
    start_vare = vary * (1 - start_h2)

    if method_variant == "BayesR"
        start_sigma_sq = vary * start_h2 / (length(dataset.marker_ids) * sum(ANNOT_BENCH_BAYESR_GAMMA .* ANNOT_BENCH_BAYESR_START_PI))
        global geno = get_genotypes(
            geno_path, start_sigma_sq;
            separator=',',
            method="BayesR",
            Pi=copy(ANNOT_BENCH_BAYESR_START_PI),
            estimatePi=true,
            G_is_marker_variance=true,
            estimate_variance=true,
            estimate_scale=false,
            quality_control=false,
            center=false,
        )
    elseif method_variant == "Annotated_BayesR"
        start_sigma_sq = vary * start_h2 / (length(dataset.marker_ids) * sum(ANNOT_BENCH_BAYESR_GAMMA .* ANNOT_BENCH_BAYESR_START_PI))
        global geno = get_genotypes(
            geno_path, start_sigma_sq;
            separator=',',
            method="BayesR",
            Pi=copy(ANNOT_BENCH_BAYESR_START_PI),
            estimatePi=true,
            G_is_marker_variance=true,
            estimate_variance=true,
            estimate_scale=false,
            quality_control=false,
            center=false,
            annotations=dataset.annotations,
        )
    elseif method_variant == "Annotated_BayesC"
        start_genetic_variance = vary * start_h2
        global geno = get_genotypes(
            geno_path, start_genetic_variance;
            separator=',',
            method="BayesC",
            Pi=0.0,
            estimatePi=true,
            estimate_variance=true,
            estimate_scale=false,
            quality_control=false,
            center=false,
            annotations=dataset.annotations,
        )
    else
        error("Unsupported method variant: $method_variant")
    end

    model = build_model("y1 = intercept + geno", start_vare)
    outputEBV(model, geno.obsID)

    output_dir = joinpath(case_outdir, "run")
    isdir(output_dir) && rm(output_dir; recursive=true, force=true)
    elapsed = @elapsed output = runMCMC(
        model,
        phenotypes;
        chain_length=chain_length,
        burnin=burnin,
        output_samples_frequency=output_samples_frequency,
        output_folder=output_dir,
        seed=mcmc_seed,
        outputEBV=true,
        output_heritability=false,
        printout_model_info=false,
        printout_frequency=chain_length + 1,
    )

    return output, model, elapsed, phenotypes
end

function main(argv=ARGS)
    length(argv) > 1 && error("Usage: julia --project=. --startup-file=no benchmarks/annotated_bayesr_comparison.jl [outdir]")

    outdir = length(argv) == 1 ? abspath(argv[1]) : joinpath(@__DIR__, "out", "annotated_bayesr_comparison")
    mkpath(outdir)

    n_obs = env_int("JWAS_ANNOT_BENCH_N_OBS", 200)
    n_markers = env_int("JWAS_ANNOT_BENCH_N_MARKERS", 1000)
    chain_length = env_int("JWAS_ANNOT_BENCH_CHAIN_LENGTH", 10000)
    burnin = env_int("JWAS_ANNOT_BENCH_BURNIN", 2000)
    output_samples_frequency = env_int("JWAS_ANNOT_BENCH_OUTPUT_FREQ", 10)
    seeds = parse_seed_list(get(ENV, "JWAS_ANNOT_BENCH_SEEDS", "2026,2027,2028,2029,2030"))
    data_seed = env_int("JWAS_ANNOT_BENCH_DATA_SEED", 20260328)
    target_h2 = env_float("JWAS_ANNOT_BENCH_H2", 0.45)
    scenario = env_string("JWAS_ANNOT_BENCH_SCENARIO", "sparse_upper_classes")

    dataset = build_annotated_benchmark_dataset(
        seed=data_seed,
        n_obs=n_obs,
        n_markers=n_markers,
        h2=target_h2,
        scenario=scenario,
    )

    datadir = joinpath(outdir, "data")
    write_benchmark_dataset(datadir;
                            ids=dataset.ids,
                            marker_ids=dataset.marker_ids,
                            X=dataset.X,
                            y=dataset.y)
    geno_path = joinpath(datadir, "genotypes.csv")
    pheno_path = joinpath(datadir, "phenotypes.csv")
    CSV.write(joinpath(outdir, "truth_metadata.csv"), dataset.truth_metadata)

    runs = NamedTuple[]
    pip_group_rows = NamedTuple[]
    coefficient_rows = NamedTuple[]

    for method_variant in ("BayesR", "Annotated_BayesR", "Annotated_BayesC")
        for seed in seeds
            case_outdir = joinpath(outdir, String(method_variant), "seed_$(seed)")
            output, model, elapsed, phenotypes = run_method_case(
                case_outdir;
                dataset=dataset,
                geno_path=geno_path,
                pheno_path=pheno_path,
                method_variant=method_variant,
                mcmc_seed=seed,
                chain_length=chain_length,
                burnin=burnin,
                output_samples_frequency=output_samples_frequency,
                start_h2=target_h2,
            )

            run_row, merged = method_case_row(
                method_variant,
                seed,
                dataset.scenario,
                output,
                phenotypes,
                dataset.truth_metadata,
                model,
                elapsed,
            )
            push!(runs, run_row)
            append!(pip_group_rows, marker_summary_rows(method_variant, seed, dataset.scenario, merged))
            append!(coefficient_rows, annotation_coefficient_rows(method_variant, seed, dataset.scenario, output))
        end
    end

    runs_df = sort!(DataFrame(runs), [:method_variant, :seed])
    summary_df = sort!(summarize_runs(runs_df), [:scenario, :method_variant])
    pip_group_df = sort!(DataFrame(pip_group_rows), [:scenario, :method_variant, :seed, :group_type, :group_label])
    coeff_df = DataFrame(coefficient_rows)
    if !isempty(coeff_df)
        sort!(coeff_df, [:scenario, :method_variant, :seed, :Annotation, :Step])
    end

    CSV.write(joinpath(outdir, "comparison_runs.csv"), runs_df)
    CSV.write(joinpath(outdir, "comparison_summary.csv"), summary_df)
    CSV.write(joinpath(outdir, "pip_group_summary.csv"), pip_group_df)
    CSV.write(joinpath(outdir, "annotation_coefficients.csv"), coeff_df)

    println("Wrote annotated BayesR benchmark outputs to ", outdir)
end

main()
