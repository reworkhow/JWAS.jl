using CSV
using DataFrames
using Statistics
using JWAS
using JWAS.Datasets

const MT_START_PI = Dict(
    [0.0, 0.0] => 0.96,
    [1.0, 0.0] => 0.015,
    [0.0, 1.0] => 0.015,
    [1.0, 1.0] => 0.01,
)
const ST_BAYESC_PI = 0.98
const ST_BAYESR_PI = Float64[0.99, 0.006, 0.003, 0.001]

env_int(name, default) = parse(Int, get(ENV, name, string(default)))
env_float(name, default) = parse(Float64, get(ENV, name, string(default)))
env_string(name, default) = String(get(ENV, name, default))
env_bool(name, default) = lowercase(env_string(name, string(default))) in ("1", "true", "yes")

function parse_seed_list(raw::AbstractString)
    isempty(strip(raw)) && error("Seed list must not be empty.")
    return parse.(Int, strip.(split(raw, ",")))
end

safe_cor(x, y) = (length(x) == 0 || length(y) == 0 || std(x) == 0 || std(y) == 0) ? NaN : cor(Float64.(x), Float64.(y))

function normalize_marker_id(x)
    s = strip(replace(String(x), "\"" => ""))
    if startswith(s, "m")
        return s
    end
    try
        return "m" * string(parse(Int, s))
    catch
        return s
    end
end

normalize_id(x) = strip(replace(String(x), "\"" => ""))

function top_k_recall(scores::AbstractVector{<:Real}, truth::AbstractVector{Bool}, k::Integer)
    k <= 0 && return NaN
    order = sortperm(Float64.(scores); rev=true)
    selected = order[1:min(k, length(order))]
    return sum(truth[selected]) / k
end

function top_k_count(scores::AbstractVector{<:Real}, truth::AbstractVector{Bool}, k::Integer)
    k <= 0 && return 0
    order = sortperm(Float64.(scores); rev=true)
    selected = order[1:min(k, length(order))]
    return sum(truth[selected])
end

function dataset_paths()
    return (
        genotypes=Datasets.dataset("genotypes.csv", dataset_name="simulated_annotations"),
        phenotypes_mt=Datasets.dataset("phenotypes_mt.csv", dataset_name="simulated_annotations"),
        annotations_mt=Datasets.dataset("annotations_mt.csv", dataset_name="simulated_annotations"),
        truth_mt=Datasets.dataset("truth_mt.csv", dataset_name="simulated_annotations"),
    )
end

function method_cases(include_annotated_bayesr::Bool)
    cases = [
        (variant="MT_BayesC", method="BayesC", annotated=false, multitrait=true, trait=""),
        (variant="MT_Annotated_BayesC", method="BayesC", annotated=true, multitrait=true, trait=""),
        (variant="BayesC_y1", method="BayesC", annotated=false, multitrait=false, trait="y1"),
        (variant="Annotated_BayesC_y1", method="BayesC", annotated=true, multitrait=false, trait="y1"),
        (variant="BayesC_y2", method="BayesC", annotated=false, multitrait=false, trait="y2"),
        (variant="Annotated_BayesC_y2", method="BayesC", annotated=true, multitrait=false, trait="y2"),
        (variant="BayesR_y1", method="BayesR", annotated=false, multitrait=false, trait="y1"),
        (variant="BayesR_y2", method="BayesR", annotated=false, multitrait=false, trait="y2"),
    ]
    if include_annotated_bayesr
        push!(cases, (variant="Annotated_BayesR_y1", method="BayesR", annotated=true, multitrait=false, trait="y1"))
        push!(cases, (variant="Annotated_BayesR_y2", method="BayesR", annotated=true, multitrait=false, trait="y2"))
    end
    return cases
end

function marker_frame(output, geno_name::AbstractString, trait::AbstractString)
    df = copy(output["marker effects " * geno_name])
    sub = filter(:Trait => ==(trait), df)
    rename!(sub, :Estimate => :estimate, :Model_Frequency => :pip)
    sub.estimate = Float64.(sub.estimate)
    sub.pip = Float64.(sub.pip)
    sub.marker_id = normalize_marker_id.(sub.Marker_ID)
    return select(sub, :marker_id, :estimate, :pip)
end

function ebv_frame(output, trait::AbstractString)
    df = copy(output["EBV_" * trait])
    rename!(df, :EBV => :ebv)
    df.ID = normalize_id.(df.ID)
    return select(df, :ID, :ebv)
end

function truth_frame(paths)
    truth = CSV.read(paths.truth_mt, DataFrame)
    truth.marker_id = normalize_marker_id.(truth.marker_id)
    truth.state = string.(truth.state)
    return truth
end

function phenotype_frame(paths)
    pheno = CSV.read(paths.phenotypes_mt, DataFrame)
    pheno.ID = normalize_id.(pheno.ID)
    return pheno
end

function annotation_matrix(paths)
    ann = CSV.read(paths.annotations_mt, DataFrame)
    ann.marker_id = normalize_marker_id.(ann.marker_id)
    return ann
end

function read_joined_marker_table(path::AbstractString)
    df = CSV.read(path, DataFrame)
    :marker_id in names(df) && (df.marker_id = normalize_marker_id.(df.marker_id))
    :estimate in names(df) && (df.estimate = Float64.(df.estimate))
    :pip in names(df) && (df.pip = Float64.(df.pip))
    :pip_y1 in names(df) && (df.pip_y1 = Float64.(df.pip_y1))
    :pip_y2 in names(df) && (df.pip_y2 = Float64.(df.pip_y2))
    :any_score in names(df) && (df.any_score = Float64.(df.any_score))
    :state in names(df) && (df.state = string.(df.state))
    return df
end

function read_marker_sample_matrix(path::AbstractString)
    df = CSV.read(path, DataFrame)
    marker_ids = normalize_marker_id.(string.(names(df)))
    length(unique(marker_ids)) == length(marker_ids) || error("Marker sample file contains duplicate marker IDs: $path")
    sample_matrix = Matrix{Float64}(df)
    return (marker_ids=marker_ids, sample_matrix=sample_matrix, nsamples=nrow(df))
end

function multitrait_joint_state_posterior_table(run_dir::AbstractString, truth::DataFrame)
    y1_path = joinpath(run_dir, "MCMC_samples_marker_effects_bench_geno_y1.txt")
    y2_path = joinpath(run_dir, "MCMC_samples_marker_effects_bench_geno_y2.txt")
    isfile(y1_path) || error("Missing marker-effect sample file: $y1_path")
    isfile(y2_path) || error("Missing marker-effect sample file: $y2_path")

    y1 = read_marker_sample_matrix(y1_path)
    y2 = read_marker_sample_matrix(y2_path)

    y1.nsamples == y2.nsamples || error("Trait sample files do not have the same number of saved samples.")
    y1.marker_ids == y2.marker_ids || error("Trait sample files are not aligned by marker order.")

    on1 = abs.(y1.sample_matrix) .> 0.0
    on2 = abs.(y2.sample_matrix) .> 0.0
    nsamples = y1.nsamples

    posterior = DataFrame(
        marker_id=y1.marker_ids,
        pip_00=vec(sum((.!on1) .& (.!on2); dims=1)) ./ nsamples,
        pip_10=vec(sum(on1 .& (.!on2); dims=1)) ./ nsamples,
        pip_01=vec(sum((.!on1) .& on2; dims=1)) ./ nsamples,
        pip_11=vec(sum(on1 .& on2; dims=1)) ./ nsamples,
        n_saved_samples=fill(nsamples, length(y1.marker_ids)),
    )

    joined = innerjoin(
        posterior,
        select(copy(truth), :marker_id, :state, :is_active_y1, :is_active_y2, :is_shared);
        on=:marker_id,
    )
    assert_marker_join(posterior, select(copy(truth), :marker_id), joined)
    return joined
end

function shared_posterior_row(label::AbstractString, seed::Int, posterior::DataFrame, threshold::Float64)
    shared_truth = Bool.(posterior.is_shared)
    declared = Float64.(posterior.pip_11) .>= threshold
    shared_truth_count = sum(shared_truth)
    declared_shared_count = sum(declared)
    true_shared_declared = sum(declared .& shared_truth)
    meta = pleiotropy_metadata(label)
    return (
        label=label,
        method=meta.method,
        annotated=meta.annotated,
        multitrait=meta.multitrait,
        seed=seed,
        threshold=threshold,
        shared_truth_count=shared_truth_count,
        declared_shared_count=declared_shared_count,
        true_shared_declared_count=true_shared_declared,
        false_shared_declared_count=sum(declared .& .!shared_truth),
        shared_precision=declared_shared_count > 0 ? true_shared_declared / declared_shared_count : NaN,
        shared_recall=shared_truth_count > 0 ? true_shared_declared / shared_truth_count : NaN,
        mean_pip_11_true_shared=shared_truth_count > 0 ? mean(Float64.(posterior.pip_11[shared_truth])) : NaN,
        mean_pip_11_not_shared=sum(.!shared_truth) > 0 ? mean(Float64.(posterior.pip_11[.!shared_truth])) : NaN,
    )
end

function summarize_multitrait_shared_posteriors(output_dir::AbstractString, truth::DataFrame, cases, seeds; threshold::Float64=0.5)
    rows = NamedTuple[]
    for case in cases
        case.multitrait || continue
        for seed in seeds
            run_dir = joinpath(output_dir, case.variant, "seed_$(seed)", "run")
            posterior = multitrait_joint_state_posterior_table(run_dir, truth)
            posterior.declared_shared = Float64.(posterior.pip_11) .>= threshold
            CSV.write(joinpath(output_dir, case.variant, "seed_$(seed)", "posterior_joint_state_probabilities.csv"), posterior)
            push!(rows, shared_posterior_row(case.variant, seed, posterior, threshold))
        end
    end

    per_seed_df = DataFrame(rows)
    summary_df = combine(
        groupby(per_seed_df, [:label, :method, :annotated, :multitrait, :threshold]),
        :shared_truth_count => first => :shared_truth_count,
        :declared_shared_count => mean => :declared_shared_count_mean,
        :declared_shared_count => std => :declared_shared_count_sd,
        :true_shared_declared_count => mean => :true_shared_declared_count_mean,
        :true_shared_declared_count => std => :true_shared_declared_count_sd,
        :false_shared_declared_count => mean => :false_shared_declared_count_mean,
        :false_shared_declared_count => std => :false_shared_declared_count_sd,
        :shared_precision => mean => :shared_precision_mean,
        :shared_precision => std => :shared_precision_sd,
        :shared_recall => mean => :shared_recall_mean,
        :shared_recall => std => :shared_recall_sd,
        :mean_pip_11_true_shared => mean => :mean_pip_11_true_shared_mean,
        :mean_pip_11_true_shared => std => :mean_pip_11_true_shared_sd,
        :mean_pip_11_not_shared => mean => :mean_pip_11_not_shared_mean,
        :mean_pip_11_not_shared => std => :mean_pip_11_not_shared_sd,
    )

    CSV.write(joinpath(output_dir, "multitrait_shared_posterior_per_seed_summary.csv"), per_seed_df)
    CSV.write(joinpath(output_dir, "multitrait_shared_posterior_summary.csv"), summary_df)
    return (per_seed=per_seed_df, summary=summary_df)
end

function pair_table_from_trait_tables(y1_table::DataFrame, y2_table::DataFrame)
    pair = innerjoin(
        select(copy(y1_table), :marker_id, :pip),
        select(copy(y2_table), :marker_id, :pip);
        on=:marker_id,
        renamecols="_y1" => "_y2",
    )
    nrow(y1_table) == nrow(y2_table) || error("Trait marker tables do not have the same row count.")
    length(unique(y1_table.marker_id)) == nrow(y1_table) || error("Trait y1 marker table contains duplicate marker IDs.")
    length(unique(y2_table.marker_id)) == nrow(y2_table) || error("Trait y2 marker table contains duplicate marker IDs.")
    Set(y1_table.marker_id) == Set(y2_table.marker_id) || error("Trait marker ID sets do not match.")
    nrow(pair) == nrow(y1_table) || error("Trait marker join did not produce the expected row count.")
    pair = leftjoin(
        pair,
        select(copy(y1_table), :marker_id, :state, :is_active_y1, :is_active_y2, :is_shared);
        on=:marker_id,
    )
    pair.any_active = Bool.(pair.is_active_y1) .| Bool.(pair.is_active_y2)
    pair.any_score = max.(Float64.(pair.pip_y1), Float64.(pair.pip_y2))
    pair.shared_score = min.(Float64.(pair.pip_y1), Float64.(pair.pip_y2))
    return pair
end

function pleiotropy_metadata(label::AbstractString)
    multitrait = startswith(label, "MT_")
    annotated = occursin("Annotated", label)
    method = occursin("BayesR", label) ? "BayesR" : "BayesC"
    return (method=method, annotated=annotated, multitrait=multitrait)
end

function pleiotropy_row(label::AbstractString, seed::Int, pair_table::DataFrame)
    shared_truth = Bool.(pair_table.is_shared)
    any_truth = Bool.(pair_table.any_active)
    shared_k = sum(shared_truth)
    any_k = sum(any_truth)
    shared_top_count = top_k_count(pair_table.shared_score, shared_truth, shared_k)
    any_selected_count = top_k_count(pair_table.any_score, shared_truth, any_k)
    meta = pleiotropy_metadata(label)
    return (
        label=label,
        method=meta.method,
        annotated=meta.annotated,
        multitrait=meta.multitrait,
        seed=seed,
        shared_truth_count=shared_k,
        any_active_truth_count=any_k,
        shared_topk_count=shared_top_count,
        shared_topk_recall=shared_k > 0 ? shared_top_count / shared_k : NaN,
        shared_count_among_detected_any_active=any_selected_count,
        shared_fraction_among_detected_any_active=any_k > 0 ? any_selected_count / any_k : NaN,
        shared_recall_within_detected_any_active=shared_k > 0 ? any_selected_count / shared_k : NaN,
    )
end

function summarize_pleiotropy_from_family_tables(case_results::Dict, truth::DataFrame, output_root::AbstractString)
    rows = NamedTuple[]

    for (label, seed_map) in case_results
        if startswith(label, "MT_")
            for (seed, trait_map) in seed_map
                haskey(trait_map, "y1") && haskey(trait_map, "y2") || continue
                pair = pair_table_from_trait_tables(trait_map["y1"], trait_map["y2"])
                assert_marker_join(pair, select(copy(truth), :marker_id), pair)
                CSV.write(joinpath(output_root, label, "seed_$(seed)", "pleiotropy_pair_table.csv"), pair)
                push!(rows, pleiotropy_row(label, seed, pair))
            end
        else
            for (seed, trait_map) in seed_map
                haskey(trait_map, "y1") && haskey(trait_map, "y2") || continue
                pair = innerjoin(
                    select(copy(trait_map["y1"]), :marker_id, :pip),
                    select(copy(trait_map["y2"]), :marker_id, :pip);
                    on=:marker_id,
                    renamecols="_y1" => "_y2",
                )
                nrow(trait_map["y1"]) == nrow(trait_map["y2"]) || error("Trait marker tables do not have the same row count.")
                length(unique(trait_map["y1"].marker_id)) == nrow(trait_map["y1"]) || error("Trait y1 marker table contains duplicate marker IDs.")
                length(unique(trait_map["y2"].marker_id)) == nrow(trait_map["y2"]) || error("Trait y2 marker table contains duplicate marker IDs.")
                Set(trait_map["y1"].marker_id) == Set(trait_map["y2"].marker_id) || error("Trait marker ID sets do not match.")
                nrow(pair) == nrow(trait_map["y1"]) || error("Trait marker join did not produce the expected row count.")
                pair = leftjoin(pair, select(copy(truth), :marker_id, :state, :is_active_y1, :is_active_y2, :is_shared); on=:marker_id)
                assert_marker_join(pair, select(copy(truth), :marker_id), pair)
                pair.any_active = Bool.(pair.is_active_y1) .| Bool.(pair.is_active_y2)
                pair.any_score = max.(Float64.(pair.pip_y1), Float64.(pair.pip_y2))
                pair.shared_score = min.(Float64.(pair.pip_y1), Float64.(pair.pip_y2))
                CSV.write(joinpath(output_root, "$(label)_seed_$(seed)_pleiotropy_pair_table.csv"), pair)
                push!(rows, pleiotropy_row(label, seed, pair))
            end
        end
    end

    per_seed_df = DataFrame(rows)
    summary_df = combine(
        groupby(per_seed_df, [:label, :method, :annotated, :multitrait]),
        :shared_truth_count => first => :shared_truth_count,
        :any_active_truth_count => first => :any_active_truth_count,
        :shared_topk_count => mean => :shared_topk_count_mean,
        :shared_topk_count => std => :shared_topk_count_sd,
        :shared_topk_recall => mean => :shared_topk_recall_mean,
        :shared_topk_recall => std => :shared_topk_recall_sd,
        :shared_count_among_detected_any_active => mean => :shared_count_among_detected_any_active_mean,
        :shared_count_among_detected_any_active => std => :shared_count_among_detected_any_active_sd,
        :shared_fraction_among_detected_any_active => mean => :shared_fraction_among_detected_any_active_mean,
        :shared_fraction_among_detected_any_active => std => :shared_fraction_among_detected_any_active_sd,
        :shared_recall_within_detected_any_active => mean => :shared_recall_within_detected_any_active_mean,
        :shared_recall_within_detected_any_active => std => :shared_recall_within_detected_any_active_sd,
    )

    CSV.write(joinpath(output_root, "pleiotropy_per_seed_summary.csv"), per_seed_df)
    CSV.write(joinpath(output_root, "pleiotropy_summary.csv"), summary_df)
    return (per_seed=per_seed_df, summary=summary_df)
end

function load_family_marker_tables(output_dir::AbstractString, cases, seeds)
    family_tables = Dict{String,Dict{Int,Dict{String,DataFrame}}}()
    for case in cases
        family = family_label(case)
        get!(family_tables, family, Dict{Int,Dict{String,DataFrame}}())
        if case.multitrait
            for seed in seeds
                seed_tables = get!(family_tables[family], seed, Dict{String,DataFrame}())
                for trait in ("y1", "y2")
                    path = joinpath(output_dir, case.variant, "seed_$(seed)", "joined_markers_$(trait).csv")
                    isfile(path) || error("Missing joined marker table: $path")
                    seed_tables[trait] = read_joined_marker_table(path)
                end
            end
        else
            for seed in seeds
                path = joinpath(output_dir, case.variant, "seed_$(seed)", "joined_markers_$(case.trait).csv")
                if !isfile(path)
                    continue
                end
                seed_tables = get!(family_tables[family], seed, Dict{String,DataFrame}())
                seed_tables[case.trait] = read_joined_marker_table(path)
            end
        end
    end
    return family_tables
end

function assert_marker_join(marker_df::DataFrame, truth_df::DataFrame, joined_df::DataFrame)
    nrow(marker_df) == nrow(truth_df) || error("Marker output row count does not match truth row count.")
    length(unique(marker_df.marker_id)) == nrow(marker_df) || error("Marker output contains duplicate marker IDs.")
    length(unique(truth_df.marker_id)) == nrow(truth_df) || error("Truth table contains duplicate marker IDs.")
    Set(marker_df.marker_id) == Set(truth_df.marker_id) || error("Marker output and truth marker ID sets do not match.")
    nrow(joined_df) == nrow(truth_df) || error("Marker join did not produce the expected row count.")
end

function assert_id_join(ebv_df::DataFrame, pheno_df::DataFrame, joined_df::DataFrame)
    length(unique(ebv_df.ID)) == nrow(ebv_df) || error("EBV output contains duplicate IDs.")
    length(unique(pheno_df.ID)) == nrow(pheno_df) || error("Phenotype table contains duplicate IDs.")
    Set(ebv_df.ID) == Set(pheno_df.ID) || error("EBV output and phenotype ID sets do not match.")
    nrow(joined_df) == nrow(pheno_df) || error("EBV join did not produce the expected row count.")
end

function annotation_rows(output, case, seed)
    key = "annotation coefficients bench_geno"
    haskey(output, key) || return DataFrame(
        variant=String[],
        method=String[],
        annotated=Bool[],
        multitrait=Bool[],
        trait=String[],
        seed=Int[],
        annotation=String[],
        step=String[],
        estimate=Float64[],
    )
    ann_df = copy(output[key])
    step_values = if "Step" in names(ann_df)
        String.(ann_df.Step)
    else
        fill("step1", nrow(ann_df))
    end
    return DataFrame(
        variant=fill(case.variant, nrow(ann_df)),
        method=fill(case.method, nrow(ann_df)),
        annotated=fill(case.annotated, nrow(ann_df)),
        multitrait=fill(case.multitrait, nrow(ann_df)),
        trait=fill(case.trait, nrow(ann_df)),
        seed=fill(seed, nrow(ann_df)),
        annotation=String.(ann_df.Annotation),
        step=step_values,
        estimate=Float64.(ann_df.Estimate),
    )
end

function run_case(case::NamedTuple, seed::Int, paths;
                  chain_length::Int, burnin::Int, output_samples_frequency::Int,
                  start_h2::Float64, output_root::AbstractString)
    pheno_all = phenotype_frame(paths)
    ann_df = annotation_matrix(paths)
    ann_matrix = Matrix{Float64}(ann_df[:, Not(:marker_id)])

    geno_kwargs = (
        separator=',',
        method=case.method,
        estimatePi=true,
        estimate_scale=false,
        quality_control=false,
        center=false,
    )

    if case.multitrait
        ymat = Matrix{Float64}(pheno_all[:, [:y1, :y2]])
        start_g = Matrix(cov(ymat)) .* start_h2
        start_r = Matrix(cov(ymat)) .* (1 - start_h2)

        global bench_geno = if case.annotated
            get_genotypes(
                paths.genotypes, start_g;
                geno_kwargs...,
                Pi=copy(MT_START_PI),
                annotations=ann_matrix,
            )
        else
            get_genotypes(
                paths.genotypes, start_g;
                geno_kwargs...,
                Pi=copy(MT_START_PI),
            )
        end

        model = build_model("y1 = intercept + bench_geno\ny2 = intercept + bench_geno", start_r)
        outputEBV(model, pheno_all.ID)
    else
        trait = Symbol(case.trait)
        start_g = var(Float64.(pheno_all[!, trait])) * start_h2
        start_r = var(Float64.(pheno_all[!, trait])) * (1 - start_h2)

        if case.method == "BayesC"
            global bench_geno = if case.annotated
                get_genotypes(
                    paths.genotypes, start_g;
                    geno_kwargs...,
                    Pi=ST_BAYESC_PI,
                    annotations=ann_matrix,
                )
            else
                get_genotypes(
                    paths.genotypes, start_g;
                    geno_kwargs...,
                    Pi=ST_BAYESC_PI,
                )
            end
        elseif case.method == "BayesR"
            global bench_geno = if case.annotated
                get_genotypes(
                    paths.genotypes, start_g;
                    geno_kwargs...,
                    Pi=copy(ST_BAYESR_PI),
                    G_is_marker_variance=false,
                    annotations=ann_matrix,
                )
            else
                get_genotypes(
                    paths.genotypes, start_g;
                    geno_kwargs...,
                    Pi=copy(ST_BAYESR_PI),
                    G_is_marker_variance=false,
                )
            end
        else
            error("Unsupported method $(case.method)")
        end

        model = build_model("$(case.trait) = intercept + bench_geno", start_r)
        outputEBV(model, pheno_all.ID)
    end

    seed_dir = joinpath(output_root, case.variant, "seed_$(seed)")
    mkpath(seed_dir)
    raw_output_dir = joinpath(seed_dir, "run")
    isdir(raw_output_dir) && rm(raw_output_dir; recursive=true, force=true)

    elapsed = @elapsed output = runMCMC(
        model,
        pheno_all;
        chain_length=chain_length,
        burnin=burnin,
        output_samples_frequency=output_samples_frequency,
        output_folder=raw_output_dir,
        seed=seed,
        outputEBV=true,
        output_heritability=false,
        printout_model_info=false,
        printout_frequency=chain_length + 1,
    )

    return (output=output, runtime_seconds=elapsed, output_dir=raw_output_dir)
end

function summarize_case(case::NamedTuple, seed::Int, run_result, paths)
    output = run_result.output
    truth = truth_frame(paths)
    pheno = phenotype_frame(paths)

    per_run_rows = NamedTuple[]
    marker_tables = Dict{String,DataFrame}()

    traits = case.multitrait ? ["y1", "y2"] : [case.trait]
    for trait in traits
        marker_df = marker_frame(output, "bench_geno", trait)
        joined_marker = innerjoin(marker_df, truth; on=:marker_id)
        assert_marker_join(marker_df, truth, joined_marker)
        CSV.write(joinpath(dirname(run_result.output_dir), "joined_markers_$(trait).csv"), joined_marker)
        marker_tables[trait] = joined_marker

        active_col = Symbol("is_active_" * trait)
        effect_col = Symbol("true_effect_" * trait)
        active_truth = Bool.(joined_marker[!, active_col])
        k_trait = sum(active_truth)

        ebv_df = ebv_frame(output, trait)
        pheno_trait = select(copy(pheno), :ID, Symbol(trait))
        rename!(pheno_trait, Symbol(trait) => :phenotype)
        joined_ebv = innerjoin(ebv_df, pheno_trait; on=:ID)
        assert_id_join(ebv_df, pheno_trait, joined_ebv)

        push!(per_run_rows, (
            variant=case.variant,
            method=case.method,
            annotated=case.annotated,
            multitrait=case.multitrait,
            trait=trait,
            seed=seed,
            runtime_seconds=run_result.runtime_seconds,
            ebv_correlation=safe_cor(joined_ebv.ebv, joined_ebv.phenotype),
            marker_effect_correlation=safe_cor(joined_marker.estimate, joined_marker[!, effect_col]),
            pip_gap=mean(Float64.(joined_marker.pip[active_truth])) - mean(Float64.(joined_marker.pip[.!active_truth])),
            topk_recall=top_k_recall(joined_marker.pip, active_truth, k_trait),
            active_count=k_trait,
            any_active_topk_recall=NaN,
        ))
    end

    if case.multitrait
        any_join = innerjoin(
            select(copy(marker_tables["y1"]), :marker_id, :pip),
            select(copy(marker_tables["y2"]), :marker_id, :pip);
            on=:marker_id,
            renamecols="_y1" => "_y2",
        )
        any_join = leftjoin(any_join, select(copy(truth), :marker_id, :is_active_y1, :is_active_y2); on=:marker_id)
        any_join.any_active = Bool.(any_join.is_active_y1) .| Bool.(any_join.is_active_y2)
        any_join.any_score = max.(Float64.(any_join.pip_y1), Float64.(any_join.pip_y2))
        any_k = sum(any_join.any_active)
        any_recall = top_k_recall(any_join.any_score, any_join.any_active, any_k)
        CSV.write(joinpath(dirname(run_result.output_dir), "joined_markers_any_active.csv"), any_join)

        per_run_rows = [
            merge(row, (any_active_topk_recall=any_recall,)) for row in per_run_rows
        ]
    end

    return (
        summary=DataFrame(per_run_rows),
        annotations=annotation_rows(output, case, seed),
        marker_tables=marker_tables,
    )
end

function family_label(case::NamedTuple)
    if case.multitrait
        return case.annotated ? "MT_Annotated_BayesC" : "MT_BayesC"
    end
    if case.method == "BayesC"
        return case.annotated ? "Annotated_BayesC_single" : "BayesC_single"
    end
    return case.annotated ? "Annotated_BayesR_single" : "BayesR_single"
end

function summarize_single_family_any_active(case_results::Dict, truth::DataFrame, output_root::AbstractString)
    family_rows = NamedTuple[]
    for (family, seed_map) in case_results
        startswith(family, "MT_") && continue
        for (seed, trait_map) in seed_map
            haskey(trait_map, "y1") && haskey(trait_map, "y2") || continue
            joined = innerjoin(
                select(copy(trait_map["y1"]), :marker_id, :pip),
                select(copy(trait_map["y2"]), :marker_id, :pip);
                on=:marker_id,
                renamecols="_y1" => "_y2",
            )
            nrow(trait_map["y1"]) == nrow(trait_map["y2"]) || error("Trait marker tables do not have the same row count.")
            length(unique(trait_map["y1"].marker_id)) == nrow(trait_map["y1"]) || error("Trait y1 marker table contains duplicate marker IDs.")
            length(unique(trait_map["y2"].marker_id)) == nrow(trait_map["y2"]) || error("Trait y2 marker table contains duplicate marker IDs.")
            Set(trait_map["y1"].marker_id) == Set(trait_map["y2"].marker_id) || error("Trait marker ID sets do not match.")
            nrow(joined) == nrow(trait_map["y1"]) || error("Trait marker join did not produce the expected row count.")
            joined = leftjoin(joined, select(copy(truth), :marker_id, :is_active_y1, :is_active_y2); on=:marker_id)
            assert_marker_join(joined, select(copy(truth), :marker_id), joined)
            joined.any_active = Bool.(joined.is_active_y1) .| Bool.(joined.is_active_y2)
            joined.any_score = max.(Float64.(joined.pip_y1), Float64.(joined.pip_y2))
            any_k = sum(joined.any_active)
            recall = top_k_recall(joined.any_score, joined.any_active, any_k)
            CSV.write(joinpath(output_root, "$(family)_seed_$(seed)_any_active.csv"), joined)
            push!(family_rows, (family=family, seed=seed, any_active_topk_recall=recall))
        end
    end
    return DataFrame(family_rows)
end

function main(args=ARGS)
    output_dir = length(args) == 1 ? abspath(args[1]) :
                 length(args) == 0 ? joinpath(@__DIR__, "out", "simulated_annotations_multitrait_comparison") :
                 error("Usage: julia --project=. --startup-file=no benchmarks/simulated_annotations_multitrait_comparison.jl [OUTPUT_DIR]")
    mkpath(output_dir)

    seeds = parse_seed_list(env_string("JWAS_SIMULATED_MT_SEEDS", "101,202"))
    chain_length = env_int("JWAS_SIMULATED_MT_CHAIN_LENGTH", 2000)
    burnin = env_int("JWAS_SIMULATED_MT_BURNIN", 500)
    output_freq = env_int("JWAS_SIMULATED_MT_OUTPUT_FREQ", 20)
    start_h2 = env_float("JWAS_SIMULATED_MT_START_H2", 0.5)
    include_annotated_bayesr = env_bool("JWAS_SIMULATED_MT_INCLUDE_ANNOTATED_BAYESR", false)
    do_warmup = env_bool("JWAS_SIMULATED_MT_WARMUP", true)
    analyze_only = env_bool("JWAS_SIMULATED_MT_ANALYZE_ONLY", false)

    paths = dataset_paths()
    truth = truth_frame(paths)
    cases = method_cases(include_annotated_bayesr)

    family_marker_tables = Dict{String,Dict{Int,Dict{String,DataFrame}}}()

    if analyze_only
        family_marker_tables = load_family_marker_tables(output_dir, cases, seeds)
    else
        per_run_tables = DataFrame[]
        annotation_tables = DataFrame[]

        for case in cases
            family = family_label(case)
            get!(family_marker_tables, family, Dict{Int,Dict{String,DataFrame}}())
            if do_warmup
                warmup_root = mktempdir()
                try
                    run_case(case, seeds[1], paths;
                             chain_length=min(chain_length, 50),
                             burnin=0,
                             output_samples_frequency=max(1, min(output_freq, 10)),
                             start_h2=start_h2,
                             output_root=warmup_root)
                finally
                    rm(warmup_root; recursive=true, force=true)
                end
            end
            for seed in seeds
                run_result = run_case(case, seed, paths;
                                      chain_length=chain_length,
                                      burnin=burnin,
                                      output_samples_frequency=output_freq,
                                      start_h2=start_h2,
                                      output_root=output_dir)
                summary = summarize_case(case, seed, run_result, paths)
                push!(per_run_tables, summary.summary)
                nrow(summary.annotations) > 0 && push!(annotation_tables, summary.annotations)
                family_marker_tables[family][seed] = get(family_marker_tables[family], seed, Dict{String,DataFrame}())
                for (trait, table) in summary.marker_tables
                    family_marker_tables[family][seed][trait] = table
                end
            end
        end

        per_run_df = vcat(per_run_tables...)
        CSV.write(joinpath(output_dir, "per_run_summary.csv"), per_run_df)

        method_summary = combine(
            groupby(per_run_df, [:variant, :method, :annotated, :multitrait, :trait]),
            :runtime_seconds => mean => :runtime_seconds_mean,
            :runtime_seconds => std => :runtime_seconds_sd,
            :ebv_correlation => mean => :ebv_correlation_mean,
            :ebv_correlation => std => :ebv_correlation_sd,
            :marker_effect_correlation => mean => :marker_effect_correlation_mean,
            :marker_effect_correlation => std => :marker_effect_correlation_sd,
            :pip_gap => mean => :pip_gap_mean,
            :pip_gap => std => :pip_gap_sd,
            :topk_recall => mean => :topk_recall_mean,
            :topk_recall => std => :topk_recall_sd,
            :active_count => first => :active_count,
            :any_active_topk_recall => mean => :any_active_topk_recall_mean,
            :any_active_topk_recall => std => :any_active_topk_recall_sd,
        )
        CSV.write(joinpath(output_dir, "method_summary.csv"), method_summary)

        if !isempty(annotation_tables)
            annotation_df = vcat(annotation_tables...)
            CSV.write(joinpath(output_dir, "annotation_coefficients_by_seed.csv"), annotation_df)
            annotation_summary = combine(
                groupby(annotation_df, [:variant, :method, :annotated, :multitrait, :trait, :annotation, :step]),
                :estimate => mean => :estimate_mean,
                :estimate => std => :estimate_sd,
            )
            CSV.write(joinpath(output_dir, "annotation_coefficients_summary.csv"), annotation_summary)
        end
    end

    family_any_active = summarize_single_family_any_active(family_marker_tables, truth, output_dir)
    CSV.write(joinpath(output_dir, "single_trait_family_any_active_summary.csv"), family_any_active)

    shared_threshold = env_float("JWAS_SIMULATED_MT_SHARED_THRESHOLD", 0.5)
    summarize_multitrait_shared_posteriors(output_dir, truth, cases, seeds; threshold=shared_threshold)

    println("Wrote multi-trait simulated annotations benchmark outputs to ", output_dir)
end

main()
