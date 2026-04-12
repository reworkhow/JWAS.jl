using CSV
using DataFrames
using Random
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
safe_rmse(x, y) = length(x) == 0 || length(y) == 0 ? NaN : sqrt(mean((Float64.(x) .- Float64.(y)) .^ 2))

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

function top_k_mask(scores::AbstractVector{<:Real}, k::Integer)
    mask = falses(length(scores))
    k <= 0 && return mask
    order = sortperm(Float64.(scores); rev=true)
    selected = order[1:min(k, length(order))]
    mask[selected] .= true
    return mask
end

function dataset_paths()
    return (
        genotypes=Datasets.dataset("genotypes.csv", dataset_name="simulated_annotations"),
        phenotypes_mt=Datasets.dataset("phenotypes_mt.csv", dataset_name="simulated_annotations"),
        annotations_mt=Datasets.dataset("annotations_mt.csv", dataset_name="simulated_annotations"),
        truth_mt=Datasets.dataset("truth_mt.csv", dataset_name="simulated_annotations"),
    )
end

function method_cases()
    return [
        (variant="MT_BayesC", method="BayesC", annotated=false, multitrait=true, trait="", multi_trait_sampler=:auto, annotation_mode=:none),
        (variant="MT_BayesC_I", method="BayesC", annotated=false, multitrait=true, trait="", multi_trait_sampler=:I, annotation_mode=:none),
        (variant="MT_BayesC_II", method="BayesC", annotated=false, multitrait=true, trait="", multi_trait_sampler=:II, annotation_mode=:none),
        (variant="MT_Annotated_BayesC_I", method="BayesC", annotated=true, multitrait=true, trait="", multi_trait_sampler=:I, annotation_mode=:real),
        (variant="MT_Annotated_BayesC_II", method="BayesC", annotated=true, multitrait=true, trait="", multi_trait_sampler=:II, annotation_mode=:real),
        (variant="MT_EmptyAnnotated_BayesC_I", method="BayesC", annotated=true, multitrait=true, trait="", multi_trait_sampler=:I, annotation_mode=:empty),
        (variant="MT_EmptyAnnotated_BayesC_II", method="BayesC", annotated=true, multitrait=true, trait="", multi_trait_sampler=:II, annotation_mode=:empty),
        (variant="BayesC_y1", method="BayesC", annotated=false, multitrait=false, trait="y1", multi_trait_sampler=:auto, annotation_mode=:none),
        (variant="Annotated_BayesC_y1", method="BayesC", annotated=true, multitrait=false, trait="y1", multi_trait_sampler=:auto, annotation_mode=:real),
        (variant="BayesC_y2", method="BayesC", annotated=false, multitrait=false, trait="y2", multi_trait_sampler=:auto, annotation_mode=:none),
        (variant="Annotated_BayesC_y2", method="BayesC", annotated=true, multitrait=false, trait="y2", multi_trait_sampler=:auto, annotation_mode=:real),
        (variant="BayesR_y1", method="BayesR", annotated=false, multitrait=false, trait="y1", multi_trait_sampler=:auto, annotation_mode=:none),
        (variant="BayesR_y2", method="BayesR", annotated=false, multitrait=false, trait="y2", multi_trait_sampler=:auto, annotation_mode=:none),
        (variant="Annotated_BayesR_y1", method="BayesR", annotated=true, multitrait=false, trait="y1", multi_trait_sampler=:auto, annotation_mode=:real),
        (variant="Annotated_BayesR_y2", method="BayesR", annotated=true, multitrait=false, trait="y2", multi_trait_sampler=:auto, annotation_mode=:real),
    ]
end

function selected_method_cases(focus_mode::Symbol)
    cases = method_cases()
    if focus_mode == :all
        return cases
    elseif focus_mode == :multitrait
        return filter(case -> case.multitrait, cases)
    elseif focus_mode == :plain_empty_sampler
        selected = Set([
            "MT_BayesC_I",
            "MT_BayesC_II",
            "MT_EmptyAnnotated_BayesC_I",
            "MT_EmptyAnnotated_BayesC_II",
        ])
        return filter(case -> case.variant in selected, cases)
    elseif focus_mode == :empty_annotated_sampler
        selected = Set([
            "MT_EmptyAnnotated_BayesC_I",
            "MT_EmptyAnnotated_BayesC_II",
        ])
        return filter(case -> case.variant in selected, cases)
    else
        error("Unsupported focus mode: $focus_mode")
    end
end

function cv_method_cases()
    selected = Set([
        "MT_BayesC",
        "MT_BayesC_I",
        "MT_BayesC_II",
        "MT_Annotated_BayesC_I",
        "MT_Annotated_BayesC_II",
        "MT_EmptyAnnotated_BayesC_I",
        "MT_EmptyAnnotated_BayesC_II",
        "BayesC_y1",
        "Annotated_BayesC_y1",
        "BayesC_y2",
        "Annotated_BayesC_y2",
        "BayesR_y1",
        "BayesR_y2",
        "Annotated_BayesR_y1",
        "Annotated_BayesR_y2",
    ])
    return filter(case -> case.variant in selected, method_cases())
end

function cv_fold_assignments(ids::AbstractVector, nfolds::Integer, seed::Int)
    nfolds >= 2 || error("Cross-validation requires at least 2 folds.")
    normalized = normalize_id.(ids)
    length(unique(normalized)) == length(normalized) || error("Fold assignment IDs must be unique.")
    length(normalized) >= nfolds || error("Number of folds exceeds the number of IDs.")

    rng = MersenneTwister(seed)
    shuffled = shuffle(rng, sort(normalized))
    folds = [mod(i - 1, nfolds) + 1 for i in eachindex(shuffled)]
    df = DataFrame(ID=shuffled, fold=folds)
    sort!(df, :ID)
    return df
end

function masked_phenotype_frame(pheno_df::DataFrame, case::NamedTuple, heldout_ids)
    heldout_set = Set(normalize_id.(collect(heldout_ids)))
    df = copy(pheno_df)
    df.ID = normalize_id.(df.ID)
    allowmissing!(df)
    heldout_mask = in.(df.ID, Ref(heldout_set))
    if case.multitrait
        df[heldout_mask, :y1] .= missing
        df[heldout_mask, :y2] .= missing
    else
        df[heldout_mask, Symbol(case.trait)] .= missing
    end
    return df
end

function heldout_prediction_join(ebv_df::DataFrame, pheno_df::DataFrame, heldout_ids, trait::AbstractString)
    heldout_set = Set(normalize_id.(collect(heldout_ids)))
    heldout_ebv = filter(:ID => id -> normalize_id(id) in heldout_set, copy(ebv_df))
    heldout_ebv.ID = normalize_id.(heldout_ebv.ID)
    heldout_pheno = filter(:ID => id -> normalize_id(id) in heldout_set, select(copy(pheno_df), :ID, Symbol(trait)))
    heldout_pheno.ID = normalize_id.(heldout_pheno.ID)
    rename!(heldout_pheno, Symbol(trait) => :phenotype)
    joined = innerjoin(heldout_ebv, heldout_pheno; on=:ID)
    assert_id_join(heldout_ebv, heldout_pheno, joined)
    return joined
end

function heldout_prediction_row(case::NamedTuple, seed::Int, fold::Int, trait::AbstractString,
                                ebv_df::DataFrame, pheno_df::DataFrame, heldout_ids, runtime_seconds::Real)
    joined = heldout_prediction_join(ebv_df, pheno_df, heldout_ids, trait)
    return (
        variant=case.variant,
        family=family_label(case),
        method=case.method,
        annotated=case.annotated,
        multitrait=case.multitrait,
        trait=trait,
        multi_trait_sampler=String(case.multi_trait_sampler),
        seed=seed,
        fold=fold,
        heldout_count=nrow(joined),
        runtime_seconds=Float64(runtime_seconds),
        heldout_correlation=safe_cor(joined.ebv, joined.phenotype),
        heldout_rmse=safe_rmse(joined.ebv, joined.phenotype),
    )
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

function annotation_mode_matrix(ann_df::DataFrame, mode::Symbol)
    base = Matrix{Float64}(ann_df[:, Not(:marker_id)])
    if mode == :real
        return base
    elseif mode == :empty
        return zeros(Float64, size(base, 1), 0)
    else
        error("Unsupported annotation mode: $mode")
    end
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

function shared_posterior_row(label::AbstractString, seed::Int, posterior::DataFrame)
    shared_truth = Bool.(posterior.is_shared)
    shared_truth_count = sum(shared_truth)
    declared = top_k_mask(posterior.pip_11, shared_truth_count)
    declared_shared_count = sum(declared)
    true_shared_declared = sum(declared .& shared_truth)
    meta = pleiotropy_metadata(label)
    precision = declared_shared_count > 0 ? true_shared_declared / declared_shared_count : NaN
    recall = shared_truth_count > 0 ? true_shared_declared / shared_truth_count : NaN
    f1 = isnan(precision) || isnan(recall) || (precision + recall == 0) ? NaN : 2 * precision * recall / (precision + recall)
    return (
        label=label,
        method=meta.method,
        annotated=meta.annotated,
        multitrait=meta.multitrait,
        seed=seed,
        shared_score_source="P11_topk",
        shared_truth_count=shared_truth_count,
        declared_shared_count=declared_shared_count,
        true_shared_declared_count=true_shared_declared,
        false_shared_declared_count=sum(declared .& .!shared_truth),
        shared_precision=precision,
        shared_recall=recall,
        shared_f1=f1,
        mean_pip_11_true_shared=shared_truth_count > 0 ? mean(Float64.(posterior.pip_11[shared_truth])) : NaN,
        mean_pip_11_not_shared=sum(.!shared_truth) > 0 ? mean(Float64.(posterior.pip_11[.!shared_truth])) : NaN,
    )
end

function summarize_multitrait_shared_posteriors(output_dir::AbstractString, truth::DataFrame, cases, seeds)
    rows = NamedTuple[]
    for case in cases
        case.multitrait || continue
        for seed in seeds
            run_dir = joinpath(output_dir, case.variant, "seed_$(seed)", "run")
            posterior = multitrait_joint_state_posterior_table(run_dir, truth)
            shared_truth_count = sum(Bool.(posterior.is_shared))
            posterior.declared_shared_topk = top_k_mask(posterior.pip_11, shared_truth_count)
            CSV.write(joinpath(output_dir, case.variant, "seed_$(seed)", "posterior_joint_state_probabilities.csv"), posterior)
            push!(rows, shared_posterior_row(case.variant, seed, posterior))
        end
    end

    per_seed_df = DataFrame(rows)
    if nrow(per_seed_df) == 0
        empty_per_seed = DataFrame(
            label=String[],
            method=String[],
            annotated=Bool[],
            multitrait=Bool[],
            seed=Int[],
            shared_score_source=String[],
            shared_truth_count=Int[],
            declared_shared_count=Int[],
            true_shared_declared_count=Int[],
            false_shared_declared_count=Int[],
            shared_precision=Float64[],
            shared_recall=Float64[],
            shared_f1=Float64[],
            mean_pip_11_true_shared=Float64[],
            mean_pip_11_not_shared=Float64[],
        )
        empty_summary = DataFrame(
            label=String[],
            method=String[],
            annotated=Bool[],
            multitrait=Bool[],
            shared_score_source=String[],
            shared_truth_count=Int[],
            declared_shared_count_mean=Float64[],
            declared_shared_count_sd=Float64[],
            true_shared_declared_count_mean=Float64[],
            true_shared_declared_count_sd=Float64[],
            false_shared_declared_count_mean=Float64[],
            false_shared_declared_count_sd=Float64[],
            shared_precision_mean=Float64[],
            shared_precision_sd=Float64[],
            shared_recall_mean=Float64[],
            shared_recall_sd=Float64[],
            shared_f1_mean=Float64[],
            shared_f1_sd=Float64[],
            mean_pip_11_true_shared_mean=Float64[],
            mean_pip_11_true_shared_sd=Float64[],
            mean_pip_11_not_shared_mean=Float64[],
            mean_pip_11_not_shared_sd=Float64[],
        )
        CSV.write(joinpath(output_dir, "multitrait_shared_posterior_per_seed_summary.csv"), empty_per_seed)
        CSV.write(joinpath(output_dir, "multitrait_shared_posterior_summary.csv"), empty_summary)
        return (per_seed=empty_per_seed, summary=empty_summary)
    end
    summary_df = combine(
        groupby(per_seed_df, [:label, :method, :annotated, :multitrait, :shared_score_source]),
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
        :shared_f1 => mean => :shared_f1_mean,
        :shared_f1 => std => :shared_f1_sd,
        :mean_pip_11_true_shared => mean => :mean_pip_11_true_shared_mean,
        :mean_pip_11_true_shared => std => :mean_pip_11_true_shared_sd,
        :mean_pip_11_not_shared => mean => :mean_pip_11_not_shared_mean,
        :mean_pip_11_not_shared => std => :mean_pip_11_not_shared_sd,
    )

    CSV.write(joinpath(output_dir, "multitrait_shared_posterior_per_seed_summary.csv"), per_seed_df)
    CSV.write(joinpath(output_dir, "multitrait_shared_posterior_summary.csv"), summary_df)
    return (per_seed=per_seed_df, summary=summary_df)
end

function summarize_multitrait_mixing(method_summary::DataFrame, multitrait_shared_summary::DataFrame)
    multitrait_method_summary = filter(:multitrait => identity, method_summary)
    isempty(multitrait_method_summary) && return DataFrame()
    isempty(multitrait_shared_summary) && return DataFrame()

    mt_shared = copy(multitrait_shared_summary)
    rename!(mt_shared, :label => :variant)

    mt_shared = select(
        mt_shared,
        :variant,
        :method,
        :annotated,
        :multitrait,
        :shared_score_source,
        :shared_truth_count,
        :declared_shared_count_mean,
        :declared_shared_count_sd,
        :true_shared_declared_count_mean,
        :true_shared_declared_count_sd,
        :false_shared_declared_count_mean,
        :false_shared_declared_count_sd,
        :shared_precision_mean,
        :shared_precision_sd,
        :shared_recall_mean,
        :shared_recall_sd,
        :shared_f1_mean,
        :shared_f1_sd,
        :mean_pip_11_true_shared_mean,
        :mean_pip_11_true_shared_sd,
        :mean_pip_11_not_shared_mean,
        :mean_pip_11_not_shared_sd,
    )

    mt_method = select(
        multitrait_method_summary,
        :variant,
        :method,
        :annotated,
        :multitrait,
        :multi_trait_sampler,
        :runtime_seconds_mean,
        :runtime_seconds_sd,
        :ebv_correlation_mean,
        :ebv_correlation_sd,
        :marker_effect_correlation_mean,
        :marker_effect_correlation_sd,
        :pip_gap_mean,
        :pip_gap_sd,
        :topk_recall_mean,
        :topk_recall_sd,
        :active_count,
        :any_active_topk_recall_mean,
        :any_active_topk_recall_sd,
    )

    mix = innerjoin(mt_method, mt_shared; on=[:variant, :method, :annotated, :multitrait])
    sort!(mix, [:variant, :multi_trait_sampler])
    return mix
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
    return pair
end

function pleiotropy_metadata(label::AbstractString)
    multitrait = startswith(label, "MT_")
    annotated = occursin("Annotated", label)
    method = occursin("BayesR", label) ? "BayesR" : "BayesC"
    return (method=method, annotated=annotated, multitrait=multitrait)
end

function single_trait_shared_overlap_row(label::AbstractString, seed::Int, pair_table::DataFrame)
    shared_truth = Bool.(pair_table.is_shared)
    k_y1 = sum(Bool.(pair_table.is_active_y1))
    k_y2 = sum(Bool.(pair_table.is_active_y2))
    detected_y1 = top_k_mask(pair_table.pip_y1, k_y1)
    detected_y2 = top_k_mask(pair_table.pip_y2, k_y2)
    declared_shared = detected_y1 .& detected_y2
    shared_truth_count = sum(shared_truth)
    declared_shared_count = sum(declared_shared)
    true_shared_declared = sum(declared_shared .& shared_truth)
    meta = pleiotropy_metadata(label)
    precision = declared_shared_count > 0 ? true_shared_declared / declared_shared_count : NaN
    recall = shared_truth_count > 0 ? true_shared_declared / shared_truth_count : NaN
    f1 = isnan(precision) || isnan(recall) || (precision + recall == 0) ? NaN : 2 * precision * recall / (precision + recall)
    return (
        label=label,
        method=meta.method,
        annotated=meta.annotated,
        multitrait=meta.multitrait,
        seed=seed,
        shared_score_source="single_trait_overlap_topk",
        shared_truth_count=shared_truth_count,
        declared_shared_count=declared_shared_count,
        true_shared_declared_count=true_shared_declared,
        false_shared_declared_count=sum(declared_shared .& .!shared_truth),
        shared_precision=precision,
        shared_recall=recall,
        shared_f1=f1,
        mean_pip_11_true_shared=NaN,
        mean_pip_11_not_shared=NaN,
    )
end

function summarize_single_trait_shared_overlaps(case_results::Dict, truth::DataFrame, output_root::AbstractString)
    rows = NamedTuple[]
    for (family, seed_map) in case_results
        startswith(family, "MT_") && continue
        for (seed, trait_map) in seed_map
            haskey(trait_map, "y1") && haskey(trait_map, "y2") || continue
            pair = pair_table_from_trait_tables(trait_map["y1"], trait_map["y2"])
            assert_marker_join(pair, select(copy(truth), :marker_id), pair)
            k_y1 = sum(Bool.(pair.is_active_y1))
            k_y2 = sum(Bool.(pair.is_active_y2))
            pair.detected_y1_topk = top_k_mask(pair.pip_y1, k_y1)
            pair.detected_y2_topk = top_k_mask(pair.pip_y2, k_y2)
            pair.declared_shared_overlap = pair.detected_y1_topk .& pair.detected_y2_topk
            CSV.write(joinpath(output_root, "$(family)_seed_$(seed)_shared_overlap_table.csv"), pair)
            push!(rows, single_trait_shared_overlap_row(family, seed, pair))
        end
    end

    per_seed_df = DataFrame(rows)
    if nrow(per_seed_df) == 0
        empty_per_seed = DataFrame(
            label=String[],
            method=String[],
            annotated=Bool[],
            multitrait=Bool[],
            seed=Int[],
            shared_score_source=String[],
            shared_truth_count=Int[],
            declared_shared_count=Int[],
            true_shared_declared_count=Int[],
            false_shared_declared_count=Int[],
            shared_precision=Float64[],
            shared_recall=Float64[],
            shared_f1=Float64[],
            mean_pip_11_true_shared=Float64[],
            mean_pip_11_not_shared=Float64[],
        )
        empty_summary = DataFrame(
            label=String[],
            method=String[],
            annotated=Bool[],
            multitrait=Bool[],
            shared_score_source=String[],
            shared_truth_count=Int[],
            declared_shared_count_mean=Float64[],
            declared_shared_count_sd=Float64[],
            true_shared_declared_count_mean=Float64[],
            true_shared_declared_count_sd=Float64[],
            false_shared_declared_count_mean=Float64[],
            false_shared_declared_count_sd=Float64[],
            shared_precision_mean=Float64[],
            shared_precision_sd=Float64[],
            shared_recall_mean=Float64[],
            shared_recall_sd=Float64[],
            shared_f1_mean=Float64[],
            shared_f1_sd=Float64[],
        )
        CSV.write(joinpath(output_root, "single_trait_shared_overlap_per_seed_summary.csv"), empty_per_seed)
        CSV.write(joinpath(output_root, "single_trait_shared_overlap_summary.csv"), empty_summary)
        return (per_seed=empty_per_seed, summary=empty_summary)
    end
    summary_df = combine(
        groupby(per_seed_df, [:label, :method, :annotated, :multitrait, :shared_score_source]),
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
        :shared_f1 => mean => :shared_f1_mean,
        :shared_f1 => std => :shared_f1_sd,
    )

    CSV.write(joinpath(output_root, "single_trait_shared_overlap_per_seed_summary.csv"), per_seed_df)
    CSV.write(joinpath(output_root, "single_trait_shared_overlap_summary.csv"), summary_df)
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
        multi_trait_sampler=String[],
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
        multi_trait_sampler=fill(String(case.multi_trait_sampler), nrow(ann_df)),
        seed=fill(seed, nrow(ann_df)),
        annotation=String.(ann_df.Annotation),
        step=step_values,
        estimate=Float64.(ann_df.Estimate),
    )
end

function run_case(case::NamedTuple, seed::Int, paths;
                  chain_length::Int, burnin::Int, output_samples_frequency::Int,
                  start_h2::Float64, output_root::AbstractString,
                  pheno_input::Union{Nothing,DataFrame}=nothing,
                  fold::Union{Nothing,Int}=nothing)
    pheno_reference = phenotype_frame(paths)
    pheno_run = pheno_input === nothing ? copy(pheno_reference) : copy(pheno_input)
    pheno_run.ID = normalize_id.(pheno_run.ID)
    ann_df = annotation_matrix(paths)
    ann_matrix = case.annotated ? annotation_mode_matrix(ann_df, case.annotation_mode) : zeros(0, 0)

    geno_kwargs = (
        separator=',',
        method=case.method,
        estimatePi=true,
        estimate_scale=false,
        quality_control=false,
        center=false,
    )

    if case.multitrait
        ymat = Matrix{Float64}(pheno_reference[:, [:y1, :y2]])
        start_g = Matrix(cov(ymat)) .* start_h2
        start_r = Matrix(cov(ymat)) .* (1 - start_h2)

        global bench_geno = if case.annotated
            get_genotypes(
                paths.genotypes, start_g;
                geno_kwargs...,
                multi_trait_sampler=case.multi_trait_sampler,
                Pi=copy(MT_START_PI),
                annotations=ann_matrix,
            )
        else
            get_genotypes(
                paths.genotypes, start_g;
                geno_kwargs...,
                multi_trait_sampler=case.multi_trait_sampler,
                Pi=copy(MT_START_PI),
            )
        end

        model = build_model("y1 = intercept + bench_geno\ny2 = intercept + bench_geno", start_r)
        outputEBV(model, pheno_run.ID)
    else
        trait = Symbol(case.trait)
        start_g = var(Float64.(pheno_reference[!, trait])) * start_h2
        start_r = var(Float64.(pheno_reference[!, trait])) * (1 - start_h2)

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
        outputEBV(model, pheno_run.ID)
    end

    seed_dir = fold === nothing ?
        joinpath(output_root, case.variant, "seed_$(seed)") :
        joinpath(output_root, case.variant, "seed_$(seed)", "fold_$(fold)")
    mkpath(seed_dir)
    raw_output_dir = joinpath(seed_dir, "run")
    isdir(raw_output_dir) && rm(raw_output_dir; recursive=true, force=true)

    elapsed = @elapsed output = runMCMC(
        model,
        pheno_run;
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
            multi_trait_sampler=String(case.multi_trait_sampler),
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

function summarize_case_cv(case::NamedTuple, seed::Int, fold::Int, run_result, pheno_reference::DataFrame, heldout_ids)
    output = run_result.output
    traits = case.multitrait ? ["y1", "y2"] : [case.trait]
    rows = NamedTuple[]
    for trait in traits
        ebv_df = ebv_frame(output, trait)
        joined = heldout_prediction_join(ebv_df, pheno_reference, heldout_ids, trait)
        CSV.write(joinpath(dirname(run_result.output_dir), "heldout_ebv_$(trait).csv"), joined)
        push!(rows, heldout_prediction_row(case, seed, fold, trait, ebv_df, pheno_reference, heldout_ids, run_result.runtime_seconds))
    end
    return DataFrame(rows)
end

function family_label(case::NamedTuple)
    if case.multitrait
        return case.variant
    end
    if case.method == "BayesC"
        return case.annotated ? "Annotated_BayesC_single" : "BayesC_single"
    end
    return case.annotated ? "Annotated_BayesR_single" : "BayesR_single"
end

function summarize_cv_results(per_fold_df::DataFrame)
    method_summary = combine(
        groupby(per_fold_df, [:variant, :family, :method, :annotated, :multitrait, :trait, :multi_trait_sampler]),
        :heldout_count => sum => :heldout_count_total,
        :runtime_seconds => mean => :runtime_seconds_mean,
        :runtime_seconds => std => :runtime_seconds_sd,
        :heldout_correlation => mean => :heldout_correlation_mean,
        :heldout_correlation => std => :heldout_correlation_sd,
        :heldout_rmse => mean => :heldout_rmse_mean,
        :heldout_rmse => std => :heldout_rmse_sd,
    )

    family_summary = combine(
        groupby(method_summary, [:family, :method, :annotated, :multitrait]),
        :heldout_count_total => sum => :heldout_count_total,
        :runtime_seconds_mean => mean => :runtime_seconds_mean,
        :runtime_seconds_mean => std => :runtime_seconds_sd,
        :heldout_correlation_mean => mean => :heldout_correlation_trait_mean,
        :heldout_correlation_mean => std => :heldout_correlation_trait_sd,
        :heldout_rmse_mean => mean => :heldout_rmse_trait_mean,
        :heldout_rmse_mean => std => :heldout_rmse_trait_sd,
    )

    sort!(method_summary, [:variant, :trait])
    sort!(family_summary, [:family])
    return (method_summary=method_summary, family_summary=family_summary)
end

function run_cv_benchmark(output_dir::AbstractString, paths, cases, seeds;
                          chain_length::Int, burnin::Int, output_freq::Int,
                          start_h2::Float64, do_warmup::Bool, nfolds::Int)
    pheno_reference = phenotype_frame(paths)
    assignment_tables = DataFrame[]
    fold_maps = Dict{Int,DataFrame}()

    for seed in seeds
        fold_df = cv_fold_assignments(pheno_reference.ID, nfolds, seed)
        fold_maps[seed] = fold_df
        push!(assignment_tables, DataFrame(seed=fill(seed, nrow(fold_df)), ID=fold_df.ID, fold=fold_df.fold))
    end

    fold_assignments = vcat(assignment_tables...)
    CSV.write(joinpath(output_dir, "cv_fold_assignments.csv"), fold_assignments)

    per_fold_tables = DataFrame[]

    for case in cases
        if do_warmup
            warmup_root = mktempdir()
            try
                warmup_seed = first(seeds)
                warmup_fold = 1
                heldout_ids = Set(fold_maps[warmup_seed].ID[fold_maps[warmup_seed].fold .== warmup_fold])
                warmup_pheno = masked_phenotype_frame(pheno_reference, case, heldout_ids)
                run_case(case, warmup_seed, paths;
                         chain_length=min(chain_length, 50),
                         burnin=0,
                         output_samples_frequency=max(1, min(output_freq, 10)),
                         start_h2=start_h2,
                         output_root=warmup_root,
                         pheno_input=warmup_pheno,
                         fold=warmup_fold)
            finally
                rm(warmup_root; recursive=true, force=true)
            end
        end

        for seed in seeds
            fold_df = fold_maps[seed]
            for fold in 1:nfolds
                heldout_ids = Set(fold_df.ID[fold_df.fold .== fold])
                masked_pheno = masked_phenotype_frame(pheno_reference, case, heldout_ids)
                run_result = run_case(case, seed, paths;
                                      chain_length=chain_length,
                                      burnin=burnin,
                                      output_samples_frequency=output_freq,
                                      start_h2=start_h2,
                                      output_root=output_dir,
                                      pheno_input=masked_pheno,
                                      fold=fold)
                push!(per_fold_tables, summarize_case_cv(case, seed, fold, run_result, pheno_reference, heldout_ids))
            end
        end
    end

    per_fold_df = vcat(per_fold_tables...)
    CSV.write(joinpath(output_dir, "cv_per_fold_summary.csv"), per_fold_df)
    cv_summaries = summarize_cv_results(per_fold_df)
    CSV.write(joinpath(output_dir, "cv_method_summary.csv"), cv_summaries.method_summary)
    CSV.write(joinpath(output_dir, "cv_family_summary.csv"), cv_summaries.family_summary)
    return (fold_assignments=fold_assignments, per_fold=per_fold_df,
            method_summary=cv_summaries.method_summary, family_summary=cv_summaries.family_summary)
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
    do_warmup = env_bool("JWAS_SIMULATED_MT_WARMUP", true)
    analyze_only = env_bool("JWAS_SIMULATED_MT_ANALYZE_ONLY", false)
    benchmark_mode = Symbol(env_string("JWAS_SIMULATED_MT_BENCHMARK_MODE", "standard"))
    cv_folds = env_int("JWAS_SIMULATED_MT_CV_FOLDS", 5)
    focused_multitrait_only = env_bool("JWAS_SIMULATED_MT_FOCUSED_MULTITRAIT", false)
    focus_mode = Symbol(env_string("JWAS_SIMULATED_MT_FOCUS_MODE", focused_multitrait_only ? "multitrait" : "all"))

    paths = dataset_paths()
    truth = truth_frame(paths)
    cases = benchmark_mode == :cv ? cv_method_cases() : selected_method_cases(focus_mode)

    if benchmark_mode == :cv
        analyze_only && error("Analyze-only mode is not supported for the CV benchmark.")
        run_cv_benchmark(output_dir, paths, cases, seeds;
                         chain_length=chain_length,
                         burnin=burnin,
                         output_freq=output_freq,
                         start_h2=start_h2,
                         do_warmup=do_warmup,
                         nfolds=cv_folds)
        println("Wrote simulated annotations CV benchmark outputs to ", output_dir)
        return
    end

    family_marker_tables = Dict{String,Dict{Int,Dict{String,DataFrame}}}()
    method_summary = DataFrame()
    multitrait_shared = nothing

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
            groupby(per_run_df, [:variant, :method, :annotated, :multitrait, :trait, :multi_trait_sampler]),
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
                groupby(annotation_df, [:variant, :method, :annotated, :multitrait, :trait, :multi_trait_sampler, :annotation, :step]),
                :estimate => mean => :estimate_mean,
                :estimate => std => :estimate_sd,
            )
            CSV.write(joinpath(output_dir, "annotation_coefficients_summary.csv"), annotation_summary)
        end

        multitrait_shared = summarize_multitrait_shared_posteriors(output_dir, truth, cases, seeds)
        mixing_summary = summarize_multitrait_mixing(method_summary, multitrait_shared.summary)
        if nrow(mixing_summary) > 0
            CSV.write(joinpath(output_dir, "multitrait_mixing_summary.csv"), mixing_summary)
        end

    end

    family_any_active = summarize_single_family_any_active(family_marker_tables, truth, output_dir)
    CSV.write(joinpath(output_dir, "single_trait_family_any_active_summary.csv"), family_any_active)

    multitrait_shared === nothing && (multitrait_shared = summarize_multitrait_shared_posteriors(output_dir, truth, cases, seeds))
    single_trait_shared = summarize_single_trait_shared_overlaps(family_marker_tables, truth, output_dir)
    CSV.write(joinpath(output_dir, "pleiotropy_per_seed_summary.csv"), vcat(multitrait_shared.per_seed, single_trait_shared.per_seed; cols=:union))
    CSV.write(joinpath(output_dir, "pleiotropy_summary.csv"), vcat(multitrait_shared.summary, single_trait_shared.summary; cols=:union))

    println("Wrote multi-trait simulated annotations benchmark outputs to ", output_dir)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
