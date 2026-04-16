# Unit tests for miscellaneous coverage improvements
# - describe() with data
# - update_priors_frequency
# - outputMCMCsamples for multiple terms
# - Datasets error handling
using Test, JWAS, DataFrames, CSV, JWAS.Datasets

phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
genofile = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstring=["NA"])

@testset "Datasets module" begin
    @testset "Valid dataset access" begin
        f = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
        @test isfile(f)
    end

    @testset "Simulated annotations dataset access" begin
        genofile = Datasets.dataset("genotypes.csv", dataset_name="simulated_annotations")
        phenofile = Datasets.dataset("phenotypes.csv", dataset_name="simulated_annotations")
        annofile = Datasets.dataset("annotations.csv", dataset_name="simulated_annotations")
        truthfile = Datasets.dataset("truth.csv", dataset_name="simulated_annotations")
        phenofile_mt = Datasets.dataset("phenotypes_mt.csv", dataset_name="simulated_annotations")
        annofile_mt = Datasets.dataset("annotations_mt.csv", dataset_name="simulated_annotations")
        truthfile_mt = Datasets.dataset("truth_mt.csv", dataset_name="simulated_annotations")
        rawgenofile = Datasets.dataset("raw_genotypes.txt", dataset_name="simulated_annotations")
        generatorscript = Datasets.dataset("generate_dataset.R", dataset_name="simulated_annotations")

        @test isfile(genofile)
        @test isfile(phenofile)
        @test isfile(annofile)
        @test isfile(truthfile)
        @test isfile(phenofile_mt)
        @test isfile(annofile_mt)
        @test isfile(truthfile_mt)
        @test isfile(rawgenofile)
        @test isfile(generatorscript)

        phen_mt = CSV.read(phenofile_mt, DataFrame)
        ann_mt = CSV.read(annofile_mt, DataFrame)
        truth_mt = CSV.read(truthfile_mt, DataFrame)
        @test names(phen_mt) == ["ID", "y1", "y2"]
        @test names(ann_mt) == ["marker_id", "active_signal", "pleiotropy_signal", "direction_signal", "random_signal"]
        @test names(truth_mt) == ["marker_id", "state", "is_active_y1", "is_active_y2", "is_shared", "true_effect_y1", "true_effect_y2"]
    end

    @testset "Invalid dataset errors" begin
        @test_throws Exception Datasets.dataset("nonexistent.txt", dataset_name="fake")
    end
end

@testset "Reproducibility across methods" begin
    global geno = get_genotypes(genofile, 1.0, separator=',', method="RR-BLUP")
    model = build_model("y1 = intercept + geno", 1.0)

    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    output_folder="test_repro",
                    seed=123)

    @test haskey(output, "location parameters")
    @test haskey(output, "residual variance")
    rm("test_repro", recursive=true)
end

@testset "outputMCMCsamples for multiple terms" begin
    global geno = get_genotypes(genofile, 1.0, separator=',', method="BayesC")
    model = build_model("y1 = intercept + x1 + geno", 1.0)
    set_covariate(model, "x1")
    outputMCMCsamples(model, "intercept", "x1")

    output = runMCMC(model, phenotypes,
                    chain_length=50,
                    output_samples_frequency=10,
                    output_folder="test_multi_samples",
                    seed=123)

    @test isfile("test_multi_samples/MCMC_samples_y1.intercept.txt")
    @test isfile("test_multi_samples/MCMC_samples_y1.x1.txt")
    rm("test_multi_samples", recursive=true)
end

@testset "Covariate in model" begin
    model = build_model("y1 = intercept + x1 + x2", 1.0)
    set_covariate(model, "x1", "x2")
    @test :x1 in model.covVec
    @test :x2 in model.covVec
end

@testset "showMME" begin
    model = build_model("y1 = intercept + x1", 1.0)
    # showMME should not error (requires data to build MME first)
    @test describe(model) === nothing
end

@testset "Multi-trait describe" begin
    R = [1.0 0.5; 0.5 1.0]
    G = [1.0 0.5; 0.5 1.0]
    global geno = get_genotypes(genofile, G, separator=',', method="BayesC")
    model = build_model("y1 = intercept + geno\ny2 = intercept + geno", R)

    output = runMCMC(model, phenotypes,
                    chain_length=50,
                    output_samples_frequency=10,
                    output_folder="test_mt_describe",
                    seed=123)

    # describe() should work after MCMC (MCMCinfo is set)
    @test describe(model) === nothing
    rm("test_mt_describe", recursive=true)
end

@testset "Independent blocks API and explicit block starts" begin
    global ib_geno_exact = get_genotypes(genofile, 1.0, separator=',', method="BayesC")
    exact_model = build_model("y1 = intercept + ib_geno_exact", 1.0)
    exact_output = runMCMC(
        exact_model,
        phenotypes;
        chain_length=6,
        output_samples_frequency=10,
        output_folder="test_independent_blocks_exact",
        seed=123,
        fast_blocks=true,
        independent_blocks=false,
    )
    @test haskey(exact_output, "marker effects ib_geno_exact")
    @test exact_model.MCMCinfo.independent_blocks == false
    isdir("test_independent_blocks_exact") && rm("test_independent_blocks_exact", recursive=true, force=true)

    global ib_geno_independent = get_genotypes(genofile, 1.0, separator=',', method="BayesC")
    independent_model = build_model("y1 = intercept + ib_geno_independent", 1.0)
    independent_output = runMCMC(
        independent_model,
        phenotypes;
        chain_length=6,
        output_samples_frequency=10,
        output_folder="test_independent_blocks_enabled",
        seed=123,
        fast_blocks=true,
        independent_blocks=true,
    )
    @test haskey(independent_output, "marker effects ib_geno_independent")
    @test independent_model.MCMCinfo.independent_blocks == true
    isdir("test_independent_blocks_enabled") && rm("test_independent_blocks_enabled", recursive=true, force=true)

    global ib_geno_requires_blocks = get_genotypes(genofile, 1.0, separator=',', method="BayesC")
    requires_blocks_model = build_model("y1 = intercept + ib_geno_requires_blocks", 1.0)
    @test_throws ErrorException runMCMC(
        requires_blocks_model,
        phenotypes;
        chain_length=6,
        output_folder="test_independent_blocks_requires_blocks",
        independent_blocks=true,
    )
    isdir("test_independent_blocks_requires_blocks") && rm("test_independent_blocks_requires_blocks", recursive=true, force=true)

    global ib_geno_explicit = get_genotypes(genofile, 1.0, separator=',', method="BayesC")
    explicit_model = build_model("y1 = intercept + ib_geno_explicit", 1.0)
    explicit_output = runMCMC(
        explicit_model,
        phenotypes;
        chain_length=6,
        output_samples_frequency=10,
        output_folder="test_independent_blocks_explicit",
        seed=123,
        fast_blocks=[1, 3, 5],
        independent_blocks=true,
    )
    @test haskey(explicit_output, "marker effects ib_geno_explicit")
    @test explicit_model.MCMCinfo.fast_blocks == [1, 3, 5]
    @test explicit_model.MCMCinfo.chain_length == 6
    isdir("test_independent_blocks_explicit") && rm("test_independent_blocks_explicit", recursive=true, force=true)

    for bad_blocks in ([2, 4], [1, 3, 3], [3, 1], [1, 10])
        global ib_geno_bad_blocks = get_genotypes(genofile, 1.0, separator=',', method="BayesC")
        bad_blocks_model = build_model("y1 = intercept + ib_geno_bad_blocks", 1.0)
        @test_throws ErrorException runMCMC(
            bad_blocks_model,
            phenotypes;
            chain_length=6,
            output_folder="test_independent_blocks_bad",
            fast_blocks=bad_blocks,
        )
        isdir("test_independent_blocks_bad") && rm("test_independent_blocks_bad", recursive=true, force=true)
    end
end

@testset "Simulated annotations multitrait benchmark contract" begin
    benchmark_path = abspath(joinpath(@__DIR__, "..", "..", "benchmarks", "simulated_annotations_multitrait_comparison.jl"))
    @test isfile(benchmark_path)

    benchmark_src = read(benchmark_path, String)
    benchmark_mod = Module(:BenchmarkHarnessContract)
    Base.include(benchmark_mod, benchmark_path)

    cases = getfield(benchmark_mod, :method_cases)()
    @test length(cases) == 15
    focused_cases = getfield(benchmark_mod, :selected_method_cases)(:multitrait)
    @test length(focused_cases) == 7
    @test all(case -> case.multitrait, focused_cases)
    focused_plain_empty_cases = getfield(benchmark_mod, :selected_method_cases)(:plain_empty_sampler)
    @test length(focused_plain_empty_cases) == 4
    focused_empty_cases = getfield(benchmark_mod, :selected_method_cases)(:empty_annotated_sampler)
    @test length(focused_empty_cases) == 2

    case_lookup = Dict(case.variant => case for case in cases)
    @test Set(keys(case_lookup)) == Set([
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

    @test case_lookup["MT_BayesC"] == (variant="MT_BayesC", method="BayesC", annotated=false, multitrait=true, trait="", multi_trait_sampler=:auto, annotation_mode=:none)
    @test case_lookup["MT_BayesC_I"] == (variant="MT_BayesC_I", method="BayesC", annotated=false, multitrait=true, trait="", multi_trait_sampler=:I, annotation_mode=:none)
    @test case_lookup["MT_BayesC_II"] == (variant="MT_BayesC_II", method="BayesC", annotated=false, multitrait=true, trait="", multi_trait_sampler=:II, annotation_mode=:none)
    @test case_lookup["MT_Annotated_BayesC_I"] == (variant="MT_Annotated_BayesC_I", method="BayesC", annotated=true, multitrait=true, trait="", multi_trait_sampler=:I, annotation_mode=:real)
    @test case_lookup["MT_Annotated_BayesC_II"] == (variant="MT_Annotated_BayesC_II", method="BayesC", annotated=true, multitrait=true, trait="", multi_trait_sampler=:II, annotation_mode=:real)
    @test case_lookup["MT_EmptyAnnotated_BayesC_I"] == (variant="MT_EmptyAnnotated_BayesC_I", method="BayesC", annotated=true, multitrait=true, trait="", multi_trait_sampler=:I, annotation_mode=:empty)
    @test case_lookup["MT_EmptyAnnotated_BayesC_II"] == (variant="MT_EmptyAnnotated_BayesC_II", method="BayesC", annotated=true, multitrait=true, trait="", multi_trait_sampler=:II, annotation_mode=:empty)
    @test case_lookup["Annotated_BayesR_y1"] == (variant="Annotated_BayesR_y1", method="BayesR", annotated=true, multitrait=false, trait="y1", multi_trait_sampler=:auto, annotation_mode=:real)
    @test case_lookup["Annotated_BayesR_y2"] == (variant="Annotated_BayesR_y2", method="BayesR", annotated=true, multitrait=false, trait="y2", multi_trait_sampler=:auto, annotation_mode=:real)
    @test Set(case.variant for case in focused_cases) == Set([
        "MT_BayesC",
        "MT_BayesC_I",
        "MT_BayesC_II",
        "MT_Annotated_BayesC_I",
        "MT_Annotated_BayesC_II",
        "MT_EmptyAnnotated_BayesC_I",
        "MT_EmptyAnnotated_BayesC_II",
    ])
    @test Set(case.variant for case in focused_plain_empty_cases) == Set([
        "MT_BayesC_I",
        "MT_BayesC_II",
        "MT_EmptyAnnotated_BayesC_I",
        "MT_EmptyAnnotated_BayesC_II",
    ])
    @test Set(case.variant for case in focused_empty_cases) == Set([
        "MT_EmptyAnnotated_BayesC_I",
        "MT_EmptyAnnotated_BayesC_II",
    ])

    non_auto_cases = filter(case -> case.multi_trait_sampler != :auto, cases)
    @test Set(case.variant for case in non_auto_cases) == Set([
        "MT_BayesC_I",
        "MT_BayesC_II",
        "MT_Annotated_BayesC_I",
        "MT_Annotated_BayesC_II",
        "MT_EmptyAnnotated_BayesC_I",
        "MT_EmptyAnnotated_BayesC_II",
    ])

    multi_row = getfield(benchmark_mod, :shared_posterior_row)(
        "MT_Annotated_BayesC_I",
        11,
        DataFrame(
            marker_id=["m1", "m2", "m3"],
            pip_11=[0.2, 0.9, 0.1],
            is_shared=[false, true, false],
        ),
    )
    @test multi_row.shared_score_source == "P11_topk"
    @test multi_row.declared_shared_count == 1
    @test multi_row.true_shared_declared_count == 1
    @test multi_row.shared_precision == 1.0
    @test multi_row.shared_recall == 1.0

    y1_table = DataFrame(
        marker_id=["m1", "m2", "m3"],
        pip=[0.95, 0.10, 0.80],
        state=["10", "01", "11"],
        is_active_y1=[true, false, true],
        is_active_y2=[false, true, true],
        is_shared=[false, false, true],
    )
    y2_table = DataFrame(
        marker_id=["m1", "m2", "m3"],
        pip=[0.10, 0.97, 0.85],
        state=["10", "01", "11"],
        is_active_y1=[true, false, true],
        is_active_y2=[false, true, true],
        is_shared=[false, false, true],
    )
    pair = getfield(benchmark_mod, :pair_table_from_trait_tables)(y1_table, y2_table)
    single_row = getfield(benchmark_mod, :single_trait_shared_overlap_row)("Annotated_BayesC_single", 11, pair)
    @test single_row.shared_score_source == "single_trait_overlap_topk"
    @test single_row.declared_shared_count == 1
    @test single_row.true_shared_declared_count == 1
    @test single_row.shared_precision == 1.0
    @test single_row.shared_recall == 1.0

    multitrait_shared_summary = DataFrame(
        label=["MT_BayesC", "MT_Annotated_BayesC_I"],
        method=["BayesC", "BayesC"],
        annotated=[false, true],
        multitrait=[true, true],
        shared_score_source=["P11_topk", "P11_topk"],
        shared_truth_count=[8, 8],
        declared_shared_count_mean=[4.0, 5.0],
        declared_shared_count_sd=[0.0, 0.5],
        true_shared_declared_count_mean=[3.0, 3.5],
        true_shared_declared_count_sd=[0.0, 0.5],
        false_shared_declared_count_mean=[1.0, 1.5],
        false_shared_declared_count_sd=[0.0, 0.5],
        shared_precision_mean=[0.75, 0.7],
        shared_precision_sd=[0.0, 0.05],
        shared_recall_mean=[0.375, 0.4375],
        shared_recall_sd=[0.0, 0.05],
        shared_f1_mean=[0.50, 0.55],
        shared_f1_sd=[0.0, 0.05],
        mean_pip_11_true_shared_mean=[0.4, 0.45],
        mean_pip_11_true_shared_sd=[0.0, 0.02],
        mean_pip_11_not_shared_mean=[0.01, 0.02],
        mean_pip_11_not_shared_sd=[0.0, 0.01],
    )
    method_summary = DataFrame(
        variant=["MT_BayesC", "MT_Annotated_BayesC_I"],
        method=["BayesC", "BayesC"],
        annotated=[false, true],
        multitrait=[true, true],
        trait=["", ""],
        multi_trait_sampler=["auto", "I"],
        runtime_seconds_mean=[10.0, 12.0],
        runtime_seconds_sd=[1.0, 1.5],
        ebv_correlation_mean=[0.8, 0.81],
        ebv_correlation_sd=[0.01, 0.02],
        marker_effect_correlation_mean=[0.5, 0.55],
        marker_effect_correlation_sd=[0.02, 0.03],
        pip_gap_mean=[0.2, 0.25],
        pip_gap_sd=[0.01, 0.02],
        topk_recall_mean=[0.3, 0.35],
        topk_recall_sd=[0.02, 0.03],
        active_count=[20, 20],
        any_active_topk_recall_mean=[0.4, 0.45],
        any_active_topk_recall_sd=[0.01, 0.02],
    )
    mixing_summary = getfield(benchmark_mod, :summarize_multitrait_mixing)(method_summary, multitrait_shared_summary)
    @test Set(Symbol.(names(mixing_summary))) == Set([
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
    ])
    @test nrow(mixing_summary) == 2
    @test mixing_summary[1, :variant] == "MT_Annotated_BayesC_I" || mixing_summary[2, :variant] == "MT_Annotated_BayesC_I"

    @test occursin("multitrait_shared = summarize_multitrait_shared_posteriors(output_dir, truth, cases, seeds)", benchmark_src)
    @test occursin("single_trait_shared = summarize_single_trait_shared_overlaps(family_marker_tables, truth, output_dir)", benchmark_src)
    @test occursin("mixing_summary = summarize_multitrait_mixing(method_summary, multitrait_shared.summary)", benchmark_src)

    cv_cases = getfield(benchmark_mod, :cv_method_cases)()
    @test Set(case.variant for case in cv_cases) == Set([
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

    ids = ["id4", "id1", "id3", "id2", "id5", "id6"]
    folds_a = getfield(benchmark_mod, :cv_fold_assignments)(ids, 3, 17)
    folds_b = getfield(benchmark_mod, :cv_fold_assignments)(ids, 3, 17)
    @test folds_a == folds_b
    @test sort(folds_a.ID) == sort(ids)
    @test Set(folds_a.fold) == Set(1:3)
    @test all(fold -> sum(folds_a.fold .== fold) > 0, 1:3)
    @test length(unique(folds_a.ID)) == length(ids)

    pheno_df = DataFrame(
        ID=["id1", "id2", "id3"],
        y1=[1.0, 2.0, 3.0],
        y2=[10.0, 20.0, 30.0],
    )
    mt_case = case_lookup["MT_Annotated_BayesC_I"]
    st_case = case_lookup["Annotated_BayesC_y1"]
    heldout_ids = Set(["id2"])

    masked_mt = getfield(benchmark_mod, :masked_phenotype_frame)(pheno_df, mt_case, heldout_ids)
    @test ismissing(masked_mt[masked_mt.ID .== "id2", :y1][1])
    @test ismissing(masked_mt[masked_mt.ID .== "id2", :y2][1])
    @test masked_mt[masked_mt.ID .== "id1", :y1][1] == 1.0
    @test masked_mt[masked_mt.ID .== "id1", :y2][1] == 10.0

    masked_st = getfield(benchmark_mod, :masked_phenotype_frame)(pheno_df, st_case, heldout_ids)
    @test ismissing(masked_st[masked_st.ID .== "id2", :y1][1])
    @test masked_st[masked_st.ID .== "id2", :y2][1] == 20.0

    ebv_df = DataFrame(ID=["id2"], ebv=[2.0])
    heldout_summary = getfield(benchmark_mod, :heldout_prediction_row)(st_case, 101, 2, "y1", ebv_df, pheno_df, heldout_ids, 12.5)
    @test heldout_summary.variant == "Annotated_BayesC_y1"
    @test heldout_summary.seed == 101
    @test heldout_summary.fold == 2
    @test heldout_summary.heldout_count == 1
    @test heldout_summary.runtime_seconds == 12.5
    @test isnan(heldout_summary.heldout_correlation)
    @test heldout_summary.heldout_rmse == 0.0

    cv_per_fold = DataFrame(
        variant=["MT_BayesC", "MT_BayesC", "BayesC_y1", "BayesC_y2"],
        family=["MT_BayesC", "MT_BayesC", "BayesC_single", "BayesC_single"],
        method=["BayesC", "BayesC", "BayesC", "BayesC"],
        annotated=[false, false, false, false],
        multitrait=[true, true, false, false],
        trait=["y1", "y2", "y1", "y2"],
        multi_trait_sampler=["auto", "auto", "auto", "auto"],
        seed=[101, 101, 101, 101],
        fold=[1, 1, 1, 1],
        heldout_count=[10, 10, 10, 10],
        runtime_seconds=[5.0, 5.0, 4.0, 4.0],
        heldout_correlation=[0.8, 0.7, 0.6, 0.5],
        heldout_rmse=[1.0, 1.2, 1.4, 1.6],
    )
    cv_summaries = getfield(benchmark_mod, :summarize_cv_results)(cv_per_fold)
    @test Set(Symbol.(names(cv_summaries.method_summary))) == Set([
        :variant,
        :family,
        :method,
        :annotated,
        :multitrait,
        :trait,
        :multi_trait_sampler,
        :heldout_count_total,
        :runtime_seconds_mean,
        :runtime_seconds_sd,
        :heldout_correlation_mean,
        :heldout_correlation_sd,
        :heldout_rmse_mean,
        :heldout_rmse_sd,
    ])
    @test nrow(cv_summaries.method_summary) == 4
    @test nrow(cv_summaries.family_summary) == 2
    bayesc_single = only(filter(:family => ==("BayesC_single"), cv_summaries.family_summary))
    @test bayesc_single.heldout_correlation_trait_mean == 0.55
    @test bayesc_single.heldout_rmse_trait_mean == 1.5
end
