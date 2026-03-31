using Test, CSV, DataFrames

include(joinpath(@__DIR__, "..", "..", "benchmarks", "bayesr_parity_common.jl"))

@testset "BayesR parity summary helpers" begin
    outdir = mktempdir()
    fake_output = Dict(
        "pi_geno" => DataFrame(π=["class1", "class2", "class3", "class4"],
                               Estimate=[0.95, 0.03, 0.015, 0.005],
                               SD=zeros(4)),
        "marker effects geno" => DataFrame(Trait=fill("y1", 3),
                                           Marker_ID=["m1", "m2", "m3"],
                                           Estimate=[0.1, -0.2, 0.0],
                                           SD=zeros(3),
                                           Model_Frequency=[1.0, 0.5, 0.0]),
        "residual variance" => DataFrame(Covariance=["y1_y1"], Estimate=[1.25], SD=[0.0])
    )

    write_jwas_parity_summary(fake_output, outdir; sigma_sq=0.8)

    scalar_metrics = CSV.read(joinpath(outdir, "scalar_metrics.csv"), DataFrame)
    pi_summary = CSV.read(joinpath(outdir, "pi.csv"), DataFrame)
    marker_effects = CSV.read(joinpath(outdir, "marker_effects.csv"), DataFrame)

    @test "sigmaSq" in scalar_metrics.metric
    @test nrow(pi_summary) == 4
    @test nrow(marker_effects) == 3
end

@testset "BayesR parity pi labels normalize to class names" begin
    outdir = mktempdir()
    fake_output = Dict(
        "pi_geno" => DataFrame(π=[1.0, 2.0, 3.0, 4.0],
                               Estimate=[0.25, 0.25, 0.25, 0.25],
                               SD=zeros(4)),
        "marker effects geno" => DataFrame(Trait=["y1"], Marker_ID=["m1"], Estimate=[0.0], SD=[0.0], Model_Frequency=[0.0]),
        "residual variance" => DataFrame(Covariance=["y1_y1"], Estimate=[1.0], SD=[0.0])
    )

    write_jwas_parity_summary(fake_output, outdir; sigma_sq=0.8)
    pi_summary = CSV.read(joinpath(outdir, "pi.csv"), DataFrame)

    @test pi_summary.class == ["class1", "class2", "class3", "class4"]
end

@testset "BayesR parity dataset export" begin
    outdir = mktempdir()
    dataset = build_bayesr_parity_dataset(seed=2026, n_obs=12, n_markers=8)
    write_parity_dataset(outdir;
                         ids=dataset.ids,
                         marker_ids=dataset.marker_ids,
                         X=dataset.X,
                         y=dataset.y,
                         gamma=[0.0, 0.01, 0.1, 1.0],
                         start_pi=[0.95, 0.03, 0.015, 0.005],
                         estimate_pi=false,
                         chain_length=40,
                         burnin=10,
                         start_h2=0.5,
                         start_sigma_sq=0.8,
                         start_vare=1.2,
                         seed=2026)
    @test isfile(joinpath(outdir, "genotypes.csv"))
    @test isfile(joinpath(outdir, "phenotypes.csv"))
    @test isfile(joinpath(outdir, "config.csv"))
end

@testset "BayesR parity initial state and trace schema" begin
    outdir = mktempdir()
    dataset = build_bayesr_parity_dataset(seed=2026, n_obs=8, n_markers=5)
    initial_state = build_bayesr_parity_initial_state(dataset.y, length(dataset.marker_ids))

    write_parity_dataset(outdir;
                         ids=dataset.ids,
                         marker_ids=dataset.marker_ids,
                         X=dataset.X,
                         y=dataset.y,
                         gamma=[0.0, 0.01, 0.1, 1.0],
                         start_pi=[0.95, 0.03, 0.015, 0.005],
                         estimate_pi=false,
                         chain_length=100,
                         burnin=0,
                         start_h2=0.5,
                         start_sigma_sq=0.8,
                         start_vare=1.2,
                         seed=2026,
                         initial_state=initial_state)

    @test isfile(joinpath(outdir, "initial_state.csv"))
    @test isfile(joinpath(outdir, "initial_scalars.csv"))

    trace = DataFrame(iter=1:3,
                      sigmaSq=[0.8, 0.7, 0.9],
                      ssq=[0.1, 0.2, 0.3],
                      nnz=[0, 1, 2],
                      vare=[1.2, 1.1, 1.0])
    write_parity_trace(joinpath(outdir, "trace_fixed_pi.csv"), trace)
    loaded = CSV.read(joinpath(outdir, "trace_fixed_pi.csv"), DataFrame)

    @test names(loaded) == ["iter", "sigmaSq", "ssq", "nnz", "vare"]
    @test nrow(loaded) == 3
end

@testset "BayesR parity comparator" begin
    ref_scalars = DataFrame(metric=["sigmaSq", "residual_variance"], value=[0.8, 1.2])
    jwas_scalars = DataFrame(metric=["sigmaSq", "residual_variance"], value=[0.82, 1.19])
    ref_pi = DataFrame(class=["class1", "class2", "class3", "class4"], estimate=[0.94, 0.03, 0.02, 0.01])
    jwas_pi = DataFrame(class=["class1", "class2", "class3", "class4"], estimate=[0.93, 0.04, 0.02, 0.01])
    ref_markers = DataFrame(marker_id=["m1", "m2"], estimate=[0.5, -0.4], model_frequency=[1.0, 0.6])
    jwas_markers = DataFrame(marker_id=["m1", "m2"], estimate=[0.48, -0.39], model_frequency=[0.95, 0.58])

    report = compare_parity_summaries(jwas_scalars, ref_scalars, jwas_pi, ref_pi, jwas_markers, ref_markers)

    @test "sigmaSq" in report.scalar_report.metric
    @test report.marker_correlation > 0.99
end

@testset "BayesR fixed-pi trace comparator" begin
    jwas_trace = DataFrame(iter=1:3, sigmaSq=[1.0, 1.1, 1.2], ssq=[0.4, 0.5, 0.6], nnz=[1, 2, 2], vare=[0.9, 0.8, 0.7])
    ref_trace = DataFrame(iter=1:3, sigmaSq=[1.0, 1.0, 1.3], ssq=[0.4, 0.55, 0.6], nnz=[1, 2, 3], vare=[0.9, 0.85, 0.7])

    report = compare_trace_metrics(jwas_trace, ref_trace)

    @test "sigmaSq_abs_diff" in names(report)
    @test nrow(report) == 3
end

@testset "BayesR generic trace comparator" begin
    dense_trace = DataFrame(
        iter=1:2,
        residual_variance=[0.95, 0.92],
        ycorr_norm=[4.0, 3.8],
        alpha_norm=[0.6, 0.8],
        alpha_abs_mean=[0.05, 0.06],
        nnz=[2, 3],
        max_abs_alpha=[0.3, 0.4],
    )
    block_trace = DataFrame(
        iter=1:2,
        residual_variance=[0.96, 0.88],
        ycorr_norm=[4.0, 3.7],
        alpha_norm=[0.61, 0.79],
        alpha_abs_mean=[0.05, 0.055],
        nnz=[2, 3],
        max_abs_alpha=[0.31, 0.41],
    )

    report = compare_trace_metrics(block_trace, dense_trace)

    @test "residual_variance_abs_diff" in names(report)
    @test "ycorr_norm_abs_diff" in names(report)
    @test "alpha_norm_abs_diff" in names(report)
    @test report.residual_variance_abs_diff[2] ≈ 0.04
end

@testset "BayesR replay draw export and comparison schema" begin
    outdir = mktempdir()
    draws = build_bayesr_replay_draws(5; seed=20260318)
    write_bayesr_replay_draws(joinpath(outdir, "replay_draws_iteration1.csv"), draws)
    loaded = read_bayesr_replay_draws(joinpath(outdir, "replay_draws_iteration1.csv"))

    @test names(loaded) == ["kind", "index", "value"]
    @test count(==("marker_class_uniform"), loaded.kind) == 5
    @test count(==("marker_beta_normal"), loaded.kind) == 5

    marker_cmp = compare_replay_marker_tables(
        DataFrame(marker_id=["m1"], rhs=[1.0], chosen_class=[2], new_alpha=[0.3]),
        DataFrame(marker_id=["m1"], rhs=[1.1], chosen_class=[2], new_alpha=[0.35]),
    )
    @test "rhs_abs_diff" in names(marker_cmp)
end

@testset "BayesR multiseed parity summary helper" begin
    runs = DataFrame(
        seed=[2026, 2030],
        sigma_rel_diff=[0.03, 0.01],
        vare_rel_diff=[0.012, 0.008],
        nonzero_abs_diff=[0.0015, 0.0030],
        max_pi_abs_diff=[0.014, 0.016],
    )

    summary = summarize_multiseed_parity(runs)

    @test "sigma_rel_diff" in summary.metric
    @test summary[summary.metric .== "sigma_rel_diff", :mean_value][1] ≈ 0.02
    @test summary[summary.metric .== "max_pi_abs_diff", :max_value][1] ≈ 0.016
end

@testset "BayesR fixed-hyperparameter summary helper" begin
    outdir = mktempdir()
    fake_output = Dict(
        "marker effects geno" => DataFrame(Trait=fill("y1", 2),
                                           Marker_ID=["m1", "m2"],
                                           Estimate=[0.25, -0.15],
                                           SD=zeros(2),
                                           Model_Frequency=[0.8, 0.3]),
        "residual variance" => DataFrame(Covariance=["y1_y1"], Estimate=[0.9], SD=[0.0])
    )

    write_jwas_parity_summary(fake_output, outdir;
                              sigma_sq=0.75,
                              pi_values=[0.95, 0.03, 0.015, 0.005],
                              fixed_hyperparameters=true)

    scalar_metrics = CSV.read(joinpath(outdir, "scalar_metrics.csv"), DataFrame)
    @test !("sigmaSq" in scalar_metrics.metric)
    @test "residual_variance" in scalar_metrics.metric
    @test "mean_nonzero_frequency" in scalar_metrics.metric
end

@testset "BayesR runtime metadata helper" begin
    runtime_df = bayesr_runtime_report(
        dense_seconds=9.5,
        block_seconds=1.5,
        requested_burnin=700,
        dense_burnin=700,
        block_burnin=100,
        dense_block_size=1,
        block_block_size=7,
    )

    @test names(runtime_df) == ["run", "seconds", "requested_burnin", "effective_burnin", "block_size", "speedup_vs_dense"]
    @test runtime_df[runtime_df.run .== "fast_blocks", :effective_burnin][1] == 100
    @test runtime_df[runtime_df.run .== "fast_blocks", :block_size][1] == 7
end

@testset "BayesR within-method multiseed summary helper" begin
    runs = DataFrame(
        seed=[2026, 2027, 2028],
        residual_variance=[0.91, 0.88, 0.89],
        mean_nonzero_frequency=[0.52, 0.56, 0.54],
        marker_abs_mean=[0.031, 0.028, 0.030],
    )

    summary = summarize_within_method_multiseed(runs)

    @test "residual_variance" in summary.metric
    @test summary[summary.metric .== "mean_nonzero_frequency", :mean_value][1] ≈ mean(runs.mean_nonzero_frequency)
    @test summary[summary.metric .== "marker_abs_mean", :max_value][1] ≈ maximum(runs.marker_abs_mean)
end

@testset "BayesR debug default-blocks single-rep benchmark mode" begin
    outdir = mktempdir()
    repo_root = normpath(joinpath(@__DIR__, "..", ".."))
    benchmark_script = joinpath(repo_root, "benchmarks", "debug", "bayesr_fast_blocks_debug.jl")
    cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) --startup-file=no $benchmark_script $outdir`
    env = copy(ENV)
    env["JWAS_BAYESR_BLOCK_MODE"] = "default_blocks_single_rep"
    env["JWAS_BAYESR_BLOCK_CHAIN_LENGTH"] = "40"
    env["JWAS_BAYESR_BLOCK_BURNIN"] = "10"
    env["JWAS_BAYESR_BLOCK_N_OBS"] = "20"
    env["JWAS_BAYESR_BLOCK_N_MARKERS"] = "12"

    @test success(pipeline(setenv(cmd, env), stdout=devnull, stderr=devnull))
    @test isfile(joinpath(outdir, "comparison_scalar_metrics.csv"))
    @test isfile(joinpath(outdir, "comparison_pi.csv"))
    @test isfile(joinpath(outdir, "comparison_marker_effects_top.csv"))
    @test isfile(joinpath(outdir, "runtime.csv"))
end

@testset "BayesR long-chain schedule comparison benchmark mode" begin
    outdir = mktempdir()
    repo_root = normpath(joinpath(@__DIR__, "..", ".."))
    benchmark_script = joinpath(repo_root, "benchmarks", "bayesr_fast_blocks_parity.jl")
    cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) --startup-file=no $benchmark_script $outdir`
    env = copy(ENV)
    env["JWAS_BAYESR_BLOCK_MODE"] = "long_chain_schedule_comparison"
    env["JWAS_BAYESR_BLOCK_CHAIN_LENGTH"] = "60"
    env["JWAS_BAYESR_BLOCK_BURNIN"] = "20"
    env["JWAS_BAYESR_BLOCK_N_OBS"] = "20"
    env["JWAS_BAYESR_BLOCK_N_MARKERS"] = "12"
    env["JWAS_BAYESR_BLOCK_SEEDS"] = "2026,2027"
    env["JWAS_BAYESR_BLOCK_SETTING"] = "2"

    @test success(pipeline(setenv(cmd, env), stdout=devnull, stderr=devnull))
    @test isfile(joinpath(outdir, "schedule_runs.csv"))
    @test isfile(joinpath(outdir, "schedule_pairwise_summary.csv"))

    runs = CSV.read(joinpath(outdir, "schedule_runs.csv"), DataFrame)
    @test sort!(collect(unique(runs.method))) == ["burnin_gated", "dense"]
    @test all(runs[runs.method .!= "dense", :block_size] .== 2)
end

@testset "BayesC and BayesR fast-block benchmark mode" begin
    outdir = mktempdir()
    repo_root = normpath(joinpath(@__DIR__, "..", ".."))
    benchmark_script = joinpath(repo_root, "benchmarks", "bayesc_bayesr_fast_blocks_comparison.jl")
    cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) --startup-file=no $benchmark_script $outdir`
    env = copy(ENV)
    env["JWAS_METHOD_BLOCK_CHAIN_LENGTH"] = "60"
    env["JWAS_METHOD_BLOCK_BURNIN"] = "20"
    env["JWAS_METHOD_BLOCK_N_OBS"] = "20"
    env["JWAS_METHOD_BLOCK_N_MARKERS"] = "12"
    env["JWAS_METHOD_BLOCK_SEEDS"] = "2026,2027"

    @test success(pipeline(setenv(cmd, env), stdout=devnull, stderr=devnull))
    @test isfile(joinpath(outdir, "comparison_runs.csv"))
    @test isfile(joinpath(outdir, "comparison_pairwise_summary.csv"))

    runs = CSV.read(joinpath(outdir, "comparison_runs.csv"), DataFrame)
    summary = CSV.read(joinpath(outdir, "comparison_pairwise_summary.csv"), DataFrame)
    @test sort!(collect(unique(runs.method_variant))) == [
        "BayesC_dense",
        "BayesC_fast_blocks_1",
        "BayesC_fast_blocks_default",
        "BayesR_dense",
        "BayesR_fast_blocks_1",
        "BayesR_fast_blocks_default",
    ]
    @test "phenotype_ebv_correlation" in names(runs)
    @test "phenotype_ebv_correlation" in summary.metric
end

@testset "Annotated BayesR production benchmark mode" begin
    outdir = mktempdir()
    repo_root = normpath(joinpath(@__DIR__, "..", ".."))
    benchmark_script = joinpath(repo_root, "benchmarks", "annotated_bayesr_comparison.jl")
    cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) --startup-file=no $benchmark_script $outdir`
    env = copy(ENV)
    env["JWAS_ANNOT_BENCH_CHAIN_LENGTH"] = "60"
    env["JWAS_ANNOT_BENCH_BURNIN"] = "20"
    env["JWAS_ANNOT_BENCH_N_OBS"] = "30"
    env["JWAS_ANNOT_BENCH_N_MARKERS"] = "40"
    env["JWAS_ANNOT_BENCH_SEEDS"] = "2026,2027"
    env["JWAS_ANNOT_BENCH_SCENARIO"] = "stepwise_annotation_signal"

    @test success(pipeline(setenv(cmd, env), stdout=devnull, stderr=devnull))
    @test isfile(joinpath(outdir, "comparison_runs.csv"))
    @test isfile(joinpath(outdir, "comparison_summary.csv"))
    @test isfile(joinpath(outdir, "pip_group_summary.csv"))
    @test isfile(joinpath(outdir, "annotation_coefficients.csv"))
    @test isfile(joinpath(outdir, "truth_metadata.csv"))

    runs = CSV.read(joinpath(outdir, "comparison_runs.csv"), DataFrame)
    summary = CSV.read(joinpath(outdir, "comparison_summary.csv"), DataFrame)
    coeffs = CSV.read(joinpath(outdir, "annotation_coefficients.csv"), DataFrame)
    truth = CSV.read(joinpath(outdir, "truth_metadata.csv"), DataFrame)

    @test sort!(collect(unique(runs.method_variant))) == [
        "Annotated_BayesC",
        "Annotated_BayesR",
        "BayesR",
    ]
    @test "phenotype_ebv_correlation" in names(runs)
    @test "scenario" in names(summary)
    @test all(summary.scenario .== "stepwise_annotation_signal")
    @test "mean_pip_causal" in names(summary)
    @test "Step" in names(coeffs)
    @test "scenario" in names(truth)
    @test all(truth.scenario .== "stepwise_annotation_signal")
end
