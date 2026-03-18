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
