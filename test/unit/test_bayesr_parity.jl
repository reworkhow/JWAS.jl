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
