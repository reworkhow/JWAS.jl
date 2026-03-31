using Test
using JWAS
using CSV
using DataFrames
using Distributions
using JWAS.Datasets
using Random

@testset "Annotated BayesR API and validation" begin
    genofile = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
    phenotypes = CSV.read(Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals"), DataFrame, delim=',', missingstring=["NA"])

    function bayesr_annotation_run_error(model, phenotypes; kwargs...)
        outdir = tempname()
        err = try
            runMCMC(
                model,
                phenotypes;
                chain_length=10,
                burnin=0,
                output_samples_frequency=5,
                output_folder=outdir,
                seed=123,
                outputEBV=false,
                output_heritability=false,
                printout_model_info=false,
                kwargs...,
            )
            nothing
        catch exc
            exc
        end
        isdir(outdir) && rm(outdir, recursive=true)
        return err
    end

    function bayesr_annotation_intercepts(default_pi)
        p1 = default_pi[2] + default_pi[3] + default_pi[4]
        p2 = (default_pi[3] + default_pi[4]) / p1
        p3 = default_pi[4] / (default_pi[3] + default_pi[4])
        return (p1 = p1, p2 = p2, p3 = p3,
                b01 = quantile(Normal(), p1),
                b02 = quantile(Normal(), p2),
                b03 = quantile(Normal(), p3))
    end

    @testset "accepts dense single-trait annotations" begin
        annotations = rand(Float64, 5, 2)
        intercepts = bayesr_annotation_intercepts(Float64[0.95, 0.03, 0.015, 0.005])
        geno = get_genotypes(
            genofile, 1.0;
            method="BayesR",
            annotations=annotations,
            separator=',',
            quality_control=false,
        )

        @test geno.method == "BayesR"
        @test geno.annotations !== false
        @test size(geno.annotations.coefficients) == (size(annotations, 2) + 1, 3)
        @test size(geno.annotations.snp_pi) == (geno.nMarkers, 4)
        @test geno.annotations.design_matrix[:, 1] == ones(size(annotations, 1))
        @test isapprox(geno.annotations.coefficients[1, 1], intercepts.b01; atol=1e-12)
        @test isapprox(geno.annotations.coefficients[1, 2], intercepts.b02; atol=1e-12)
        @test isapprox(geno.annotations.coefficients[1, 3], intercepts.b03; atol=1e-12)
    end

    @testset "rejects unsupported BayesR annotation modes" begin
        annotations = rand(Float64, 5, 2)

        begin
            global annotated_bayesr_fast = get_genotypes(
                genofile, 1.0;
                method="BayesR",
                annotations=annotations,
                separator=',',
                quality_control=false,
            )
            model_fast = build_model("y1 = intercept + annotated_bayesr_fast", 1.0)
            err_fast = bayesr_annotation_run_error(model_fast, phenotypes; fast_blocks=true)
            @test err_fast isa Exception
        end

        mktempdir() do tmpdir
            cd(tmpdir) do
                open("annotated_bayesr_stream.csv", "w") do io
                    println(io, "ID,m1,m2,m3,m4,m5")
                    println(io, "a1,0,1,2,1,0")
                    println(io, "a2,1,0,1,2,1")
                    println(io, "a3,2,1,0,1,2")
                    println(io, "a4,0,2,1,0,1")
                end
                prefix = JWAS.prepare_streaming_genotypes(
                    "annotated_bayesr_stream.csv";
                    separator=',',
                    header=true,
                    quality_control=false,
                    center=true,
                )
                err_stream = try
                    global annotated_bayesr_stream = get_genotypes(
                        prefix, 1.0;
                        method="BayesR",
                        storage=:stream,
                        annotations=reshape(rand(Float64, 5, 2), 5, 2),
                    )
                    stream_pheno = DataFrame(
                        ID=["a1", "a2", "a3", "a4"],
                        y1=Float32[1.0, -0.2, 0.5, -0.8],
                    )
                    model_stream = build_model("y1 = intercept + annotated_bayesr_stream", 1.0)
                    bayesr_annotation_run_error(model_stream, stream_pheno; memory_guard=:off)
                catch exc
                    exc
                end
                @test err_stream isa Exception
            end
        end

        begin
            global annotated_bayesr_mt = get_genotypes(
                genofile, [1.0 0.0; 0.0 1.0];
                method="BayesR",
                annotations=annotations,
                separator=',',
                quality_control=false,
            )
            model_mt = build_model(
                "y1 = intercept + annotated_bayesr_mt\ny2 = intercept + annotated_bayesr_mt",
                [1.0 0.0; 0.0 1.0],
            )
            phenotypes_mt = DataFrame(
                ID=copy(phenotypes.ID),
                y1=copy(phenotypes.y1),
                y2=coalesce.(phenotypes.y1, 0.0),
            )
            err_mt = bayesr_annotation_run_error(model_mt, phenotypes_mt)
            @test err_mt isa Exception
        end

        begin
            global annotated_bayesr_rrm = get_genotypes(
                genofile, 1.0;
                method="BayesR",
                annotations=annotations,
                separator=',',
                quality_control=false,
            )
            model_rrm = build_model("y1 = intercept + annotated_bayesr_rrm", 1.0)
            err_rrm = bayesr_annotation_run_error(model_rrm, phenotypes; RRM=ones(Float64, 1, 1))
            @test err_rrm isa Exception
        end
    end

    @testset "rejects degenerate annotated BayesR Pi splits" begin
        annotations = rand(Float64, 5, 2)

        err_p3 = try
            get_genotypes(
                genofile, 1.0;
                method="BayesR",
                Pi=Float64[0.95, 0.03, 0.02, 0.0],
                annotations=annotations,
                separator=',',
                quality_control=false,
            )
            nothing
        catch exc
            exc
        end
        @test err_p3 isa Exception
        @test occursin("classes 3 and 4", sprint(showerror, err_p3))

        err_p2 = try
            get_genotypes(
                genofile, 1.0;
                method="BayesR",
                Pi=Float64[0.95, 0.0, 0.0, 0.05],
                annotations=annotations,
                separator=',',
                quality_control=false,
            )
            nothing
        catch exc
            exc
        end
        @test err_p2 isa Exception
        @test occursin("classes 2 versus 3/4", sprint(showerror, err_p2))
    end

    @testset "default annotation intercept conversion" begin
        default_pi = Float64[0.95, 0.03, 0.015, 0.005]
        intercepts = bayesr_annotation_intercepts(default_pi)

        @test isapprox(intercepts.p1, 0.05; atol=1e-12)
        @test isapprox(intercepts.p2, 0.4; atol=1e-12)
        @test isapprox(intercepts.p3, 0.25; atol=1e-12)
        @test isapprox(intercepts.b01, quantile(Normal(), 0.05); atol=1e-12)
        @test isapprox(intercepts.b02, quantile(Normal(), 0.4); atol=1e-12)
        @test isapprox(intercepts.b03, quantile(Normal(), 0.25); atol=1e-12)
    end

    @testset "initialized BayesR annotation container shape" begin
        annotations = rand(Float64, 5, 2)
        geno = get_genotypes(
            genofile, 1.0;
            method="BayesR",
            annotations=annotations,
            separator=',',
            quality_control=false,
        )

        @test size(geno.annotations.coefficients) == (size(annotations, 2) + 1, 3)
        @test size(geno.annotations.mean_coefficients) == (size(annotations, 2) + 1, 3)
        @test size(geno.annotations.mean_coefficients2) == (size(annotations, 2) + 1, 3)
        @test size(geno.annotations.liability) == (geno.nMarkers, 3)
        @test size(geno.annotations.mu) == (geno.nMarkers, 3)
        @test size(geno.annotations.lower_bound) == (geno.nMarkers, 3)
        @test size(geno.annotations.upper_bound) == (geno.nMarkers, 3)
        @test size(geno.annotations.snp_pi) == (geno.nMarkers, 4)
    end

    @testset "BayesR annotation helper initialization uses four classes" begin
        annotations = rand(Float64, 5, 2)
        geno = get_genotypes(
            genofile, 1.0;
            method="BayesR",
            annotations=annotations,
            separator=',',
            quality_control=false,
        )
        geno.δ = [ones(Int, geno.nMarkers)]

        Random.seed!(2026)
        JWAS.initialize_annotation_indicators!(geno)

        @test length(geno.δ[1]) == geno.nMarkers
        @test all(1 .<= geno.δ[1] .<= 4)
    end

    @testset "binary annotation bounds helper" begin
        lower = fill(999.0, 4)
        upper = fill(999.0, 4)

        JWAS.annotation_binary_bounds!(lower, upper, Int[0, 1, 1, 0])

        @test lower == [-Inf, 0.0, 0.0, -Inf]
        @test upper == [0.0, Inf, Inf, 0.0]
    end

    @testset "BayesR annotation prior refresh builds 4-class snp_pi" begin
        annotations = rand(Float64, 5, 2)
        geno = get_genotypes(
            genofile, 1.0;
            method="BayesR",
            annotations=annotations,
            separator=',',
            quality_control=false,
        )
        geno.δ = [Int[1, 2, 3, 4, 2]]

        JWAS.update_annotation_priors!(geno)

        @test size(geno.annotations.snp_pi) == (geno.nMarkers, 4)
        @test all(isfinite, geno.annotations.snp_pi)
        @test all(abs.(sum(geno.annotations.snp_pi, dims=2) .- 1.0) .< 1e-8)
        @test all((0.0 .< geno.annotations.snp_pi) .& (geno.annotations.snp_pi .< 1.0))
    end

    @testset "dense BayesR sweep accepts per-SNP class priors" begin
        x1 = Float64[0.0, 1.0, 2.0, 1.0]
        x2 = Float64[2.0, 1.0, 0.0, 1.0]
        xArray = [x1, x2]
        xRinvArray = copy(xArray)
        xpRinvx = Float64[sum(abs2, x1), sum(abs2, x2)]
        yCorr = Float64[0.8, -0.1, 0.3, 0.5]
        α = zeros(Float64, 2)
        δ = ones(Int, 2)
        gamma = Float64[0.0, 0.01, 0.1, 1.0]
        snp_pi = Float64[
            0.0 1.0 0.0 0.0
            1.0 0.0 0.0 0.0
        ]

        Random.seed!(20260327)
        JWAS.BayesR!(xArray, xRinvArray, xpRinvx, yCorr, α, δ, 1.0, 0.2, snp_pi, gamma)

        @test δ == Int[2, 1]
    end

    @testset "Annotated BayesR dense run" begin
        annotations = [
            0.0 1.0
            1.0 0.0
            1.0 1.0
            0.0 0.0
            0.5 0.5
        ]

        mktempdir() do tmpdir
            cd(tmpdir) do
                global annotated_bayesr_dense = get_genotypes(
                    genofile, 1.0;
                    separator=',',
                    method="BayesR",
                    quality_control=false,
                    annotations=annotations,
                )
                local model = build_model("y1 = intercept + annotated_bayesr_dense", 1.0)
                local output = runMCMC(
                    model,
                    phenotypes,
                    chain_length=30,
                    burnin=10,
                    output_samples_frequency=10,
                    printout_frequency=31,
                    seed=2026,
                    outputEBV=false,
                    output_heritability=false,
                )

                @test haskey(output, "marker effects annotated_bayesr_dense")
                @test "Model_Frequency" in names(output["marker effects annotated_bayesr_dense"])
                @test haskey(output, "annotation coefficients annotated_bayesr_dense")
                @test names(output["annotation coefficients annotated_bayesr_dense"]) == ["Annotation", "Step", "Estimate", "SD"]
                @test Set(output["annotation coefficients annotated_bayesr_dense"][!, :Step]) ==
                      Set(["step1_zero_vs_nonzero", "step2_small_vs_larger", "step3_medium_vs_large"])
                @test size(model.M[1].annotations.snp_pi, 2) == 4
                @test all(abs.(sum(model.M[1].annotations.snp_pi, dims=2) .- 1.0) .< 1e-8)
            end
        end
    end
end
