using Test
using JWAS
using CSV
using DataFrames
using Distributions
using JWAS.Datasets
using Random

@testset "Annotated BayesC API and validation" begin
    genofile = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")

    @testset "rejects unsupported methods" begin
        annotations = rand(Float64, 5, 2)
        err = nothing
        try
            get_genotypes(
                genofile, 1.0;
                method="RR-BLUP",
                annotations=annotations,
                separator=',',
                quality_control=false,
            )
        catch e
            err = e
        end
        @test err isa ErrorException
        @test occursin("annotations", sprint(showerror, err))
    end

    @testset "rejects annotation row mismatches" begin
        annotations = rand(Float64, 4, 2)
        err = nothing
        try
            get_genotypes(
                genofile, 1.0;
                method="BayesC",
                annotations=annotations,
                separator=',',
                quality_control=false,
            )
        catch e
            err = e
        end
        @test err isa ErrorException
        @test occursin("rows", sprint(showerror, err))
    end

    @testset "filters annotations with QC and prepends intercept" begin
        geno_df = DataFrame(
            ID=["a1", "a2", "a3", "a4"],
            m1=[0.0, 1.0, 2.0, 1.0],
            m2=[1.0, 1.0, 1.0, 1.0],
            m3=[2.0, 1.0, 0.0, 1.0],
        )
        annotations = reshape([10.0, 20.0, 30.0], 3, 1)

        geno = get_genotypes(
            geno_df, 1.0;
            method="BayesC",
            annotations=annotations,
            quality_control=true,
            MAF=0.01,
        )

        @test geno.nMarkers == 2
        @test geno.annotations !== false
        @test size(geno.annotations.design_matrix) == (2, 2)
        @test all(geno.annotations.design_matrix[:, 1] .== 1.0)
        @test Matrix(geno.annotations.design_matrix[:, 2:2]) == annotations[[1, 3], :]
    end

    @testset "forces estimatePi for annotations" begin
        annotations = rand(Float64, 5, 2)
        geno = @test_logs (:warn, r"estimatePi=false is ignored when annotations are provided") get_genotypes(
            genofile, 1.0;
            method="BayesC",
            annotations=annotations,
            estimatePi=false,
            separator=',',
            quality_control=false,
        )

        @test geno.estimatePi == true
        @test hasproperty(geno.annotations, :variance)
        @test !hasproperty(geno.annotations, :sd)
    end

    @testset "rejects constant annotation columns" begin
        annotations = [
            1.0 0.0
            1.0 1.0
            1.0 0.0
            1.0 1.0
            1.0 0.5
        ]
        err = nothing
        try
            get_genotypes(
                genofile, 1.0;
                method="BayesC",
                annotations=annotations,
                separator=',',
                quality_control=false,
            )
        catch e
            err = e
        end
        @test err isa ErrorException
        @test occursin("constant", sprint(showerror, err))
    end

    @testset "rejects collinear annotation columns" begin
        annotations = [
            0.0 0.0
            1.0 1.0
            0.0 0.0
            1.0 1.0
            0.5 0.5
        ]
        err = nothing
        try
            get_genotypes(
                genofile, 1.0;
                method="BayesC",
                annotations=annotations,
                separator=',',
                quality_control=false,
            )
        catch e
            err = e
        end
        @test err isa ErrorException
        @test occursin("collinear", sprint(showerror, err))
    end

    @testset "streaming filters raw-marker annotations with backend QC" begin
        mktempdir() do tmpdir
            cd(tmpdir) do
                open("annotated_stream_qc.csv", "w") do io
                    println(io, "ID,m1,m2,m3")
                    println(io, "a1,0,1,2")
                    println(io, "a2,1,1,1")
                    println(io, "a3,2,1,0")
                    println(io, "a4,1,1,2")
                end
                annotations = reshape([10.0, 20.0, 30.0], 3, 1)
                prefix = JWAS.prepare_streaming_genotypes(
                    "annotated_stream_qc.csv";
                    separator=',',
                    header=true,
                    quality_control=true,
                    MAF=0.01,
                    center=true,
                )

                geno = get_genotypes(
                    prefix, 1.0;
                    method="BayesC",
                    storage=:stream,
                    annotations=annotations,
                )

                @test geno.nMarkers == 2
                @test geno.annotations !== false
                @test size(geno.annotations.design_matrix) == (2, 2)
                @test all(geno.annotations.design_matrix[:, 1] .== 1.0)
                @test Matrix(geno.annotations.design_matrix[:, 2:2]) == annotations[[1, 3], :]
            end
        end
    end

    @testset "annotated streaming rejects legacy backends without raw marker mapping" begin
        mktempdir() do tmpdir
            cd(tmpdir) do
                open("annotated_stream_legacy.csv", "w") do io
                    println(io, "ID,m1,m2,m3")
                    println(io, "a1,0,1,2")
                    println(io, "a2,1,1,1")
                    println(io, "a3,2,1,0")
                    println(io, "a4,1,1,2")
                end
                prefix = JWAS.prepare_streaming_genotypes(
                    "annotated_stream_legacy.csv";
                    separator=',',
                    header=true,
                    quality_control=true,
                    MAF=0.01,
                    center=true,
                )
                meta_path = prefix * ".meta"
                legacy_lines = filter(line -> !startswith(line, "selected_path\t") && !startswith(line, "nMarkersAll\t"), readlines(meta_path))
                open(meta_path, "w") do io
                    for line in legacy_lines
                        println(io, line)
                    end
                end

                err = nothing
                try
                    get_genotypes(
                        prefix, 1.0;
                        method="BayesC",
                        storage=:stream,
                        annotations=reshape([10.0, 20.0, 30.0], 3, 1),
                    )
                catch e
                    err = e
                end

                @test err isa ErrorException
                @test occursin("rebuild", sprint(showerror, err))
            end
        end
    end

    @testset "annotated BayesC startup matches Jian-style zero-alpha initialization" begin
        marker_names = Symbol.("m" .* string.(1:20))
        geno_df = DataFrame(ID=["a1", "a2"])
        for (j, marker) in enumerate(marker_names)
            geno_df[!, marker] = Float64[mod(j, 3), mod(j + 1, 3)]
        end
        annotations = reshape(collect(1.0:20.0), 20, 1)

        Random.seed!(2026)
        geno_zero = get_genotypes(
            geno_df, 1.0;
            method="BayesC",
            Pi=0.0,
            annotations=annotations,
            quality_control=false,
        )
        Random.seed!(2027)
        geno_zero_repeat = get_genotypes(
            geno_df, 1.0;
            method="BayesC",
            Pi=0.0,
            annotations=annotations,
            quality_control=false,
        )
        geno_one = get_genotypes(
            geno_df, 1.0;
            method="BayesC",
            Pi=1.0,
            annotations=annotations,
            quality_control=false,
        )
        geno_mid = get_genotypes(
            geno_df, 1.0;
            method="BayesC",
            Pi=0.3,
            annotations=annotations,
            quality_control=false,
        )

        @test geno_zero.π isa AbstractVector
        @test length(geno_zero.π) == geno_zero.nMarkers
        @test all(geno_zero.π .== 0.0)
        @test all(geno_zero.annotations.coefficients .== 0.0)
        @test all(geno_zero.annotations.mu .== 0.0)
        @test geno_zero_repeat.π == geno_zero.π
        @test geno_zero_repeat.annotations.coefficients == geno_zero.annotations.coefficients
        @test geno_zero_repeat.annotations.mu == geno_zero.annotations.mu

        @test geno_one.π isa AbstractVector
        @test length(geno_one.π) == geno_one.nMarkers
        @test all(geno_one.π .== 1.0)
        @test all(geno_one.annotations.coefficients .== 0.0)
        @test all(geno_one.annotations.mu .== 0.0)

        @test geno_mid.π isa AbstractVector
        @test length(geno_mid.π) == geno_mid.nMarkers
        @test all(geno_mid.π .== 0.3)
        @test all(geno_mid.annotations.coefficients .== 0.0)
        @test all(geno_mid.annotations.mu .== 0.0)
    end

    @testset "annotation sampler uses standard probit latent variance" begin
        geno = get_genotypes(
            DataFrame(ID=["a1", "a2"], m1=[0.0, 1.0]),
            1.0;
            method="BayesC",
            quality_control=false,
        )
        geno.δ = [Float64[1.0]]
        geno.π = 0.5
        geno.annotations = JWAS.MarkerAnnotations(reshape([1.0], 1, 1); variance=4.0)
        geno.annotations.mu[1] = 0.3

        Random.seed!(20260309)
        JWAS.update_annotation_priors!(geno)
        actual_liability = copy(geno.annotations.liability)
        actual_coefficients = copy(geno.annotations.coefficients)
        actual_mu = copy(geno.annotations.mu)
        actual_pi = copy(geno.π)

        expected = JWAS.MarkerAnnotations(reshape([1.0], 1, 1); variance=4.0)
        expected.mu[1] = 0.3
        Random.seed!(20260309)
        expected.liability[1] = rand(truncated(Normal(expected.mu[1], 1.0), 0.0, Inf))
        rhs = expected.design_matrix' * expected.liability
        JWAS.Gibbs(expected.lhs, expected.coefficients, rhs, expected.variance)
        expected.mu .= expected.design_matrix * expected.coefficients
        expected_pi = clamp.(1 .- cdf.(Normal(), expected.mu), eps(Float64), 1 - eps(Float64))

        @test actual_liability == expected.liability
        @test actual_coefficients == expected.coefficients
        @test actual_mu == expected.mu
        @test actual_pi isa AbstractVector
        @test actual_pi == expected_pi
    end

    @testset "genetic-to-marker variance setup accepts marker-level BayesC priors" begin
        geno = get_genotypes(
            DataFrame(
                ID=["a1", "a2", "a3"],
                m1=[0.0, 1.0, 2.0],
                m2=[1.0, 1.0, 0.0],
                m3=[2.0, 1.0, 1.0],
            ),
            2.5;
            method="BayesC",
            Pi=0.3,
            annotations=reshape([0.0, 1.0, 0.5], 3, 1),
            quality_control=false,
        )

        JWAS.genetic2marker(geno, geno.π)

        expected = geno.genetic_variance.val / ((1 - 0.3) * geno.sum2pq)
        @test geno.G.val ≈ expected atol=1e-6 rtol=1e-6
    end

    @testset "BayesABC rejects mismatched pi vector length" begin
        xArray = [Float64[1.0, 0.0], Float64[0.0, 1.0]]
        xRinvArray = copy(xArray)
        xpRinvx = Float64[1.0, 1.0]
        yCorr = Float64[0.2, -0.1]
        α = zeros(Float64, 2)
        β = zeros(Float64, 2)
        δ = ones(Float64, 2)
        err = nothing
        try
            JWAS.BayesABC!(xArray, xRinvArray, xpRinvx, yCorr, α, β, δ, 1.0, ones(Float64, 2), Float64[0.5])
        catch e
            err = e
        end
        @test err isa ErrorException
        @test occursin("length", sprint(showerror, err))
        @test occursin("pi", sprint(showerror, err))
    end
end

@testset "Standard BayesC preserves scalar pi output" begin
    phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
    genofile = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
    phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstring=["NA"])

    mktempdir() do tmpdir
        cd(tmpdir) do
            global plain_geno = get_genotypes(
                genofile, 1.0;
                separator=',',
                method="BayesC",
                quality_control=false,
            )
            local model = build_model("y1 = intercept + plain_geno", 1.0)
            local output = runMCMC(
                model,
                phenotypes,
                chain_length=12,
                burnin=2,
                output_samples_frequency=1,
                printout_frequency=99,
                outputEBV=false,
                output_heritability=false,
                seed=2026,
            )

            @test model.M[1].π isa AbstractVector
            @test nrow(output["pi_plain_geno"]) == 1
            @test output["pi_plain_geno"][1, :π] == "π"

            pi_sample_lines = filter(!isempty, readlines(joinpath("results", "MCMC_samples_pi_plain_geno.txt")))
            @test length(pi_sample_lines) == 10
        end
    end
end

@testset "Annotated BayesC dense run" begin
    phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
    genofile = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
    phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstring=["NA"])
    annotations = [
        0.0 1.0
        1.0 0.0
        1.0 1.0
        0.0 0.0
        0.5 0.5
    ]

    mktempdir() do tmpdir
        cd(tmpdir) do
            global annotated_geno = get_genotypes(
                genofile, 1.0;
                separator=',',
                method="BayesC",
                quality_control=false,
                annotations=annotations,
            )
            local model = build_model("y1 = intercept + annotated_geno", 1.0)
            local output = runMCMC(
                model,
                phenotypes,
                chain_length=30,
                burnin=10,
                output_samples_frequency=10,
                printout_frequency=31,
                seed=2026,
            )

            @test model.M[1].π isa AbstractVector
            @test length(model.M[1].π) == model.M[1].nMarkers
            @test haskey(output, "annotation coefficients annotated_geno")
            @test nrow(output["annotation coefficients annotated_geno"]) == size(annotations, 2) + 1

            summary_text = mktemp() do path, io
                redirect_stdout(io) do
                    describe(model)
                end
                flush(io)
                seekstart(io)
                read(io, String)
            end
            @test occursin("π_j (min/mean/max)", summary_text)
            @test !occursin(r"π\s+\[", summary_text)
        end
    end
end

@testset "Annotated BayesC fast_blocks run" begin
    phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
    genofile = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
    phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstring=["NA"])
    annotations = [
        0.0 1.0
        1.0 0.0
        1.0 1.0
        0.0 0.0
        0.5 0.5
    ]

    mktempdir() do tmpdir
        cd(tmpdir) do
            global annotated_blocks = get_genotypes(
                genofile, 1.0;
                separator=',',
                method="BayesC",
                quality_control=false,
                annotations=annotations,
            )
            local model = build_model("y1 = intercept + annotated_blocks", 1.0)
            local output = runMCMC(
                model,
                phenotypes,
                chain_length=30,
                burnin=10,
                output_samples_frequency=10,
                printout_frequency=31,
                seed=2026,
                fast_blocks=true,
            )

            @test model.M[1].π isa AbstractVector
            @test length(model.M[1].π) == model.M[1].nMarkers
            @test haskey(output, "annotation coefficients annotated_blocks")
            @test nrow(output["annotation coefficients annotated_blocks"]) == size(annotations, 2) + 1
        end
    end
end

@testset "Annotated BayesC streaming run" begin
    annotations = [
        0.0 1.0
        1.0 0.0
        1.0 1.0
        0.0 0.0
        0.5 0.5
    ]

    mktempdir() do tmpdir
        cd(tmpdir) do
            open("annotated_stream.csv", "w") do io
                println(io, "ID,m1,m2,m3,m4,m5")
                println(io, "a1,0,1,2,0,1")
                println(io, "a2,1,0,1,2,0")
                println(io, "a3,2,1,0,1,2")
                println(io, "a4,0,2,1,0,1")
                println(io, "a5,1,1,2,2,0")
                println(io, "a6,2,0,0,1,2")
            end
            local prefix = JWAS.prepare_streaming_genotypes(
                "annotated_stream.csv";
                separator=',',
                header=true,
                quality_control=false,
                center=true,
            )
            global annotated_stream = get_genotypes(
                prefix, 1.0;
                method="BayesC",
                storage=:stream,
                annotations=annotations,
            )
            local phenotypes_stream = DataFrame(
                ID=copy(annotated_stream.obsID),
                y1=Float32[1.1, -0.3, 0.8, -0.9, 0.5, -0.1],
            )
            local model = build_model("y1 = intercept + annotated_stream", 1.0)
            local output = runMCMC(
                model,
                phenotypes_stream,
                chain_length=30,
                burnin=10,
                output_samples_frequency=10,
                printout_frequency=31,
                outputEBV=false,
                output_heritability=false,
                seed=2026,
                memory_guard=:off,
            )

            @test model.M[1].π isa AbstractVector
            @test length(model.M[1].π) == model.M[1].nMarkers
            @test haskey(output, "annotation coefficients annotated_stream")
            @test nrow(output["annotation coefficients annotated_stream"]) == size(annotations, 2) + 1
        end
    end
end
