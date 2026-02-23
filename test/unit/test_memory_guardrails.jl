using Test
using JWAS
using DataFrames
using CSV
using JWAS.Datasets

@testset "Memory guardrails: estimator and policy" begin
    @testset "estimate_marker_memory non-block unit weights" begin
        est = JWAS.estimate_marker_memory(10, 20;
                                          element_bytes=4,
                                          has_nonunit_weights=false,
                                          block_starts=false)
        @test est.bytes_X == 10 * 20 * 4
        @test est.bytes_xRinvArray == 0
        @test est.bytes_XRinvArray == 0
        @test est.bytes_XpRinvX == 0
        @test est.bytes_xpRinvx == 20 * 4
        @test est.bytes_total == est.bytes_X + est.bytes_xpRinvx
    end

    @testset "estimate_marker_memory non-block non-unit weights" begin
        est = JWAS.estimate_marker_memory(10, 20;
                                          element_bytes=8,
                                          has_nonunit_weights=true,
                                          block_starts=false)
        @test est.bytes_X == 10 * 20 * 8
        @test est.bytes_xRinvArray == 10 * 20 * 8
        @test est.bytes_XRinvArray == 0
        @test est.bytes_XpRinvX == 0
        @test est.bytes_xpRinvx == 20 * 8
    end

    @testset "estimate_marker_memory block mode extras" begin
        est = JWAS.estimate_marker_memory(10, 20;
                                          element_bytes=4,
                                          has_nonunit_weights=false,
                                          block_starts=[1, 6, 11, 16])
        @test est.bytes_XRinvArray == 10 * 20 * 4
        # block sizes are 5,5,5,5 => sum(s_i^2) = 100
        @test est.bytes_XpRinvX == 100 * 4
    end

    @testset "check_marker_memory_guard! modes" begin
        @test_throws ErrorException JWAS.check_marker_memory_guard!(;
            mode=:error,
            ratio=0.5,
            estimated_bytes=600,
            total_memory_bytes=1000,
            context_string="test")

        @test JWAS.check_marker_memory_guard!(;
            mode=:warn,
            ratio=0.5,
            estimated_bytes=600,
            total_memory_bytes=1000,
            context_string="test") == :warned

        @test JWAS.check_marker_memory_guard!(;
            mode=:off,
            ratio=0.5,
            estimated_bytes=600,
            total_memory_bytes=1000,
            context_string="test") == :skipped
    end

    @testset "check_marker_memory_guard! invalid arguments" begin
        @test_throws ErrorException JWAS.check_marker_memory_guard!(;
            mode=:unknown,
            ratio=0.5,
            estimated_bytes=1,
            total_memory_bytes=1000,
            context_string="test")

        @test_throws ErrorException JWAS.check_marker_memory_guard!(;
            mode=:error,
            ratio=0.0,
            estimated_bytes=1,
            total_memory_bytes=1000,
            context_string="test")
    end
end

@testset "Memory guardrails: runMCMC integration" begin
    phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
    genofile  = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
    phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstring=["NA"])

    @testset "runMCMC throws early in :error mode when guard threshold is tiny" begin
        global geno = get_genotypes(genofile, 1.0, separator=',', method="BayesC")
        model = build_model("y1 = intercept + geno", 1.0)
        @test_throws ErrorException runMCMC(model, phenotypes,
                                            chain_length=10,
                                            burnin=0,
                                            output_samples_frequency=10,
                                            output_folder="guardrail_error_mode",
                                            seed=123,
                                            memory_guard=:error,
                                            memory_guard_ratio=1e-12)
        isdir("guardrail_error_mode") && rm("guardrail_error_mode", recursive=true)
    end

    @testset "runMCMC proceeds in :off mode with same threshold" begin
        global geno = get_genotypes(genofile, 1.0, separator=',', method="BayesC")
        model = build_model("y1 = intercept + geno", 1.0)
        output = runMCMC(model, phenotypes,
                         chain_length=10,
                         burnin=0,
                         output_samples_frequency=10,
                         output_folder="guardrail_off_mode",
                         seed=123,
                         memory_guard=:off,
                         memory_guard_ratio=1e-12)
        @test haskey(output, "location parameters")
        isdir("guardrail_off_mode") && rm("guardrail_off_mode", recursive=true)
    end
end
