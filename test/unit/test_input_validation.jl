# Unit tests for input data validation
using Test, JWAS, DataFrames, CSV, JWAS.Datasets

@testset "Input Validation" begin
    phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
    genofile = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
    phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstrings=["NA"])

    @testset "Invalid Bayesian method" begin
        global geno = get_genotypes(genofile, 1.0, separator=',', method="BayesC")
        # Manually set invalid method to trigger validation
        geno.method = "InvalidMethod"
        model = build_model("y1 = intercept + geno", 1.0)
        @test_throws ErrorException runMCMC(model, phenotypes,
                                           chain_length=10,
                                           output_folder="test_invalid_method",
                                           seed=123)
        isdir("test_invalid_method") && rm("test_invalid_method", recursive=true)
    end

    @testset "Single-step requires genotypes" begin
        model = build_model("y1 = intercept", 1.0)
        @test_throws ErrorException runMCMC(model, phenotypes,
                                           chain_length=10,
                                           single_step_analysis=true,
                                           output_folder="test_ss_no_geno",
                                           seed=123)
        isdir("test_ss_no_geno") && rm("test_ss_no_geno", recursive=true)
    end

    @testset "output_samples_frequency validation" begin
        global geno = get_genotypes(genofile, 1.0, separator=',', method="BayesC")
        model = build_model("y1 = intercept + geno", 1.0)
        @test_throws ErrorException runMCMC(model, phenotypes,
                                           chain_length=10,
                                           output_samples_frequency=0,
                                           output_folder="test_freq_zero",
                                           seed=123)
        isdir("test_freq_zero") && rm("test_freq_zero", recursive=true)
    end

    @testset "Describe model" begin
        model = build_model("y1 = intercept + x1", 1.0)
        # describe() should run without error
        @test describe(model) === nothing
    end

    @testset "Multiple Bayesian methods load correctly" begin
        for method in ["BayesA", "BayesB", "BayesC", "BayesL", "RR-BLUP", "GBLUP"]
            global geno = get_genotypes(genofile, 1.0, separator=',', method=method)
            model = build_model("y1 = intercept + geno", 1.0)
            @test model.nModels == 1
        end
    end
end
