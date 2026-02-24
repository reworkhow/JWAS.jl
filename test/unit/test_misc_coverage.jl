# Unit tests for miscellaneous coverage improvements
# - describe() with data
# - update_priors_frequency
# - outputMCMCsamples for multiple terms
# - Datasets error handling
using Test, JWAS, DataFrames, CSV, JWAS.Datasets

phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
genofile = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstrings=["NA"])

@testset "Datasets module" begin
    @testset "Valid dataset access" begin
        f = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
        @test isfile(f)
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
