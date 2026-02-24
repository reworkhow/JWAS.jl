# Unit tests for BayesB and BayesA specific code paths
using Test, JWAS, DataFrames, CSV, JWAS.Datasets

phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
genofile = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstrings=["NA"])

@testset "BayesB single-trait" begin
    global geno = get_genotypes(genofile, 1.0, separator=',', method="BayesB")
    model = build_model("y1 = intercept + geno", 1.0)

    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    output_folder="test_bayesb",
                    seed=123)

    @test haskey(output, "marker effects geno")
    @test haskey(output, "location parameters")
    rm("test_bayesb", recursive=true)
end

@testset "BayesA single-trait (converted to BayesB with π=0)" begin
    global geno = get_genotypes(genofile, 1.0, separator=',', method="BayesA")
    model = build_model("y1 = intercept + geno", 1.0)

    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    output_folder="test_bayesa",
                    seed=123)

    @test haskey(output, "marker effects geno")
    @test haskey(output, "location parameters")
    rm("test_bayesa", recursive=true)
end

@testset "BayesL (Lasso) single-trait" begin
    global geno = get_genotypes(genofile, 1.0, separator=',', method="BayesL")
    model = build_model("y1 = intercept + geno", 1.0)

    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    output_folder="test_bayesl",
                    seed=123)

    @test haskey(output, "marker effects geno")
    @test haskey(output, "location parameters")
    rm("test_bayesl", recursive=true)
end

@testset "GBLUP" begin
    global geno = get_genotypes(genofile, 1.0, separator=',', method="GBLUP")
    model = build_model("y1 = intercept + geno", 1.0)

    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    output_folder="test_gblup",
                    seed=123)

    @test haskey(output, "location parameters")
    rm("test_gblup", recursive=true)
end
