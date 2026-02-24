# Unit tests for multi-trait MCMC with different Bayesian methods
using Test, JWAS, DataFrames, CSV, JWAS.Datasets

phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
genofile = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
pedfile = Datasets.dataset("pedigree.txt", dataset_name="demo_7animals")
phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstrings=["NA"])

R = [1.0 0.5; 0.5 1.0]

@testset "Multi-trait BayesC" begin
    G = [1.0 0.5; 0.5 1.0]
    global geno = get_genotypes(genofile, G, separator=',', method="BayesC")
    model = build_model("y1 = intercept + geno\ny2 = intercept + geno", R)

    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    output_folder="test_mt_bayesc",
                    seed=123)

    @test haskey(output, "location parameters")
    @test haskey(output, "residual variance")
    @test haskey(output, "marker effects geno")
    @test size(output["residual variance"], 1) == 4  # 2x2 covariance flattened
    rm("test_mt_bayesc", recursive=true)
end

@testset "Multi-trait RR-BLUP" begin
    G = [1.0 0.5; 0.5 1.0]
    global geno = get_genotypes(genofile, G, separator=',', method="RR-BLUP")
    model = build_model("y1 = intercept + geno\ny2 = intercept + geno", R)

    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    output_folder="test_mt_rrblup",
                    seed=123)

    @test haskey(output, "location parameters")
    @test haskey(output, "residual variance")
    rm("test_mt_rrblup", recursive=true)
end

@testset "Multi-trait with pedigree" begin
    G = [1.0 0.5; 0.5 1.0]
    ped = get_pedigree(pedfile, separator=",", header=true)
    model = build_model("y1 = intercept + ID\ny2 = intercept + ID", R)
    set_random(model, "ID", ped, G)

    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    output_folder="test_mt_ped",
                    seed=123)

    @test haskey(output, "polygenic effects covariance matrix")
    @test haskey(output, "location parameters")
    @test haskey(output, "residual variance")
    rm("test_mt_ped", recursive=true)
end

@testset "Multi-trait with i.i.d. random effect" begin
    G_geno = [1.0 0.5; 0.5 1.0]
    global geno = get_genotypes(genofile, G_geno, separator=',', method="RR-BLUP")
    model = build_model("y1 = intercept + x2 + geno\ny2 = intercept + x2 + geno", R)
    G = [0.5 0.1; 0.1 0.5]
    set_random(model, "x2", G)

    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    output_folder="test_mt_iid",
                    seed=123)

    @test haskey(output, "location parameters")
    @test haskey(output, "residual variance")
    rm("test_mt_iid", recursive=true)
end
