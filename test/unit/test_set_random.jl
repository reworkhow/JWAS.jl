# Unit tests for set_random() — random effects configuration
using Test, JWAS, DataFrames, CSV, JWAS.Datasets

@testset "set_random with i.i.d. effects" begin
    phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
    phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstrings=["NA"])

    @testset "Single-trait i.i.d. random effect" begin
        model = build_model("y1 = intercept + x1", 1.0)
        set_random(model, "x1", 0.5)
        @test length(model.rndTrmVec) == 1
        @test model.rndTrmVec[1].randomType == "I"
    end

    @testset "Multi-trait i.i.d. random effect" begin
        R = [1.0 0.5; 0.5 1.0]
        model = build_model("y1 = intercept + x1\ny2 = intercept + x1", R)
        G = [0.5 0.1; 0.1 0.5]
        set_random(model, "x1", G)
        @test length(model.rndTrmVec) == 1
        @test model.rndTrmVec[1].randomType == "I"
        @test length(model.rndTrmVec[1].term_array) == 2
    end

    @testset "Invalid covariance matrix for random effect" begin
        model = build_model("y1 = intercept + x1", 1.0)
        @test_throws ErrorException set_random(model, "x1", -0.5)
    end
end

@testset "set_random with pedigree" begin
    pedfile = Datasets.dataset("pedigree.txt", dataset_name="demo_7animals")
    ped = get_pedigree(pedfile, separator=",", header=true)

    @testset "Single-trait polygenic effect" begin
        model = build_model("y1 = intercept + animal", 1.0)
        set_random(model, "animal", ped, 1.6)
        @test model.ped !== nothing
        @test model.ped != 0
        @test length(model.rndTrmVec) == 1
        @test model.rndTrmVec[1].randomType == "A"
        @test length(model.pedTrmVec) > 0
    end

    @testset "Cannot add pedigree twice" begin
        model = build_model("y1 = intercept + animal", 1.0)
        set_random(model, "animal", ped, 1.6)
        @test_throws ErrorException set_random(model, "animal", ped, 1.6)
    end

    @testset "Pedigree-based MCMC runs" begin
        phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
        phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstrings=["NA"])

        model = build_model("y1 = intercept + ID", 1.0)
        set_random(model, "ID", ped, 1.6)

        output = runMCMC(model, phenotypes,
                        chain_length=50,
                        output_folder="test_set_random_ped",
                        seed=123)
        @test haskey(output, "polygenic effects covariance matrix")
        @test haskey(output, "location parameters")
        rm("test_set_random_ped", recursive=true)
    end
end
