# Unit tests for outputEBV(), getEBV(), and output functions
using Test, JWAS, DataFrames, CSV, JWAS.Datasets

phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
genofile = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstrings=["NA"])

@testset "outputEBV and EBV results" begin
    @testset "EBV output with genotypes" begin
        global geno = get_genotypes(genofile, 1.0, separator=',', method="BayesC")
        model = build_model("y1 = intercept + geno", 1.0)
        outputEBV(model, geno.obsID)

        output = runMCMC(model, phenotypes,
                        chain_length=100,
                        burnin=20,
                        output_samples_frequency=10,
                        outputEBV=true,
                        output_folder="test_ebv_output",
                        seed=123)

        @test haskey(output, "EBV_y1")
        ebv_df = output["EBV_y1"]
        @test :ID in propertynames(ebv_df)
        @test :EBV in propertynames(ebv_df)
        @test :PEV in propertynames(ebv_df)
        @test size(ebv_df, 1) > 0
        @test all(ebv_df.PEV .>= 0)  # PEV should be non-negative

        rm("test_ebv_output", recursive=true)
    end

    @testset "EBV output with heritability" begin
        global geno = get_genotypes(genofile, 1.0, separator=',', method="RR-BLUP")
        model = build_model("y1 = intercept + geno", 1.0)

        output = runMCMC(model, phenotypes,
                        chain_length=100,
                        burnin=20,
                        output_samples_frequency=10,
                        outputEBV=true,
                        output_heritability=true,
                        output_folder="test_ebv_h2",
                        seed=123)

        @test haskey(output, "EBV_y1")
        @test haskey(output, "heritability")
        @test haskey(output, "genetic_variance")

        h2 = output["heritability"]
        @test :Estimate in propertynames(h2)
        @test all(0 .<= h2.Estimate .<= 1)

        rm("test_ebv_h2", recursive=true)
    end
end

@testset "outputMCMCsamples for location parameters" begin
    global geno = get_genotypes(genofile, 1.0, separator=',', method="BayesC")
    model = build_model("y1 = intercept + geno", 1.0)
    outputMCMCsamples(model, "intercept")

    output = runMCMC(model, phenotypes,
                    chain_length=50,
                    output_samples_frequency=10,
                    output_folder="test_mcmc_samples",
                    seed=123)

    @test isfile("test_mcmc_samples/MCMC_samples_y1.intercept.txt")
    rm("test_mcmc_samples", recursive=true)
end
