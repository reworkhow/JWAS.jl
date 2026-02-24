# Unit tests for GWAS window-based analysis
using Test, JWAS, DataFrames, CSV, JWAS.Datasets

phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
genofile = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
mapfile = Datasets.dataset("map.txt", dataset_name="demo_7animals")
phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstrings=["NA"])

# Run a short MCMC to generate marker effect samples
global geno = get_genotypes(genofile, 1.0, separator=',', method="BayesC")
model = build_model("y1 = intercept + geno", 1.0)
outputEBV(model, geno.obsID)

output = runMCMC(model, phenotypes,
                chain_length=100,
                burnin=20,
                output_samples_frequency=10,
                outputEBV=true,
                output_folder="test_gwas_win",
                seed=123)

marker_file = "test_gwas_win/MCMC_samples_marker_effects_geno_y1.txt"

@testset "GWAS window-based analysis" begin
    @testset "Non-overlapping windows with map file" begin
        result = GWAS(model, mapfile, marker_file,
                     window_size="1 Mb", header=true)
        @test length(result) >= 1
        gwas_df = result[1]
        @test "WPPA" in names(gwas_df)
        @test "chr" in names(gwas_df)
        @test "numSNP" in names(gwas_df)
        @test "estimateGenVar" in names(gwas_df)
        @test all(0 .<= gwas_df.WPPA .<= 1)
    end

    @testset "Sliding windows" begin
        result = GWAS(model, mapfile, marker_file,
                     window_size="1 Mb", sliding_window=true, header=true)
        @test length(result) >= 1
        gwas_df = result[1]
        @test "WPPA" in names(gwas_df)
    end

    @testset "Window with output_winVarProps" begin
        result, winVarProps = GWAS(model, mapfile, marker_file,
                                  window_size="1 Mb", output_winVarProps=true, header=true)
        @test length(result) >= 1
        @test length(winVarProps) >= 1
    end

    @testset "Different threshold" begin
        result = GWAS(model, mapfile, marker_file,
                     window_size="1 Mb", threshold=0.01, header=true)
        @test length(result) >= 1
        gwas_df = result[1]
        @test "WPPA" in names(gwas_df)
    end
end

# Cleanup
rm("test_gwas_win", recursive=true)
for f in filter(x -> startswith(x, "GWAS_") || startswith(x, "MCMC_samples_local"), readdir())
    rm(f)
end
