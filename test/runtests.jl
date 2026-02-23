# Comprehensive Test Suite for JWAS.jl
# Combines unit tests with conditional integration tests
# Uses unique test folder for isolated, reproducible testing

using Test
using JWAS, DataFrames, CSV, JWAS.Datasets
using Statistics
using Random

# ============================================================================
# SETUP: Cleanup and Create Unique Test Environment
# ============================================================================

println("="^70)
println("JWAS.jl Comprehensive Test Suite")
println("="^70)

# Clean up old test artifacts from previous runs
println("\nCleaning up old test artifacts...")
for old_folder in ["results", "results1", "results2", "results3", "results4", "results5",
                   "test_results_temp", "test_results_temp1", "test_results_temp2", "test_results_temp3",
                   "temp1", "temp11", "temp12", "temp13", "temp2",
                   "gwas_test", "gwas_test1", "gwas_test2", "gwas_test3", "gwas_test4",
                   "missing_test", "mytest"]
    if isdir(old_folder)
        rm(old_folder, recursive=true)
        println("  Removed folder: $old_folder")
    end
end

# Remove old ID files
for old_file in ["IDs_for_individuals_with_genotypes.txt", 
                 "IDs_for_individuals_with_pedigree.txt",
                 "IDs_for_individuals_with_phenotypes.txt"]
    if isfile(old_file)
        rm(old_file)
        println("  Removed file: $old_file")
    end
end
println("Cleanup complete.\n")

# Create unique test directory for this run
test_home = "test_run_" * string(rand(10000:99999))
original_dir = pwd()
mkdir(test_home)
cd(test_home)

println("="^70)
println("Running tests in: $test_home")
println("="^70)
println()

# Configuration flags from environment variables
RUN_INTEGRATION_TESTS = get(ENV, "RUN_INTEGRATION_TESTS", "false") == "true"
println("Test Configuration:")
println("  Unit Tests:        Always run ✓")
println("  Integration Tests: $(RUN_INTEGRATION_TESTS ? "Enabled ✓" : "Disabled (set RUN_INTEGRATION_TESTS=true)")")
println("="^70)
println()

# ============================================================================
# MAIN TEST SUITE
# ============================================================================

@testset "JWAS.jl Full Test Suite" begin
    
    # ========================================================================
    # UNIT TESTS (Always run - fast feedback)
    # ========================================================================
    @testset "Unit Tests" begin
        
        # ====================================================================
        # Test 1: Model Building
        # ====================================================================
        @testset "Model Building" begin
            @testset "Single-trait model" begin
                model = build_model("y = intercept + x1", 1.0)
                @test model.nModels == 1
                @test length(model.lhsVec) == 1
                @test model.lhsVec[1] == :y
                @test model.R.val == 1.0
            end
            
            @testset "Multi-trait model" begin
                R = [1.0 0.5; 0.5 1.0]
                model = build_model("y1 = intercept\ny2 = intercept", R)
                @test model.nModels == 2
                @test length(model.lhsVec) == 2
                @test model.lhsVec[1] == :y1
                @test model.lhsVec[2] == :y2
                @test size(model.R.val) == (2, 2)
            end
            
            @testset "Set covariate" begin
                model = build_model("y = intercept + x1 + x2", 1.0)
                set_covariate(model, "x1")
                @test :x1 in model.covVec
                @test :x2 ∉ model.covVec
            end

            @testset "NNMM keywords removed" begin
                @test_throws MethodError build_model("y = intercept + x1", 1.0; nonlinear_function="tanh")
                @test_throws MethodError build_model("y = intercept + x1", 1.0; num_hidden_nodes=3)
                @test_throws MethodError build_model("y = intercept + x1", 1.0; latent_traits=["z1"])
            end
        end
        
        # ====================================================================
        # Test 2: Genotype Loading
        # ====================================================================
        @testset "Genotype Loading" begin
            genofile = Datasets.dataset("genotypes.txt", dataset_name="example")
            
            @testset "Load from file with header" begin
                geno = get_genotypes(genofile, 1.0, separator=',', header=true, 
                                    method="BayesC")
                @test geno.nMarkers > 0
                @test geno.nObs > 0
                @test geno.centered == true
                @test geno.method == "BayesC"
                @test length(geno.obsID) == geno.nObs
                @test length(geno.markerID) == geno.nMarkers
            end
            
            @testset "Load with different methods" begin
                methods = ["BayesA", "BayesB", "BayesC", "RR-BLUP", "BayesL", "GBLUP"]
                for method in methods
                    geno = get_genotypes(genofile, 1.0, separator=',', 
                                        method=method)
                    @test geno.method == method
                    @test geno.nMarkers > 0
                end
            end
            
            @testset "Quality control" begin
                geno_with_qc = get_genotypes(genofile, 1.0, separator=',',
                                            quality_control=true, MAF=0.01)
                geno_no_qc = get_genotypes(genofile, 1.0, separator=',',
                                          quality_control=false)
                # QC should remove some markers or keep same
                @test geno_with_qc.nMarkers <= geno_no_qc.nMarkers
            end
        end
        
        # ====================================================================
        # Test 3: Pedigree Module
        # ====================================================================
        @testset "Pedigree Module" begin
            pedfile = Datasets.dataset("pedigree.txt", dataset_name="example")
            
            @testset "Load pedigree" begin
                ped = get_pedigree(pedfile, separator=",", header=true)
                @test length(ped.idMap) > 0
                @test length(ped.IDs) == length(ped.idMap)
                @test ped.currentID > 1
            end
            
            @testset "Pedigree information" begin
                ped = get_pedigree(pedfile, separator=",", header=true)
                # Check that all IDs are unique
                @test length(unique(ped.IDs)) == length(ped.IDs)
                # Check that inbreeding coefficients are calculated
                for (id, node) in ped.idMap
                    @test node.f >= 0.0  # Inbreeding coefficient should be non-negative
                end
            end
        end
        
        # ====================================================================
        # Test 4: MCMC Functionality (Short runs)
        # ====================================================================
        @testset "MCMC Functionality" begin
            phenofile = Datasets.dataset("phenotypes.txt", dataset_name="example")
            genofile = Datasets.dataset("genotypes.txt", dataset_name="example")
            phenotypes = CSV.read(phenofile, DataFrame, delim=',', 
                                 missingstrings=["NA"])
            
            @testset "Single-trait BayesC" begin
                global geno = get_genotypes(genofile, 1.0, separator=',', method="BayesC")
                model = build_model("y1 = intercept + geno", 1.0)
                
                output = runMCMC(model, phenotypes, 
                               chain_length=50,
                               burnin=10,
                               output_samples_frequency=10,
                               output_folder="results",
                               seed=123)
                
                @test haskey(output, "location parameters")
                @test haskey(output, "residual variance")
                @test haskey(output, "marker effects geno")
                @test size(output["location parameters"], 1) > 0
                
                # Cleanup
                rm("results", recursive=true)
            end
            
            @testset "Output folder creation" begin
                global geno = get_genotypes(genofile, 1.0, separator=',', method="RR-BLUP")
                model = build_model("y1 = intercept + geno", 1.0)
                
                output = runMCMC(model, phenotypes,
                               chain_length=50,
                               output_folder="test_results_temp",
                               seed=123)
                
                @test isdir("test_results_temp")
                @test isfile("test_results_temp/location_parameters.txt")
                @test isfile("test_results_temp/residual_variance.txt")
                
                # Cleanup
                rm("test_results_temp", recursive=true)
            end
            
            @testset "Reproducibility with seed" begin
                global geno = get_genotypes(genofile, 1.0, separator=',', method="RR-BLUP")
                model1 = build_model("y1 = intercept + geno", 1.0)
                out1 = runMCMC(model1, phenotypes, chain_length=50, 
                              output_folder="temp1", seed=999)
                
                global geno = get_genotypes(genofile, 1.0, separator=',', method="RR-BLUP")
                model2 = build_model("y1 = intercept + geno", 1.0)
                out2 = runMCMC(model2, phenotypes, chain_length=50,
                              output_folder="temp2", seed=999)
                
                # Results should be identical with same seed
                @test out1["residual variance"][1, :Estimate] ≈ 
                      out2["residual variance"][1, :Estimate] atol=1e-10
                
                # Cleanup
                rm("temp1", recursive=true)
                rm("temp2", recursive=true)
            end
        end
        
        # ====================================================================
        # Test 5: GWAS Module
        # ====================================================================
        @testset "GWAS Module" begin
            phenofile = Datasets.dataset("phenotypes.txt", dataset_name="example")
            genofile = Datasets.dataset("genotypes.txt", dataset_name="example")
            mapfile = Datasets.dataset("map.txt", dataset_name="example")
            phenotypes = CSV.read(phenofile, DataFrame, delim=',', 
                                 missingstrings=["NA"])
            
            @testset "Model frequency calculation" begin
                global geno = get_genotypes(genofile, 1.0, separator=',', method="BayesC")
                model = build_model("y1 = intercept + geno", 1.0)
                
                output = runMCMC(model, phenotypes, chain_length=50,
                               output_folder="gwas_test", seed=123)
                
                gwas_result = GWAS("gwas_test/MCMC_samples_marker_effects_geno_y1.txt")
                
                @test size(gwas_result, 1) > 0
                @test "marker_ID" in names(gwas_result)
                @test "modelfrequency" in names(gwas_result)
                @test all(0 .<= gwas_result.modelfrequency .<= 1)
                
                # Cleanup
                rm("gwas_test", recursive=true)
            end
        end
        
        # ====================================================================
        # Test 6: Edge Cases and Error Handling
        # ====================================================================
        @testset "Edge Cases" begin
            @testset "Invalid residual variance" begin
                # Negative variance should error
                @test_throws ErrorException build_model("y = intercept", -1.0)
            end
            
            @testset "Invalid covariance matrix" begin
                # Non-positive definite matrix should error
                R = [1.0 2.0; 2.0 1.0]  # Not positive definite
                @test_throws ErrorException build_model("y1 = intercept\ny2 = intercept", R)
            end
            
            @testset "Missing phenotype handling" begin
                phenofile = Datasets.dataset("phenotypes.txt", dataset_name="example")
                phenotypes = CSV.read(phenofile, DataFrame, delim=',',
                                     missingstrings=["NA"])
                
                # Should handle missing values without error
                model = build_model("y1 = intercept", 1.0)
                output = runMCMC(model, phenotypes, chain_length=10,
                                output_folder="missing_test", seed=123)
                @test haskey(output, "location parameters")
                
                # Cleanup
                rm("missing_test", recursive=true)
            end
        end
        
        # ====================================================================
        # Test 7: Data Types and Precision
        # ====================================================================
        @testset "Data Types" begin
            genofile = Datasets.dataset("genotypes.txt", dataset_name="example")
            
            @testset "Single precision (default)" begin
                geno = get_genotypes(genofile, 1.0, separator=',',
                                    double_precision=false)
                @test eltype(geno.genotypes) == Float32
            end
            
            @testset "Double precision" begin
                geno = get_genotypes(genofile, 1.0, separator=',',
                                    double_precision=true)
                @test eltype(geno.genotypes) == Float64
            end
        end
        
    end  # End of Unit Tests
    
    # ========================================================================
    # INTEGRATION TESTS (Conditional - your existing test files)
    # ========================================================================
    if RUN_INTEGRATION_TESTS
        println("\n" * "="^70)
        println("Running Integration Tests...")
        println("="^70)
        
        @testset "Integration Tests" begin

            @testset "Simulated Omics Data Tests" begin
                println("\n→ Running integration/test_simulated_omics.jl")
                integration_file = joinpath(original_dir, "test", "integration", "test_simulated_omics.jl")
                include(integration_file)
            end

        end
    else
        @info "Skipping integration tests (set RUN_INTEGRATION_TESTS=true to run)"
    end
    
end  # End of main testset

# ============================================================================
# CLEANUP: Return to original directory
# ============================================================================
cd(original_dir)

# Remove test directory (set to false to keep for debugging)
CLEANUP_ON_SUCCESS = true

println("\n" * "="^70)
println("Test Run Summary")
println("="^70)
println("  Test directory: $test_home")

if CLEANUP_ON_SUCCESS
    println("  Cleaning up test files...")
    if isdir(test_home)
        try
            rm(test_home, recursive=true, force=true)
            println("  ✓ Test directory cleaned up")
        catch e
            println("  ⚠ Could not remove test directory (files may be locked)")
            println("  Directory: $test_home")
            # Don't fail tests because of cleanup issues
        end
    else
        println("  ✓ Test directory already cleaned")
    end
else
    println("  Test files kept for inspection")
    println("  To clean up manually: rm(\"$test_home\", recursive=true)")
end

println("="^70)
println("✓ All tests completed successfully!")
println("="^70)
println()

# ============================================================================
# Usage Instructions
# ============================================================================
"""
# How to Run This Test Suite:

## Quick unit tests only (default - ~40 seconds):
julia --project=. test/runtests.jl

## With integration tests (~minutes):
RUN_INTEGRATION_TESTS=true julia --project=. test/runtests.jl

## Keep test files for debugging:
Edit CLEANUP_ON_SUCCESS = false in this file

# Features:
- ✓ Unique test folder per run (test_run_XXXXX)
- ✓ Automatic cleanup before and after
- ✓ No file conflicts between runs
- ✓ Conditional test execution
- ✓ Fast feedback for development
- ✓ Comprehensive testing for releases
"""
