# Common Testing Patterns for JWAS.jl

This guide shows specific testing patterns for JWAS functionality.

## Table of Contents
1. [Model Testing](#model-testing)
2. [MCMC Testing](#mcmc-testing)
3. [Genotype Testing](#genotype-testing)
4. [Pedigree Testing](#pedigree-testing)
5. [Output Testing](#output-testing)
6. [Error Testing](#error-testing)
7. [Performance Testing](#performance-testing)

---

## Model Testing

### Test Model Construction
```julia
@testset "Model Construction" begin
    # Single trait
    model = build_model("y = intercept + x1", 1.0)
    @test model.nModels == 1
    @test :y in model.lhsVec
    
    # Multi-trait
    R = [1.0 0.5; 0.5 1.0]
    model = build_model("y1 = intercept\ny2 = intercept", R)
    @test model.nModels == 2
    @test size(model.R.val) == (2, 2)
end
```

### Test Covariate Assignment
```julia
@testset "Covariates" begin
    model = build_model("y = intercept + age + sex", 1.0)
    set_covariate(model, "age")
    
    @test :age in model.covVec
    @test :sex ∉ model.covVec  # sex is still a factor
end
```

### Test Random Effects
```julia
@testset "Random Effects" begin
    model = build_model("y = intercept + herd", 1.0)
    set_random(model, "herd", G1)
    
    @test length(model.rndTrmVec) == 1
    @test model.rndTrmVec[1].term_array == ["1:herd"]
end
```

### Test Pedigree Random Effects
```julia
@testset "Pedigree Effects" begin
    ped = get_pedigree("pedigree.txt", separator=",")
    model = build_model("y = intercept + animal", 1.0)
    set_random(model, "animal", ped, G)
    
    @test model.pedTrmVec != 0
    @test length(model.rndTrmVec) > 0
end
```

---

## MCMC Testing

### Test Basic MCMC Run
```julia
@testset "Basic MCMC" begin
    model = build_model("y = intercept", 1.0)
    output = runMCMC(model, data, 
                    chain_length=50,
                    output_folder="test_temp",
                    seed=123)
    
    @test haskey(output, "location parameters")
    @test haskey(output, "residual variance")
    @test size(output["location parameters"], 1) > 0
    
    rm("test_temp", recursive=true)
end
```

### Test MCMC with Genotypes
```julia
@testset "MCMC with Genotypes" begin
    geno = get_genotypes("genotypes.txt", 1.0, method="BayesC")
    model = build_model("y = intercept + geno", 1.0)
    
    output = runMCMC(model, data,
                    chain_length=50,
                    output_folder="test_geno",
                    seed=123)
    
    @test haskey(output, "marker effects geno")
    @test size(output["marker effects geno"], 1) > 0
    
    rm("test_geno", recursive=true)
end
```

### Test Single-Step Analysis
```julia
@testset "Single-Step SSGBLUP" begin
    ped = get_pedigree("pedigree.txt")
    geno = get_genotypes("genotypes.txt", 1.0, method="BayesC")
    model = build_model("y = intercept + animal + geno", 1.0)
    set_random(model, "animal", ped, G)
    
    output = runMCMC(model, data,
                    single_step_analysis=true,
                    pedigree=ped,
                    chain_length=50,
                    output_folder="test_ss",
                    seed=123)
    
    @test output isa Dict
    rm("test_ss", recursive=true)
end
```

### Test Reproducibility
```julia
@testset "MCMC Reproducibility" begin
    model1 = build_model("y = intercept", 1.0)
    model2 = build_model("y = intercept", 1.0)
    
    out1 = runMCMC(model1, data, chain_length=50,
                  output_folder="temp1", seed=999)
    out2 = runMCMC(model2, data, chain_length=50,
                  output_folder="temp2", seed=999)
    
    # Should get identical results with same seed
    @test out1["residual variance"][1, :Estimate] ≈ 
          out2["residual variance"][1, :Estimate] atol=1e-10
    
    rm("temp1", recursive=true)
    rm("temp2", recursive=true)
end
```

### Test Different Methods
```julia
@testset "Bayesian Methods" begin
    methods = ["BayesA", "BayesB", "BayesC", "RR-BLUP", "BayesL", "GBLUP"]
    
    for method in methods
        @testset "$method" begin
            geno = get_genotypes("genotypes.txt", 1.0, method=method)
            @test geno.method == method
            
            model = build_model("y = intercept + geno", 1.0)
            output = runMCMC(model, data,
                           chain_length=10,
                           output_folder="test_$method",
                           seed=123)
            
            @test output isa Dict
            rm("test_$method", recursive=true)
        end
    end
end
```

---

## Genotype Testing

### Test Genotype Loading from Different Sources
```julia
@testset "Genotype Loading" begin
    @testset "From file with header" begin
        geno = get_genotypes("genotypes.txt", 1.0,
                           separator=',', header=true)
        @test geno.nMarkers > 0
        @test geno.nObs > 0
        @test length(geno.markerID) == geno.nMarkers
    end
    
    @testset "From DataFrame" begin
        df = CSV.read("genotypes.txt", DataFrame)
        geno = get_genotypes(df, 1.0, header=true)
        @test geno.nMarkers > 0
    end
    
    @testset "From Array" begin
        arr = readdlm("genotypes.txt", ',')[2:end, 2:end]
        geno = get_genotypes(arr, 1.0)
        @test geno.nObs == size(arr, 1)
        @test geno.nMarkers == size(arr, 2)
    end
end
```

### Test Quality Control
```julia
@testset "Genotype QC" begin
    geno_with_qc = get_genotypes("genotypes.txt", 1.0,
                                 quality_control=true,
                                 MAF=0.05)
    geno_no_qc = get_genotypes("genotypes.txt", 1.0,
                               quality_control=false)
    
    # QC should remove some markers
    @test geno_with_qc.nMarkers <= geno_no_qc.nMarkers
    
    # All markers should have MAF >= threshold
    for freq in geno_with_qc.alleleFreq
        @test 0.05 <= freq <= 0.95
    end
end
```

### Test Centering
```julia
@testset "Genotype Centering" begin
    geno_centered = get_genotypes("genotypes.txt", 1.0,
                                  center=true)
    geno_raw = get_genotypes("genotypes.txt", 1.0,
                             center=false)
    
    @test geno_centered.centered == true
    @test geno_raw.centered == false
    
    # Centered genotypes should have mean ≈ 0
    if geno_centered.centered
        col_means = mean(geno_centered.genotypes, dims=1)
        @test all(abs.(col_means) .< 1e-10)
    end
end
```

### Test GBLUP with GRM
```julia
@testset "GBLUP with GRM" begin
    # Test with genomic relationship matrix
    geno = get_genotypes("genotypes.txt", 1.0, method="GBLUP")
    
    @test geno.isGRM == true
    @test issymmetric(geno.genotypes)
    @test isposdef(geno.genotypes)
end
```

---

## Pedigree Testing

### Test Pedigree Loading
```julia
@testset "Pedigree Loading" begin
    ped = get_pedigree("pedigree.txt", separator=",", header=true)
    
    @test length(ped.idMap) > 0
    @test length(ped.IDs) == length(ped.idMap)
    @test all(node.seqID > 0 for node in values(ped.idMap))
end
```

### Test Inbreeding Calculation
```julia
@testset "Inbreeding Coefficients" begin
    ped = get_pedigree("pedigree.txt", separator=",")
    
    # All inbreeding coefficients should be non-negative
    for (id, node) in ped.idMap
        @test node.f >= 0.0
        @test node.f <= 1.0  # Cannot exceed 1
    end
    
    # Founders should have f = 0
    for (id, node) in ped.idMap
        if node.sire == "missing" && node.dam == "missing"
            @test node.f == 0.0
        end
    end
end
```

### Test Numerator Relationship Matrix
```julia
@testset "A Matrix" begin
    ped = get_pedigree("pedigree.txt", separator=",")
    A = PedModule.AInverse(ped)
    
    @test size(A, 1) == size(A, 2)
    @test size(A, 1) == length(ped.idMap)
    @test issymmetric(A)
end
```

---

## Output Testing

### Test Output Files Created
```julia
@testset "Output Files" begin
    model = build_model("y = intercept", 1.0)
    output = runMCMC(model, data,
                    chain_length=50,
                    output_samples_frequency=10,
                    output_folder="test_out")
    
    @test isdir("test_out")
    @test isfile("test_out/location_parameters.txt")
    @test isfile("test_out/residual_variance.txt")
    @test isfile("test_out/MCMC_samples_residual_variance.txt")
    
    rm("test_out", recursive=true)
end
```

### Test EBV Output
```julia
@testset "EBV Output" begin
    geno = get_genotypes("genotypes.txt", 1.0)
    model = build_model("y = intercept + geno", 1.0)
    outputEBV(model, geno.obsID)
    
    output = runMCMC(model, data,
                    chain_length=50,
                    outputEBV=true,
                    output_folder="test_ebv")
    
    @test haskey(output, "EBV_y")
    @test "ID" in names(output["EBV_y"])
    @test "EBV" in names(output["EBV_y"])
    @test "PEV" in names(output["EBV_y"])
    
    rm("test_ebv", recursive=true)
end
```

### Test MCMC Samples
```julia
@testset "MCMC Samples" begin
    model = build_model("y = intercept + herd", 1.0)
    set_random(model, "herd")
    outputMCMCsamples(model, "herd")
    
    output = runMCMC(model, data,
                    chain_length=100,
                    output_samples_frequency=10,
                    output_folder="test_samples")
    
    samples_file = "test_samples/MCMC_samples_1.herd.txt"
    @test isfile(samples_file)
    
    samples = readdlm(samples_file, ',', skipstart=1)
    @test size(samples, 1) == 10  # 100 / 10 = 10 samples
    
    rm("test_samples", recursive=true)
end
```

---

## Error Testing

### Test Invalid Inputs
```julia
@testset "Error Handling" begin
    @testset "Negative variance" begin
        @test_throws ErrorException build_model("y = intercept", -1.0)
    end
    
    @testset "Non-PD covariance matrix" begin
        R = [1.0 2.0; 2.0 1.0]  # Not positive definite
        @test_throws ErrorException build_model("y1=intercept\ny2=intercept", R)
    end
    
    @testset "Invalid model equation" begin
        @test_throws ErrorException build_model("", 1.0)
    end
    
    @testset "Mismatched dimensions" begin
        R = [1.0 0.5; 0.5 1.0]  # 2x2
        # But only one trait in equation
        @test_throws ErrorException build_model("y = intercept", R)
    end
end
```

### Test Warning Conditions
```julia
@testset "Warnings" begin
    # Test that genotypes out of range produce warnings
    @test_logs (:warn, r"genotype scores") begin
        # Create genotypes with values > 2
        bad_geno = [0.0 1.0 3.0; 1.0 2.0 1.0]
        get_genotypes(bad_geno, 1.0, quality_control=false)
    end
end
```

---

## Performance Testing

### Test Execution Time
```julia
@testset "Performance" begin
    geno = get_genotypes("genotypes.txt", 1.0, method="RR-BLUP")
    model = build_model("y = intercept + geno", 1.0)
    
    # Test that analysis completes in reasonable time
    time = @elapsed begin
        output = runMCMC(model, data,
                       chain_length=100,
                       output_folder="perf_test")
    end
    
    @test time < 60  # Should complete in < 60 seconds
    rm("perf_test", recursive=true)
end
```

### Test Memory Usage
```julia
@testset "Memory Efficiency" begin
    # Test single vs double precision
    geno_single = get_genotypes("genotypes.txt", 1.0,
                                double_precision=false)
    geno_double = get_genotypes("genotypes.txt", 1.0,
                                double_precision=true)
    
    @test eltype(geno_single.genotypes) == Float32
    @test eltype(geno_double.genotypes) == Float64
    
    # Single precision should use less memory
    mem_single = Base.summarysize(geno_single.genotypes)
    mem_double = Base.summarysize(geno_double.genotypes)
    @test mem_single < mem_double
end
```

### Benchmark Different Methods
```julia
using BenchmarkTools

@testset "Method Comparison" begin
    methods = ["RR-BLUP", "BayesC"]
    times = Dict{String, Float64}()
    
    for method in methods
        geno = get_genotypes("genotypes.txt", 1.0, method=method)
        model = build_model("y = intercept + geno", 1.0)
        
        t = @belapsed runMCMC($model, $data,
                             chain_length=100,
                             output_folder="bench_$method",
                             seed=123)
        times[method] = t
        rm("bench_$method", recursive=true)
    end
    
    # RR-BLUP should be faster than BayesC
    @test times["RR-BLUP"] < times["BayesC"]
end
```

---

## Tips for Writing Good Tests

### 1. Use Descriptive Test Names
```julia
# Good
@testset "BayesC with π estimation" begin ... end

# Bad
@testset "Test 1" begin ... end
```

### 2. Test One Thing at a Time
```julia
# Good - separate tests
@testset "Model has correct number of traits" begin
    @test model.nModels == 2
end

@testset "Residual variance is set correctly" begin
    @test model.R.val == R
end

# Bad - tests multiple things
@testset "Model" begin
    @test model.nModels == 2
    @test model.R.val == R
    @test length(model.lhsVec) == 2
end
```

### 3. Clean Up After Tests
```julia
@testset "My Test" begin
    # ... test code ...
    
    # Always cleanup
    if isdir("test_output")
        rm("test_output", recursive=true)
    end
end
```

### 4. Use Appropriate Tolerances
```julia
# For float comparisons
@test result ≈ expected atol=1e-6

# For exact comparisons
@test count == expected_count
```

### 5. Test Edge Cases
```julia
@testset "Edge Cases" begin
    @test my_function(0) > 0
    @test my_function(Inf) < Inf
    @test_throws DomainError my_function(-1)
end
```

---

## Quick Reference

```julia
# Basic assertions
@test x == y          # Exact equality
@test x ≈ y          # Approximate equality
@test x ≈ y atol=ε   # With tolerance
@test x > y          # Comparison
@test x isa Type     # Type checking
@test_throws Error f() # Should error
@test_nowarn f()     # Should not warn/error
@test_logs (:warn,) f() # Should produce warning

# Testset structure
@testset "Name" begin
    # tests here
end

# Skip tests
@testset skip=true "Not ready" begin
    # tests skipped
end

# Test with setup/teardown
@testset "With cleanup" begin
    setup_test()
    try
        @test ...
    finally
        cleanup_test()
    end
end
```


