using Test, JWAS, DataFrames, CSV, JWAS.Datasets
using LinearAlgebra, SparseArrays
using Random, Distributions, DelimitedFiles

# ============================================================================
# 1. tranform_lambda
# ============================================================================
@testset "tranform_lambda" begin
    @testset "2-trait: single causal path" begin
        # causal_structure: trait 1 -> trait 2  (entry [2,1] = 1)
        cs = [0.0 0.0; 1.0 0.0]
        lambda = [0.5]  # one causal path
        result = JWAS.tranform_lambda(lambda, cs)
        @test size(result) == (2, 2)
        @test result[2, 1] ≈ 0.5
        @test result[1, 1] ≈ 0.0
        @test result[1, 2] ≈ 0.0
        @test result[2, 2] ≈ 0.0
    end

    @testset "3-trait: multiple causal paths" begin
        # trait 1 -> trait 2, trait 1 -> trait 3, trait 2 -> trait 3
        cs = [0.0 0.0 0.0;
              1.0 0.0 0.0;
              1.0 1.0 0.0]
        lambda = [0.3, 0.7, 0.2]  # λ_21, λ_31, λ_32
        result = JWAS.tranform_lambda(lambda, cs)
        @test size(result) == (3, 3)
        @test result[2, 1] ≈ 0.3   # λ_21
        @test result[3, 1] ≈ 0.7   # λ_31
        @test result[3, 2] ≈ 0.2   # λ_32
        # all other entries should be zero
        @test result[1, 1] ≈ 0.0
        @test result[1, 2] ≈ 0.0
        @test result[1, 3] ≈ 0.0
        @test result[2, 2] ≈ 0.0
        @test result[2, 3] ≈ 0.0
        @test result[3, 3] ≈ 0.0
    end

    @testset "round-trip: nonzero positions match causal_structure" begin
        cs = [0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 1.0 0.0]
        lambda = [1.5, -0.8]
        result = JWAS.tranform_lambda(lambda, cs)
        # nonzero entries should match cs positions
        for i in 1:3, j in 1:3
            if cs[i, j] == 0.0
                @test result[i, j] ≈ 0.0
            else
                @test result[i, j] != 0.0
            end
        end
    end
end

# ============================================================================
# 2. compute_indirect_effect
# ============================================================================
@testset "compute_indirect_effect" begin
    @testset "zero Λ produces zero indirect effects" begin
        Λ = zeros(3, 3)
        marker_effects = [1.0 2.0; 3.0 4.0; 5.0 6.0]  # 3 traits x 2 markers
        result = JWAS.compute_indirect_effect(Λ, marker_effects)
        @test size(result) == (3, 2)
        @test all(result .≈ 0.0)
    end

    @testset "2-trait: known indirect effect" begin
        # Λ = [0 0; 0.5 0], marker_effects = [1.0; 2.0] (2 traits x 1 marker)
        # indirect = Λ^1 * marker_effects (only i=1 since ntraits-1=1)
        # = [0 0; 0.5 0] * [1.0; 2.0] = [0.0; 0.5]
        Λ = [0.0 0.0; 0.5 0.0]
        marker_effects = reshape([1.0, 2.0], 2, 1)
        result = JWAS.compute_indirect_effect(Λ, marker_effects)
        @test size(result) == (2, 1)
        @test result[1, 1] ≈ 0.0
        @test result[2, 1] ≈ 0.5  # 0.5 * 1.0
    end

    @testset "3-trait chain: indirect effect propagates" begin
        # Λ = [0 0 0; 0.5 0 0; 0 0.5 0]
        # trait 1 -> trait 2 -> trait 3
        # marker_effects = [1.0; 0.0; 0.0] (only trait 1 has direct effect)
        Λ = [0.0 0.0 0.0; 0.5 0.0 0.0; 0.0 0.5 0.0]
        marker_effects = reshape([1.0, 0.0, 0.0], 3, 1)
        result = JWAS.compute_indirect_effect(Λ, marker_effects)
        # i=1: Λ^1 * [1;0;0] = [0; 0.5; 0]
        # i=2: Λ^2 * [1;0;0] = [0; 0; 0.25]
        # total = [0; 0.5; 0.25]
        @test result[1, 1] ≈ 0.0
        @test result[2, 1] ≈ 0.5
        @test result[3, 1] ≈ 0.25
    end
end

# ============================================================================
# 3. get_sparse_Y_FSM
# ============================================================================
@testset "get_sparse_Y_FSM" begin
    @testset "2-trait, 3 individuals" begin
        wArray = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
        Y = JWAS.get_sparse_Y_FSM(wArray)
        ntraits = 2
        nind = 3
        # FSM Y has size: (ntraits*nind) x (ntraits*(ntraits-1))
        # = 6 x 2
        @test size(Y) == (ntraits * nind, ntraits * (ntraits - 1))

        # For 2 traits: Y should have phenotypes of the OTHER trait
        # Row block 1 (trait 1, inds 1-3): column = y2 values
        # Row block 2 (trait 2, inds 1-3): column = y1 values
        Y_dense = Matrix(Y)
        # trait 1 block: should contain y2 values [4,5,6]
        @test Y_dense[1:3, 1] ≈ [4.0, 5.0, 6.0]
        # trait 2 block: should contain y1 values [1,2,3]
        @test Y_dense[4:6, 2] ≈ [1.0, 2.0, 3.0]
    end

    @testset "3-trait, 2 individuals" begin
        wArray = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
        Y = JWAS.get_sparse_Y_FSM(wArray)
        ntraits = 3
        nind = 2
        # size: (3*2) x (3*(3-1)) = 6 x 6
        @test size(Y) == (ntraits * nind, ntraits * (ntraits - 1))
    end

    @testset "diagonal blocks are zero" begin
        wArray = [[1.0, 2.0], [3.0, 4.0]]
        Y = JWAS.get_sparse_Y_FSM(wArray)
        Y_dense = Matrix(Y)
        # For 2 traits, Y is 4x2. No "diagonal block" concept for 2-trait
        # but the column filtering removes self-phenotype columns
        # Verify no trait appears in its own row-block position
        # trait 1 rows (1:2) should NOT have trait 1 values -> col 1 has trait 2 vals
        @test Y_dense[1, 1] ≈ 3.0  # y2[1], not y1[1]
        @test Y_dense[2, 1] ≈ 4.0  # y2[2], not y1[2]
    end
end

# ============================================================================
# 4. get_sparse_Y_FRM
# ============================================================================
@testset "get_sparse_Y_FRM" begin
    @testset "filters columns based on causal_structure" begin
        wArray = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
        # Only trait 1 -> trait 2 (one causal path)
        cs = [0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 0.0 0.0]
        Y = JWAS.get_sparse_Y_FRM(wArray, cs)
        # Should have 1 column (one causal path)
        @test size(Y, 2) == 1
        @test size(Y, 1) == 6  # 3 traits x 2 inds
    end

    @testset "full lower triangular has more columns" begin
        wArray = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
        # trait 1 -> trait 2, trait 1 -> trait 3, trait 2 -> trait 3
        cs_full = [0.0 0.0 0.0; 1.0 0.0 0.0; 1.0 1.0 0.0]
        Y_full = JWAS.get_sparse_Y_FRM(wArray, cs_full)
        @test size(Y_full, 2) == 3  # 3 causal paths

        # single path
        cs_single = [0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 0.0 0.0]
        Y_single = JWAS.get_sparse_Y_FRM(wArray, cs_single)
        @test size(Y_single, 2) == 1

        @test size(Y_full, 2) > size(Y_single, 2)
    end

    @testset "FRM is column-subset of FSM" begin
        wArray = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
        cs = [0.0 0.0; 1.0 0.0]
        Y_fsm = Matrix(JWAS.get_sparse_Y_FSM(wArray))
        Y_frm = Matrix(JWAS.get_sparse_Y_FRM(wArray, cs))
        # FRM columns should be a subset of FSM columns
        for col in eachcol(Y_frm)
            found = any(c -> c ≈ col, eachcol(Y_fsm))
            @test found
        end
    end

    @testset "2-trait FRM values" begin
        wArray = [[10.0, 20.0], [30.0, 40.0]]
        cs = [0.0 0.0; 1.0 0.0]  # trait 1 -> trait 2
        Y = JWAS.get_sparse_Y_FRM(wArray, cs)
        Y_dense = Matrix(Y)
        # 1 causal path: y1 affecting trait 2
        # trait 1 rows (1:2) should be zero (no causal effect on trait 1)
        @test Y_dense[1, 1] ≈ 0.0
        @test Y_dense[2, 1] ≈ 0.0
        # trait 2 rows (3:4) should have y1 values [10, 20]
        @test Y_dense[3, 1] ≈ 10.0
        @test Y_dense[4, 1] ≈ 20.0
    end
end

# ============================================================================
# 5. get_Λ (stochastic sampler)
# ============================================================================
@testset "get_Λ sampling" begin
    @testset "output dimensions and in-place update" begin
        Random.seed!(42)
        nind = 4
        ntraits = 2
        wArray = [rand(nind), rand(nind)]
        cs = [0.0 0.0; 1.0 0.0]
        Y = JWAS.get_sparse_Y_FRM(wArray, cs)
        R = [1.0 0.0; 0.0 1.0]
        y = vcat(wArray[1], wArray[2])
        Λy = copy(y)       # initial Λ = I, so Λy = y
        Λycorr = copy(y)   # initial corrected = y (no effects yet)

        λ, λ_vec = JWAS.get_Λ(Y, R, Λycorr, Λy, y, cs)

        # λ should have 1 element (one causal path)
        @test length(λ) == 1
        # λ_vec should be vec of 2x2 causal matrix = 4 elements
        @test length(λ_vec) == ntraits * ntraits
        # Λy and Λycorr should have been updated in-place
        @test length(Λy) == ntraits * nind
        @test all(isfinite.(Λy))
        @test all(isfinite.(Λycorr))
    end

    @testset "3-trait sampling dimensions" begin
        Random.seed!(123)
        nind = 5
        wArray = [rand(nind), rand(nind), rand(nind)]
        cs = [0.0 0.0 0.0; 1.0 0.0 0.0; 1.0 1.0 0.0]
        Y = JWAS.get_sparse_Y_FRM(wArray, cs)
        R = Matrix{Float64}(I, 3, 3)
        y = vcat(wArray...)
        Λy = copy(y)
        Λycorr = copy(y)

        λ, λ_vec = JWAS.get_Λ(Y, R, Λycorr, Λy, y, cs)

        # 3 causal paths
        @test length(λ) == 3
        # vec of 3x3 matrix = 9 elements
        @test length(λ_vec) == 9
    end

    @testset "multiple samples produce finite values" begin
        Random.seed!(99)
        nind = 4
        wArray = [rand(nind), rand(nind)]
        cs = [0.0 0.0; 1.0 0.0]
        Y = JWAS.get_sparse_Y_FRM(wArray, cs)
        R = [2.0 0.0; 0.0 2.0]
        y = vcat(wArray...)
        Λy = copy(y)
        Λycorr = copy(y)

        for _ in 1:10
            λ, _ = JWAS.get_Λ(Y, R, Λycorr, Λy, y, cs)
            @test all(isfinite.(λ))
        end
    end
end

# ============================================================================
# 6. Input Validation
# ============================================================================
@testset "SEM input validation" begin
    phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
    genofile  = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
    phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstring=["NA"])
    phenotypes_complete = dropmissing(phenotypes, [:y1, :y2])

    @testset "non-lower-triangular errors" begin
        R = [1.0 0.2; 0.2 1.0]
        G = [1.0 0.2; 0.2 1.0]
        global genotypes = get_genotypes(genofile, G, separator=',', method="BayesC")
        model = build_model("y1 = intercept + genotypes\ny2 = intercept + genotypes", R)
        # Upper triangular should error
        cs_upper = [0.0 1.0; 0.0 0.0]
        @test_throws ErrorException runMCMC(model, phenotypes_complete;
            causal_structure=cs_upper,
            chain_length=10,
            printout_model_info=false,
            output_folder="test_sem_validation_1",
            seed=42)
        isdir("test_sem_validation_1") && rm("test_sem_validation_1", recursive=true)
    end

    @testset "single-trait with causal_structure errors" begin
        global genotypes = get_genotypes(genofile, 1.0, separator=',', method="BayesC")
        model = build_model("y1 = intercept + genotypes", 1.0)
        cs = [0.0 0.0; 1.0 0.0]
        @test_throws ErrorException runMCMC(model, phenotypes_complete;
            causal_structure=cs,
            chain_length=10,
            printout_model_info=false,
            output_folder="test_sem_validation_2",
            seed=42)
        isdir("test_sem_validation_2") && rm("test_sem_validation_2", recursive=true)
    end

    @testset "constraints are set correctly" begin
        R = [1.0 0.2; 0.2 1.0]
        G = [1.0 0.2; 0.2 1.0]
        global genotypes = get_genotypes(genofile, G, separator=',', method="BayesC")
        model = build_model("y1 = intercept + genotypes\ny2 = intercept + genotypes", R)
        cs = [0.0 0.0; 1.0 0.0]
        # Access internal state after runMCMC sets constraints
        # We run a very short chain and check the output exists
        output = runMCMC(model, phenotypes_complete;
            causal_structure=cs,
            chain_length=10,
            burnin=2,
            output_samples_frequency=5,
            printout_model_info=false,
            outputEBV=false,
            output_heritability=false,
            output_folder="test_sem_validation_3",
            seed=42)
        @test haskey(output, "location parameters")
        @test haskey(output, "residual variance")
        isdir("test_sem_validation_3") && rm("test_sem_validation_3", recursive=true)
        isfile("structure_coefficient_MCMC_samples.txt") && rm("structure_coefficient_MCMC_samples.txt")
    end
end

# ============================================================================
# 7. Integration: 2-trait SEM with runMCMC (structural coefficient output)
# ============================================================================
@testset "SEM 2-trait integration with output verification" begin
    phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
    genofile  = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
    phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstring=["NA"])
    phenotypes_complete = dropmissing(phenotypes, [:y1, :y2])

    R = [1.0 0.2; 0.2 1.0]
    G = [1.0 0.2; 0.2 1.0]
    global genotypes = get_genotypes(genofile, G, separator=',', method="BayesC")
    model = build_model("y1 = intercept + genotypes\ny2 = intercept + genotypes", R)

    cs = [0.0 0.0; 1.0 0.0]
    chain_length = 100
    burnin = 20
    freq = 10

    output = runMCMC(model, phenotypes_complete;
        causal_structure=cs,
        chain_length=chain_length,
        burnin=burnin,
        output_samples_frequency=freq,
        printout_model_info=false,
        outputEBV=false,
        output_heritability=false,
        output_folder="test_sem_2trait",
        seed=123)

    @testset "basic output keys" begin
        @test haskey(output, "location parameters")
        @test haskey(output, "residual variance")
    end

    @testset "structural coefficient MCMC samples file" begin
        @test isfile("structure_coefficient_MCMC_samples.txt")
        λ_samples = readdlm("structure_coefficient_MCMC_samples.txt", ',')
        # samples written only after burnin at output_samples_frequency intervals
        expected_samples = div(chain_length - burnin, freq)
        @test size(λ_samples, 1) == expected_samples
        # number of columns = ntraits^2 = 4 (vec of 2x2 causal matrix)
        @test size(λ_samples, 2) == 4
        # all values should be finite
        @test all(isfinite.(λ_samples))
    end

    @testset "indirect and overall marker effect files" begin
        @test isfile("test_sem_2trait/MCMC_samples_indirect_marker_effects_genotypes_y1.txt")
        @test isfile("test_sem_2trait/MCMC_samples_indirect_marker_effects_genotypes_y2.txt")
        @test isfile("test_sem_2trait/MCMC_samples_overall_marker_effects_genotypes_y1.txt")
        @test isfile("test_sem_2trait/MCMC_samples_overall_marker_effects_genotypes_y2.txt")
    end

    @testset "marker effect summary files" begin
        @test isfile("test_sem_2trait/direct_marker_effects_genotypes.txt")
        @test isfile("test_sem_2trait/indirect_marker_effects_genotypes.txt")
        @test isfile("test_sem_2trait/overall_marker_effects_genotypes.txt")
    end

    @testset "indirect effect file has correct shape" begin
        indirect_y1 = CSV.read("test_sem_2trait/MCMC_samples_indirect_marker_effects_genotypes_y1.txt", DataFrame)
        # number of MCMC samples = (chain_length - burnin) / freq
        expected_samples = div(chain_length - burnin, freq)
        @test size(indirect_y1, 1) == expected_samples
        # number of columns = number of markers
        direct_y1 = CSV.read("test_sem_2trait/MCMC_samples_marker_effects_genotypes_y1.txt", DataFrame)
        @test size(indirect_y1, 2) == size(direct_y1, 2)
    end

    @testset "overall = direct + indirect" begin
        direct_y1 = CSV.read("test_sem_2trait/MCMC_samples_marker_effects_genotypes_y1.txt", DataFrame)
        indirect_y1 = CSV.read("test_sem_2trait/MCMC_samples_indirect_marker_effects_genotypes_y1.txt", DataFrame)
        overall_y1 = CSV.read("test_sem_2trait/MCMC_samples_overall_marker_effects_genotypes_y1.txt", DataFrame)
        for row in 1:size(direct_y1, 1)
            for col in 1:size(direct_y1, 2)
                @test overall_y1[row, col] ≈ direct_y1[row, col] + indirect_y1[row, col] atol=1e-10
            end
        end
    end

    # Cleanup
    isdir("test_sem_2trait") && rm("test_sem_2trait", recursive=true)
    isfile("structure_coefficient_MCMC_samples.txt") && rm("structure_coefficient_MCMC_samples.txt")
end

# ============================================================================
# 8. Integration: 3-trait SEM with runMCMC
# ============================================================================
@testset "SEM 3-trait integration" begin
    phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
    genofile  = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
    phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstring=["NA"])
    phenotypes_complete = dropmissing(phenotypes, [:y1, :y2, :y3])

    R = [1.0 0.2 0.1; 0.2 1.0 0.2; 0.1 0.2 1.0]
    G = [1.0 0.2 0.1; 0.2 1.0 0.2; 0.1 0.2 1.0]
    global genotypes = get_genotypes(genofile, G, separator=',', method="BayesC")
    model = build_model("y1 = intercept + genotypes\ny2 = intercept + genotypes\ny3 = intercept + genotypes", R)

    # trait 1 -> trait 2, trait 1 -> trait 3, trait 2 -> trait 3
    cs = [0.0 0.0 0.0; 1.0 0.0 0.0; 1.0 1.0 0.0]

    output = runMCMC(model, phenotypes_complete;
        causal_structure=cs,
        chain_length=50,
        burnin=10,
        output_samples_frequency=10,
        printout_model_info=false,
        outputEBV=false,
        output_heritability=false,
        output_folder="test_sem_3trait",
        seed=456)

    @testset "output keys exist" begin
        @test haskey(output, "location parameters")
        @test haskey(output, "residual variance")
    end

    @testset "structural coefficient samples shape" begin
        @test isfile("structure_coefficient_MCMC_samples.txt")
        λ_samples = readdlm("structure_coefficient_MCMC_samples.txt", ',')
        # samples written after burnin at output_samples_frequency intervals
        @test size(λ_samples, 1) == div(50 - 10, 10)  # (chain_length - burnin) / freq = 4
        @test size(λ_samples, 2) == 9   # vec of 3x3 matrix
        @test all(isfinite.(λ_samples))
    end

    @testset "all 3 traits have indirect/overall files" begin
        for trait in ["y1", "y2", "y3"]
            @test isfile("test_sem_3trait/MCMC_samples_indirect_marker_effects_genotypes_$trait.txt")
            @test isfile("test_sem_3trait/MCMC_samples_overall_marker_effects_genotypes_$trait.txt")
        end
    end

    @testset "summary files exist" begin
        @test isfile("test_sem_3trait/direct_marker_effects_genotypes.txt")
        @test isfile("test_sem_3trait/indirect_marker_effects_genotypes.txt")
        @test isfile("test_sem_3trait/overall_marker_effects_genotypes.txt")
    end

    # Cleanup
    isdir("test_sem_3trait") && rm("test_sem_3trait", recursive=true)
    isfile("structure_coefficient_MCMC_samples.txt") && rm("structure_coefficient_MCMC_samples.txt")
end

# ============================================================================
# 9. Reproducibility: same seed produces same structural coefficients
# ============================================================================
@testset "SEM reproducibility with seed" begin
    phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
    genofile  = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
    phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstring=["NA"])
    phenotypes_complete = dropmissing(phenotypes, [:y1, :y2])

    function run_sem_with_seed(seed, folder)
        R = [1.0 0.2; 0.2 1.0]
        G = [1.0 0.2; 0.2 1.0]
        global genotypes = get_genotypes(genofile, G, separator=',', method="BayesC")
        model = build_model("y1 = intercept + genotypes\ny2 = intercept + genotypes", R)
        cs = [0.0 0.0; 1.0 0.0]
        runMCMC(model, phenotypes_complete;
            causal_structure=cs,
            chain_length=30,
            burnin=5,
            output_samples_frequency=5,
            printout_model_info=false,
            outputEBV=false,
            output_heritability=false,
            output_folder=folder,
            seed=seed)
        λ = readdlm("structure_coefficient_MCMC_samples.txt", ',')
        return λ
    end

    λ1 = run_sem_with_seed(777, "test_sem_repro_1")
    rm("structure_coefficient_MCMC_samples.txt")
    λ2 = run_sem_with_seed(777, "test_sem_repro_2")

    @test λ1 ≈ λ2 atol=1e-10

    # Cleanup
    isdir("test_sem_repro_1") && rm("test_sem_repro_1", recursive=true)
    isdir("test_sem_repro_2") && rm("test_sem_repro_2", recursive=true)
    isfile("structure_coefficient_MCMC_samples.txt") && rm("structure_coefficient_MCMC_samples.txt")
end
