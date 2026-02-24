# Unit tests for advanced coverage: multi-trait BayesB, constraint modes,
# missing phenotypes, heterogeneous residuals, showMME, Pi estimation
using Test, JWAS, DataFrames, CSV, JWAS.Datasets

phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
genofile = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
pedfile = Datasets.dataset("pedigree.txt", dataset_name="demo_7animals")
phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstrings=["NA"])

# ========================================================================
# Multi-trait BayesB (vector variances in multi-trait)
# ========================================================================
@testset "Multi-trait BayesB" begin
    R = [1.0 0.5; 0.5 1.0]
    G = [1.0 0.5; 0.5 1.0]
    global geno = get_genotypes(genofile, G, separator=',', method="BayesB")
    model = build_model("y1 = intercept + geno\ny2 = intercept + geno", R)

    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    output_folder="test_mt_bayesb",
                    seed=123)

    @test haskey(output, "marker effects geno")
    @test haskey(output, "location parameters")
    @test haskey(output, "residual variance")
    rm("test_mt_bayesb", recursive=true)
end

# ========================================================================
# Multi-trait BayesL (Lasso)
# ========================================================================
@testset "Multi-trait BayesL" begin
    R = [1.0 0.5; 0.5 1.0]
    G = [1.0 0.5; 0.5 1.0]
    global geno = get_genotypes(genofile, G, separator=',', method="BayesL")
    model = build_model("y1 = intercept + geno\ny2 = intercept + geno", R)

    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    output_folder="test_mt_bayesl",
                    seed=123)

    @test haskey(output, "marker effects geno")
    @test haskey(output, "location parameters")
    rm("test_mt_bayesl", recursive=true)
end

# ========================================================================
# Multi-trait GBLUP
# ========================================================================
@testset "Multi-trait GBLUP" begin
    R = [1.0 0.5; 0.5 1.0]
    G = [1.0 0.5; 0.5 1.0]
    global geno = get_genotypes(genofile, G, separator=',', method="GBLUP")
    model = build_model("y1 = intercept + geno\ny2 = intercept + geno", R)

    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    output_folder="test_mt_gblup",
                    seed=123)

    @test haskey(output, "location parameters")
    @test haskey(output, "residual variance")
    rm("test_mt_gblup", recursive=true)
end

# ========================================================================
# Constraint modes (G_constraint for mega-trait analysis)
# ========================================================================
@testset "Multi-trait BayesC with G constraint (mega-trait)" begin
    R = [1.0 0.0; 0.0 1.0]
    G = [1.0 0.0; 0.0 1.0]
    global geno = get_genotypes(genofile, G, separator=',', method="BayesC",
                                constraint=true)
    model = build_model("y1 = intercept + geno\ny2 = intercept + geno", R)

    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    output_folder="test_constraint",
                    seed=123)

    @test haskey(output, "location parameters")
    @test haskey(output, "residual variance")
    rm("test_constraint", recursive=true)
end

# ========================================================================
# Multi-trait with missing phenotypes
# ========================================================================
@testset "Multi-trait with missing phenotypes" begin
    R = [1.0 0.5; 0.5 1.0]
    G = [1.0 0.5; 0.5 1.0]
    global geno = get_genotypes(genofile, G, separator=',', method="RR-BLUP")
    model = build_model("y1 = intercept + geno\ny2 = intercept + geno", R)

    # Phenotype data already has missing values (NA)
    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    missing_phenotypes=true,
                    output_folder="test_missing_pheno",
                    seed=123)

    @test haskey(output, "location parameters")
    @test haskey(output, "residual variance")
    rm("test_missing_pheno", recursive=true)
end

# ========================================================================
# Pi estimation in multi-trait BayesC
# ========================================================================
@testset "Multi-trait BayesC with Pi estimation" begin
    R = [1.0 0.5; 0.5 1.0]
    G = [1.0 0.5; 0.5 1.0]
    global geno = get_genotypes(genofile, G, separator=',', method="BayesC",
                                estimatePi=true)
    model = build_model("y1 = intercept + geno\ny2 = intercept + geno", R)

    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    output_folder="test_pi_est",
                    seed=123)

    @test haskey(output, "pi_geno")
    @test haskey(output, "marker effects geno")
    rm("test_pi_est", recursive=true)
end

# ========================================================================
# Single-trait with estimate_scale
# ========================================================================
@testset "Single-trait BayesC with estimate_scale" begin
    global geno = get_genotypes(genofile, 1.0, separator=',', method="BayesC",
                                estimate_scale=true)
    model = build_model("y1 = intercept + geno", 1.0)

    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    output_folder="test_est_scale",
                    seed=123)

    @test haskey(output, "location parameters")
    rm("test_est_scale", recursive=true)
end

# ========================================================================
# Multi-trait with EBV output
# ========================================================================
@testset "Multi-trait EBV output" begin
    R = [1.0 0.5; 0.5 1.0]
    G = [1.0 0.5; 0.5 1.0]
    global geno = get_genotypes(genofile, G, separator=',', method="BayesC")
    model = build_model("y1 = intercept + geno\ny2 = intercept + geno", R)
    outputEBV(model, geno.obsID)

    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    outputEBV=true,
                    output_heritability=true,
                    output_folder="test_mt_ebv",
                    seed=123)

    @test haskey(output, "EBV_y1")
    @test haskey(output, "EBV_y2")
    @test haskey(output, "genetic_variance")
    @test haskey(output, "heritability")
    rm("test_mt_ebv", recursive=true)
end

# ========================================================================
# Heterogeneous residuals (weighted analysis)
# ========================================================================
@testset "Single-trait with weights (heterogeneous residuals)" begin
    global geno = get_genotypes(genofile, 1.0, separator=',', method="BayesC")
    model = build_model("y1 = intercept + geno", 1.0)

    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    heterogeneous_residuals=true,
                    output_folder="test_het_res",
                    seed=123)

    @test haskey(output, "location parameters")
    rm("test_het_res", recursive=true)
end

# ========================================================================
# Double precision MCMC
# ========================================================================
@testset "Double precision MCMC" begin
    global geno = get_genotypes(genofile, 1.0, separator=',', method="BayesC",
                                double_precision=true)
    model = build_model("y1 = intercept + geno", 1.0)

    output = runMCMC(model, phenotypes,
                    chain_length=50,
                    output_samples_frequency=10,
                    double_precision=true,
                    output_folder="test_double_prec",
                    seed=123)

    @test haskey(output, "location parameters")
    rm("test_double_prec", recursive=true)
end
