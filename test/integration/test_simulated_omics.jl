# Integration Tests Using Simulated Omics Data
# Dataset: 3,534 genotyped animals, 1,000 SNPs, 6,473 pedigree animals
# Known true effects enable numerical assertions on prediction accuracy
#
# Data structure (from simulation):
#   - trait1 = group + litter + direct_SNP_genetic + omics_contribution + residual
#   - Total genetic h² ≈ 0.25 (20% direct, 80% indirect via omics)
#   - 10 omics layers driven by 10 consecutive SNP blocks

using Test, CSV, DataFrames, JWAS, JWAS.Datasets, Statistics, Random

# ============================================================================
# DATA LOADING
# ============================================================================

phenofile = Datasets.dataset("phenotypes_sim.txt", dataset_name="simulated_omics")
genofile  = Datasets.dataset("genotypes_1000snps.txt", dataset_name="simulated_omics")
pedfile   = Datasets.dataset("pedigree.txt", dataset_name="simulated_omics")

phenotypes = CSV.read(phenofile, DataFrame, delim=',', header=true, missingstrings=["NA"])
pedigree   = get_pedigree(pedfile, separator=",", header=true)

# Convert group and litter to string for categorical treatment
phenotypes[!, :group]  = string.(phenotypes[!, :group])
phenotypes[!, :litter] = string.(phenotypes[!, :litter])

# ============================================================================
# TEST 1: Single-Trait Bayesian Methods (Complete Genomic Data)
# ============================================================================

@testset "Single-trait Bayesian methods" begin
    for test_method in ["BayesA", "BayesB", "BayesC", "RR-BLUP", "BayesL", "GBLUP"]
        @testset "$test_method" begin
            estimatePi = test_method in ["BayesC", "BayesB"]
            G = 1.0
            global geno = get_genotypes(genofile, G, header=true, separator=',',
                                        method=test_method, estimatePi=estimatePi)

            model = build_model("trait1 = intercept + geno", 1.0)

            outdir = "st_omics_$(test_method)"
            output = runMCMC(model, phenotypes,
                             chain_length=500,
                             burnin=100,
                             output_samples_frequency=50,
                             output_folder=outdir,
                             printout_frequency=500,
                             seed=314)

            # Output structure checks
            @test haskey(output, "location parameters")
            @test haskey(output, "residual variance")
            @test haskey(output, "marker effects geno")
            @test size(output["location parameters"], 1) > 0

            # Residual variance should be positive
            @test output["residual variance"][1, :Estimate] > 0.0

            # EBV prediction accuracy against true genetic_total
            if haskey(output, "EBV_trait1")
                results = innerjoin(output["EBV_trait1"], phenotypes, on=:ID)
                accuracy = cor(Float64.(results[!, :EBV]), Float64.(results[!, :genetic_total]))
                @test accuracy > 0.1  # Low bar: only 20% of genetic variance is direct
                println("  $test_method prediction accuracy (vs genetic_total): $(round(accuracy, digits=3))")
            end

            # Output files should exist
            @test isdir(outdir)
            @test isfile(joinpath(outdir, "location_parameters.txt"))
            @test isfile(joinpath(outdir, "residual_variance.txt"))

            # Cleanup
            rm(outdir, recursive=true, force=true)
        end
    end
end

# ============================================================================
# TEST 2: Single-Step Analysis (Incomplete Genomic Data)
# ============================================================================

# For single-step, we need phenotyped individuals WITHOUT genotypes.
# The pedigree has 6,473 animals but only 3,534 are genotyped.
# Create a phenotype set that uses only a subset of genotyped animals,
# plus some non-genotyped pedigree animals with simulated phenotypes.
@testset "Single-step analysis" begin
    # Build phenotype data for single-step: keep half the genotyped animals
    # and add non-genotyped pedigree animals with random phenotypes
    pedfile_df = CSV.read(pedfile, DataFrame, delim=',', header=true)
    geno_ids = Set(string.(phenotypes[!, :ID]))
    non_geno_ids = setdiff(Set(string.(pedfile_df[!, :ID])), geno_ids)

    # Take first 500 non-genotyped pedigree animals, give them phenotypes
    extra_ids = collect(non_geno_ids)[1:min(500, length(non_geno_ids))]
    Random.seed!(42)
    extra_pheno = DataFrame(
        ID = extra_ids,
        trait1 = randn(length(extra_ids)),
        genetic_total = zeros(length(extra_ids))  # unknown true BV
    )
    phenotypes_ss = vcat(
        select(phenotypes, [:ID, :trait1, :genetic_total]),
        extra_pheno
    )
    # Convert ID to string to match genotype IDs
    phenotypes_ss[!, :ID] = string.(phenotypes_ss[!, :ID])

    for test_method in ["BayesC", "RR-BLUP", "BayesA"]
        @testset "SS-$test_method" begin
            estimatePi = test_method == "BayesC"
            G = 1.0
            global geno = get_genotypes(genofile, G, header=true, separator=',',
                                        method=test_method, estimatePi=estimatePi)

            model = build_model("trait1 = intercept + geno", 1.0)

            outdir = "ss_omics_$(test_method)"
            output = runMCMC(model, phenotypes_ss,
                             chain_length=500,
                             burnin=100,
                             output_samples_frequency=50,
                             output_folder=outdir,
                             printout_frequency=500,
                             single_step_analysis=true,
                             pedigree=pedigree,
                             seed=314)

            # Output structure checks
            @test haskey(output, "location parameters")
            @test haskey(output, "residual variance")
            @test output["residual variance"][1, :Estimate] > 0.0

            # EBV accuracy (only for genotyped animals with known true BV)
            if haskey(output, "EBV_trait1")
                results = innerjoin(output["EBV_trait1"], phenotypes, on=:ID)
                ind_id = findall(x -> !ismissing(x), results[!, :genetic_total])
                accuracy = cor(Float64.(results[ind_id, :EBV]), Float64.(results[ind_id, :genetic_total]))
                @test accuracy > 0.0
                println("  SS-$test_method prediction accuracy: $(round(accuracy, digits=3))")
            end

            # Cleanup
            rm(outdir, recursive=true, force=true)
        end
    end
end

# ============================================================================
# TEST 3: Multi-Trait Analysis
# ============================================================================

@testset "Multi-trait analysis" begin
    G = [1.0 0.5; 0.5 1.0]
    global geno = get_genotypes(genofile, G, header=true, separator=',',
                                method="BayesC", estimatePi=true)

    R = [1.0 0.5; 0.5 1.0]
    model = build_model("trait1 = intercept + geno\nomic1 = intercept + geno", R)

    outdir = "mt_omics"
    output = runMCMC(model, phenotypes,
                     chain_length=500,
                     burnin=100,
                     output_samples_frequency=50,
                     output_folder=outdir,
                     printout_frequency=500,
                     seed=314)

    # Output structure checks
    @test haskey(output, "location parameters")
    @test haskey(output, "residual variance")
    @test size(output["location parameters"], 1) > 0

    # Should have EBVs for both traits
    @test haskey(output, "EBV_trait1")
    @test haskey(output, "EBV_omic1")

    # Residual variance should be a 2x2 structure
    @test size(output["residual variance"], 1) >= 2

    # Cleanup
    rm(outdir, recursive=true, force=true)
end

# ============================================================================
# TEST 4: Model with Random Effects (Group/Litter)
# ============================================================================

@testset "Model with random effects" begin
    G_marker = 1.0
    global geno = get_genotypes(genofile, G_marker, header=true, separator=',',
                                method="BayesC", estimatePi=true)

    model = build_model("trait1 = intercept + group + litter + geno", 1.0)
    set_random(model, "litter", 1.0)

    outdir = "re_omics"
    output = runMCMC(model, phenotypes,
                     chain_length=500,
                     burnin=100,
                     output_samples_frequency=50,
                     output_folder=outdir,
                     printout_frequency=500,
                     seed=314)

    # Output structure checks
    @test haskey(output, "location parameters")
    @test haskey(output, "residual variance")
    @test output["residual variance"][1, :Estimate] > 0.0
    @test size(output["location parameters"], 1) > 0

    # EBV accuracy — model with group+litter should do at least as well
    if haskey(output, "EBV_trait1")
        results = innerjoin(output["EBV_trait1"], phenotypes, on=:ID)
        accuracy = cor(Float64.(results[!, :EBV]), Float64.(results[!, :genetic_total]))
        @test accuracy > 0.1
        println("  Random effects model accuracy: $(round(accuracy, digits=3))")
    end

    # Cleanup
    rm(outdir, recursive=true, force=true)
end

# ============================================================================
# TEST 5: GWAS on Simulated Data
# ============================================================================

@testset "GWAS with simulated data" begin
    G = 1.0
    global geno = get_genotypes(genofile, G, header=true, separator=',',
                                method="BayesC", estimatePi=true)

    model = build_model("trait1 = intercept + geno", 1.0)

    outdir = "gwas_omics"
    output = runMCMC(model, phenotypes,
                     chain_length=500,
                     burnin=100,
                     output_samples_frequency=50,
                     output_folder=outdir,
                     printout_frequency=500,
                     seed=314)

    # Run GWAS on marker effect samples
    gwas_result = GWAS(joinpath(outdir, "MCMC_samples_marker_effects_geno_trait1.txt"))

    @test size(gwas_result, 1) > 0
    @test "marker_ID" in names(gwas_result)
    @test "modelfrequency" in names(gwas_result)
    @test all(0 .<= gwas_result.modelfrequency .<= 1)

    # Cleanup
    rm(outdir, recursive=true, force=true)
end

# ============================================================================
# TEST 6: Reproducibility with Seed
# ============================================================================

@testset "Reproducibility with seed" begin
    G = 1.0
    global geno = get_genotypes(genofile, G, header=true, separator=',', method="RR-BLUP")
    model1 = build_model("trait1 = intercept + geno", 1.0)
    out1 = runMCMC(model1, phenotypes, chain_length=200,
                   output_folder="repro1", printout_frequency=500, seed=999)

    global geno = get_genotypes(genofile, G, header=true, separator=',', method="RR-BLUP")
    model2 = build_model("trait1 = intercept + geno", 1.0)
    out2 = runMCMC(model2, phenotypes, chain_length=200,
                   output_folder="repro2", printout_frequency=500, seed=999)

    @test out1["residual variance"][1, :Estimate] ≈
          out2["residual variance"][1, :Estimate] atol=1e-10

    # EBVs should be identical
    if haskey(out1, "EBV_trait1") && haskey(out2, "EBV_trait1")
        ebv1 = sort(out1["EBV_trait1"], :ID)
        ebv2 = sort(out2["EBV_trait1"], :ID)
        @test ebv1[!, :EBV] ≈ ebv2[!, :EBV] atol=1e-10
    end

    # Cleanup
    rm("repro1", recursive=true, force=true)
    rm("repro2", recursive=true, force=true)
end
