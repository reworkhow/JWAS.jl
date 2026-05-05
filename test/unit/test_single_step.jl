# Unit tests for single-step analysis regressions
using Test, JWAS, DataFrames, CSV, JWAS.Datasets

@testset "Single-step analysis" begin
    phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
    pedfile = Datasets.dataset("pedigree.txt", dataset_name="demo_7animals")
    genofile = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")

    phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstring=["NA"])
    pedigree = get_pedigree(pedfile, separator=',', header=true)

    @testset "BayesC accepts CSV inline-string IDs" begin
        output_folder = "test_single_step_inline_ids"
        isdir(output_folder) && rm(output_folder, recursive=true, force=true)
        global genotypes = get_genotypes(genofile, separator=',', method="BayesC")
        model = build_model("y1 = intercept + genotypes")

        local output
        try
            output = runMCMC(
                model,
                phenotypes;
                single_step_analysis=true,
                pedigree=pedigree,
                chain_length=2,
                output_samples_frequency=1,
                output_heritability=false,
                printout_model_info=false,
                output_folder=output_folder,
                seed=160,
            )

            @test haskey(output, "EBV_y1")
            @test eltype(model.M[1].genotypes) == Float32
            @test eltype(model.M[1].output_genotypes) == Float32
        finally
            isdir(output_folder) && rm(output_folder, recursive=true, force=true)
        end
    end

    @testset "BayesC accepts multiple genotype categories" begin
        output_folder = "test_single_step_multiclass"
        isdir(output_folder) && rm(output_folder, recursive=true, force=true)
        global geno1 = get_genotypes(genofile, separator=',', method="BayesC")
        global geno2 = get_genotypes(genofile, separator=',', method="BayesC")
        model = build_model("y1 = intercept + geno1 + geno2")

        local output
        try
            output = runMCMC(
                model,
                phenotypes;
                single_step_analysis=true,
                pedigree=pedigree,
                chain_length=2,
                output_samples_frequency=1,
                output_heritability=false,
                printout_model_info=false,
                output_folder=output_folder,
                seed=161,
            )

            @test haskey(output, "EBV_y1")
            @test haskey(output, "marker effects geno1")
            @test haskey(output, "marker effects geno2")
            @test haskey(model.modelTermDict, "y1:ϵ")
            @test haskey(model.modelTermDict, "y1:J")
            @test !haskey(model.modelTermDict, "y1:epsilon_geno1")
            @test !haskey(model.modelTermDict, "y1:epsilon_geno2")
            @test count(effect -> effect.term_array == ["y1:ϵ"], model.rndTrmVec) == 1
            @test JWAS.total_single_step_genetic_variance(model.M) ≈
                  model.M[1].genetic_variance.val + model.M[2].genetic_variance.val
        finally
            isdir(output_folder) && rm(output_folder, recursive=true, force=true)
        end
    end

    @testset "Annotated BayesC and BayesR accept multiple genotype categories" begin
        annotations = [
            0.0 1.0
            1.0 0.0
            1.0 1.0
            0.0 0.0
            0.5 0.5
        ]

        for method in ["BayesC", "BayesR"]
            output_folder = "test_single_step_annotated_" * lowercase(method)
            isdir(output_folder) && rm(output_folder, recursive=true, force=true)
            global annotated_ss_geno1 = get_genotypes(
                genofile,
                1.0;
                separator=',',
                method=method,
                quality_control=false,
                annotations=annotations,
            )
            global annotated_ss_geno2 = get_genotypes(
                genofile,
                1.0;
                separator=',',
                method=method,
                quality_control=false,
                annotations=annotations,
            )
            model = build_model("y1 = intercept + annotated_ss_geno1 + annotated_ss_geno2", 1.0)

            local output
            try
                output = runMCMC(
                    model,
                    phenotypes;
                    single_step_analysis=true,
                    pedigree=pedigree,
                    chain_length=2,
                    output_samples_frequency=1,
                    output_heritability=false,
                    printout_model_info=false,
                    output_folder=output_folder,
                    seed=162,
                )

                @test haskey(output, "EBV_y1")
                @test haskey(output, "marker effects annotated_ss_geno1")
                @test haskey(output, "marker effects annotated_ss_geno2")
                @test haskey(output, "annotation coefficients annotated_ss_geno1")
                @test haskey(output, "annotation coefficients annotated_ss_geno2")
                @test haskey(model.modelTermDict, "y1:ϵ")
                @test haskey(model.modelTermDict, "y1:J")
                @test count(effect -> effect.term_array == ["y1:ϵ"], model.rndTrmVec) == 1

                if method == "BayesR"
                    @test size(model.M[1].annotations.snp_pi, 2) == 4
                    @test size(model.M[2].annotations.snp_pi, 2) == 4
                    @test "Step" in names(output["annotation coefficients annotated_ss_geno1"])
                    @test "Step" in names(output["annotation coefficients annotated_ss_geno2"])
                else
                    @test model.M[1].π isa AbstractVector
                    @test model.M[2].π isa AbstractVector
                    @test length(model.M[1].π) == model.M[1].nMarkers
                    @test length(model.M[2].π) == model.M[2].nMarkers
                end
            finally
                isdir(output_folder) && rm(output_folder, recursive=true, force=true)
            end
        end
    end
end
