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
end
