using Test, JWAS, DataFrames, CSV, JWAS.Datasets

@testset "SEM causal_structure regression (issue #162)" begin
    phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
    genofile = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
    phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstring=["NA"])
    phenotypes = dropmissing(phenotypes, [:y1, :y2])

    R = [1.0 0.2; 0.2 1.0]
    G = [1.0 0.2; 0.2 1.0]
    global genotypes = get_genotypes(genofile, G, separator=',', method="BayesC")
    model = build_model("y1 = intercept + genotypes\ny2 = intercept + genotypes", R)

    causal_structure = [0.0 0.0; 1.0 0.0]
    output = runMCMC(model, phenotypes;
                     causal_structure=causal_structure,
                     chain_length=20,
                     burnin=5,
                     output_samples_frequency=5,
                     printout_model_info=false,
                     outputEBV=false,
                     output_heritability=false,
                     output_folder="test_sem_issue162",
                     seed=123)

    @test haskey(output, "location parameters")
    @test haskey(output, "residual variance")

    rm("test_sem_issue162", recursive=true)
    if isfile("structure_coefficient_MCMC_samples.txt")
        rm("structure_coefficient_MCMC_samples.txt")
    end
end
