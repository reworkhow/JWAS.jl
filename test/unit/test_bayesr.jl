using Test, JWAS, DataFrames, CSV, JWAS.Datasets, LinearAlgebra

phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
genofile = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstring=["NA"])

@testset "BayesR validation" begin
    @testset "single-trait dense BayesR genotype loads" begin
        global geno = get_genotypes(genofile, 1.0, separator=',',
                                    method="BayesR",
                                    Pi=Float64[0.95, 0.03, 0.015, 0.005],
                                    estimatePi=true)
        model = build_model("y1 = intercept + geno", 1.0)
        @test geno.method == "BayesR"
        @test model.nModels == 1
    end

    @testset "BayesR rejects bad Pi length" begin
        global geno = get_genotypes(genofile, 1.0, separator=',',
                                    method="BayesR",
                                    Pi=Float64[0.95, 0.05, 0.0],
                                    estimatePi=true)
        model = build_model("y1 = intercept + geno", 1.0)
        err = try
            runMCMC(model, phenotypes,
                    chain_length=10,
                    burnin=0,
                    output_samples_frequency=5,
                    output_folder="test_bayesr_bad_pi",
                    seed=123,
                    printout_model_info=false,
                    outputEBV=false,
                    output_heritability=false)
            nothing
        catch exc
            exc
        end
        @test err isa ErrorException
        @test occursin("length 4", sprint(showerror, err))
        isdir("test_bayesr_bad_pi") && rm("test_bayesr_bad_pi", recursive=true)
    end

    @testset "BayesR rejects fast_blocks" begin
        global geno = get_genotypes(genofile, 1.0, separator=',',
                                    method="BayesR",
                                    Pi=Float64[0.95, 0.03, 0.015, 0.005],
                                    estimatePi=true)
        model = build_model("y1 = intercept + geno", 1.0)
        err = try
            runMCMC(model, phenotypes,
                    chain_length=10,
                    burnin=0,
                    output_samples_frequency=5,
                    output_folder="test_bayesr_fastblocks",
                    seed=123,
                    printout_model_info=false,
                    outputEBV=false,
                    output_heritability=false,
                    fast_blocks=true)
            nothing
        catch exc
            exc
        end
        @test err isa ErrorException
        @test occursin("BayesR", sprint(showerror, err))
        isdir("test_bayesr_fastblocks") && rm("test_bayesr_fastblocks", recursive=true)
    end
end

@testset "BayesR dense sampler" begin
    x1 = Float64[0, 1, 2, 1]
    x2 = Float64[2, 1, 0, 1]
    xArray = [x1, x2]
    xRinvArray = [x1, x2]
    xpRinvx = Float64[dot(x1, x1), dot(x2, x2)]
    yCorr = Float64[0.8, -0.1, 0.3, 0.5]
    α = zeros(Float64, 2)
    δ = ones(Int, 2)
    π = Float64[0.95, 0.03, 0.015, 0.005]
    gamma = Float64[0.0, 0.01, 0.1, 1.0]

    JWAS.BayesR!(xArray, xRinvArray, xpRinvx, yCorr, α, δ, 1.0, 0.2, π, gamma)

    @test all(1 .<= δ .<= 4)
    @test length(α) == 2
end
