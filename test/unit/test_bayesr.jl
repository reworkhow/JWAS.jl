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
        outdir = tempname()
        err = try
            runMCMC(model, phenotypes,
                    chain_length=10,
                    burnin=0,
                    output_samples_frequency=5,
                    output_folder=outdir,
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
        isdir(outdir) && rm(outdir, recursive=true)
    end

    @testset "BayesR rejects fast_blocks" begin
        global geno = get_genotypes(genofile, 1.0, separator=',',
                                    method="BayesR",
                                    Pi=Float64[0.95, 0.03, 0.015, 0.005],
                                    estimatePi=true)
        model = build_model("y1 = intercept + geno", 1.0)
        outdir = tempname()
        err = try
            runMCMC(model, phenotypes,
                    chain_length=10,
                    burnin=0,
                    output_samples_frequency=5,
                    output_folder=outdir,
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
        isdir(outdir) && rm(outdir, recursive=true)
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

@testset "BayesR single-trait run" begin
    global geno = get_genotypes(genofile, 1.0, separator=',',
                                method="BayesR",
                                Pi=Float64[0.95, 0.03, 0.015, 0.005],
                                estimatePi=false,
                                estimate_variance=true)
    model = build_model("y1 = intercept + geno", 1.0)
    outdir = tempname()

    outpath, io = mktemp()
    close(io)
    output_ref = Ref{Any}()
    open(outpath, "w") do outio
        redirect_stdout(outio) do
            output_ref[] = runMCMC(model, phenotypes,
                                   chain_length=20,
                                   burnin=5,
                                   output_samples_frequency=5,
                                   output_folder=outdir,
                                   seed=123,
                                   printout_model_info=true,
                                   outputEBV=false,
                                   output_heritability=false,
                                   fast_blocks=false)
        end
    end
    printed = read(outpath, String)
    rm(outpath, force=true)
    output = output_ref[]

    @test occursin("BayesR", printed)
    @test occursin("starting pi", lowercase(printed))
    @test haskey(output, "marker effects geno")
    @test haskey(output, "location parameters")
    isdir(outdir) && rm(outdir, recursive=true)
end

@testset "BayesR estimatePi output" begin
    global geno = get_genotypes(genofile, 1.0, separator=',',
                                method="BayesR",
                                Pi=Float64[0.95, 0.03, 0.015, 0.005],
                                estimatePi=true,
                                estimate_variance=true)
    model = build_model("y1 = intercept + geno", 1.0)
    outdir = tempname()
    output = runMCMC(model, phenotypes,
                     chain_length=30,
                     burnin=5,
                     output_samples_frequency=5,
                     output_folder=outdir,
                     seed=321,
                     printout_model_info=false,
                     outputEBV=false,
                     output_heritability=false,
                     fast_blocks=false)

    @test haskey(output, "pi_geno")
    @test nrow(output["pi_geno"]) == 4
    @test isapprox(sum(output["pi_geno"][!, :Estimate]), 1.0; atol=1e-6)
    @test all(0.0 .<= output["marker effects geno"][!, :Model_Frequency] .<= 1.0)
    isdir(outdir) && rm(outdir, recursive=true)
end
