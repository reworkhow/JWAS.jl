using Test, JWAS, DataFrames, CSV, JWAS.Datasets, LinearAlgebra, Random, Distributions

Base.@kwdef mutable struct DummyVarianceState
    df::Float64
    scale::Float64
    val::Float64 = 0.0
end

Base.@kwdef mutable struct DummyBayesRState
    method::String = "BayesR"
    ntraits::Int = 1
    α::Vector{Vector{Float64}}
    δ::Vector{Vector{Int}}
    G::DummyVarianceState
end

phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
genofile = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstring=["NA"])
const BAYESR_BLOCK_DISPATCH_COUNTER = Ref(0)
const BAYESR_BLOCK_DISPATCH_INSTALLED = Ref(false)

function bayesr_run_error(model, phenotypes; kwargs...)
    outdir = tempname()
    err = try
        runMCMC(model, phenotypes;
                chain_length=10,
                burnin=0,
                output_samples_frequency=5,
                output_folder=outdir,
                seed=123,
                printout_model_info=false,
                outputEBV=false,
                output_heritability=false,
                kwargs...)
        nothing
    catch exc
        exc
    end
    isdir(outdir) && rm(outdir, recursive=true)
    return err
end

function install_bayesr_block_dispatch_probe(genotype_type)
    if BAYESR_BLOCK_DISPATCH_INSTALLED[]
        return nothing
    end
    @eval function JWAS.BayesR_block!(genotypes::$genotype_type, ycorr, vare, Rinv, iter::Integer, burnin::Integer)
        Main.BAYESR_BLOCK_DISPATCH_COUNTER[] += 1
        return invoke(JWAS.BayesR_block!, Tuple{Any, Any, Any, Any, Integer, Integer},
                      genotypes, ycorr, vare, Rinv, iter, burnin)
    end
    BAYESR_BLOCK_DISPATCH_INSTALLED[] = true
    return nothing
end

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
        err = bayesr_run_error(model, phenotypes)
        @test err isa ErrorException
        @test occursin("length 4", sprint(showerror, err))
    end

    @testset "BayesR fast_blocks gets past validation" begin
        global geno = get_genotypes(genofile, 1.0, separator=',',
                                    method="BayesR",
                                    Pi=Float64[0.95, 0.03, 0.015, 0.005],
                                    estimatePi=true)
        model = build_model("y1 = intercept + geno", 1.0)
        err = bayesr_run_error(model, phenotypes; fast_blocks=true)
        @test !(err isa ErrorException &&
                occursin("BayesR v1 does not support fast_blocks", sprint(showerror, err)))
    end

    @testset "BayesR fast_blocks=1 means block size 1" begin
        global geno = get_genotypes(genofile, 1.0, separator=',',
                                    method="BayesR",
                                    Pi=Float64[0.95, 0.03, 0.015, 0.005],
                                    estimatePi=false,
                                    estimate_variance=false)
        model = build_model("y1 = intercept + geno", 1.0)
        outdir = tempname()
        runMCMC(model, phenotypes;
                chain_length=10,
                burnin=0,
                output_samples_frequency=5,
                output_folder=outdir,
                seed=123,
                printout_model_info=false,
                outputEBV=false,
                output_heritability=false,
                fast_blocks=1)
        @test model.MCMCinfo.fast_blocks == collect(1:1:geno.nMarkers)
        isdir(outdir) && rm(outdir, recursive=true)
    end

    @testset "BayesR copies caller Pi input" begin
        start_pi = Float64[0.95, 0.03, 0.015, 0.005]
        original_pi = copy(start_pi)
        global geno = get_genotypes(genofile, 1.0, separator=',',
                                    method="BayesR",
                                    Pi=start_pi,
                                    estimatePi=true)
        @test geno.π !== start_pi

        model = build_model("y1 = intercept + geno", 1.0)
        outdir = tempname()
        runMCMC(model, phenotypes;
                chain_length=10,
                burnin=0,
                output_samples_frequency=5,
                output_folder=outdir,
                seed=123,
                printout_model_info=false,
                outputEBV=false,
                output_heritability=false)
        @test start_pi == original_pi
        isdir(outdir) && rm(outdir, recursive=true)
    end

    @testset "BayesR still rejects stream, multi-trait, annotations, and RRM" begin
        @testset "stream is rejected" begin
            mktempdir() do tmpdir
                cd(tmpdir) do
                    open("geno.csv", "w") do io
                        println(io, "ID,m1,m2,m3")
                        println(io, "a1,0,1,2")
                        println(io, "a2,1,0,1")
                        println(io, "a3,2,1,0")
                        println(io, "a4,0,2,1")
                    end
                    prefix = JWAS.prepare_streaming_genotypes("geno.csv";
                                                              separator=',',
                                                              header=true,
                                                              quality_control=false,
                                                              center=true)
                    local stream_pheno = DataFrame(ID=["a1", "a2", "a3", "a4"],
                                                   y1=Float32[1.0, -0.2, 0.5, -0.8])
                    global geno = get_genotypes(prefix, 1.0;
                                                method="BayesR",
                                                storage=:stream,
                                                Pi=Float64[0.95, 0.03, 0.015, 0.005],
                                                estimatePi=true)
                    model = build_model("y1 = intercept + geno", 1.0)
                    err = bayesr_run_error(model, stream_pheno; memory_guard=:off)
                    @test err isa ErrorException
                    @test occursin("storage=:dense", sprint(showerror, err))
                end
            end
        end

        @testset "multi-trait is rejected" begin
            phenotypes_mt = DataFrame(ID=copy(phenotypes.ID),
                                      y1=copy(phenotypes.y1),
                                      y2=coalesce.(phenotypes.y1, 0.0))
            global geno = get_genotypes(genofile,
                                        [1.0 0.0; 0.0 1.0],
                                        separator=',',
                                        method="BayesR",
                                        Pi=Float64[0.95, 0.03, 0.015, 0.005],
                                        estimatePi=true)
            model = build_model("y1 = intercept + geno\ny2 = intercept + geno",
                                [1.0 0.0; 0.0 1.0])
            err = bayesr_run_error(model, phenotypes_mt)
            @test err isa ErrorException
            @test occursin("single-trait", sprint(showerror, err))
        end

        @testset "annotations are rejected" begin
            annotations = rand(Float64, 5, 1)
            err = try
                get_genotypes(genofile, 1.0,
                              separator=',',
                              method="BayesR",
                              Pi=Float64[0.95, 0.03, 0.015, 0.005],
                              estimatePi=true,
                              annotations=annotations,
                              quality_control=false)
                nothing
            catch exc
                exc
            end
            @test err isa ErrorException
            @test occursin("annotations", sprint(showerror, err))
        end

        @testset "RRM is rejected" begin
            global geno = get_genotypes(genofile, 1.0, separator=',',
                                        method="BayesR",
                                        Pi=Float64[0.95, 0.03, 0.015, 0.005],
                                        estimatePi=true)
            model = build_model("y1 = intercept + geno", 1.0)
            err = bayesr_run_error(model, phenotypes; RRM=ones(Float64, 1, 1))
            @test err isa ErrorException
            @test occursin("RRM", sprint(showerror, err))
        end
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

@testset "BayesR block sampler" begin
    X = Float64[
        0 2
        1 1
        2 0
        1 1
    ]
    XArray = [X]
    xpRinvx = Float64[dot(view(X, :, 1), view(X, :, 1)),
                       dot(view(X, :, 2), view(X, :, 2))]
    XpRinvX = [Matrix(X' * X)]
    yCorr = Float64[0.8, -0.1, 0.3, 0.5]
    α = zeros(Float64, 2)
    δ = ones(Int, 2)
    π = Float64[0.95, 0.03, 0.015, 0.005]
    gamma = Float64[0.0, 0.01, 0.1, 1.0]

    JWAS.BayesR_block!(XArray, xpRinvx, X, XpRinvX,
                       yCorr, α, δ, 1.0, 0.2, π, gamma, ones(Float64, length(yCorr)))

    @test all(1 .<= δ .<= 4)
    @test length(α) == 2
end

@testset "BayesR fast-block repetition schedule" begin
    @test JWAS.bayesr_block_nreps(1, 10, 7) == 1
    @test JWAS.bayesr_block_nreps(10, 10, 7) == 1
    @test JWAS.bayesr_block_nreps(11, 10, 7) == 7
    @test JWAS.bayesr_block_nreps(25, 0, 7) == 7
    @test JWAS.bayesr_block_nreps(3, 8, 1) == 1
end

@testset "BayesR variance sufficient statistics" begin
    α = Float64[0.0, 0.4, -0.3, 0.1, 0.0]
    δ = Int[1, 2, 4, 3, 1]
    gamma = Float64[0.0, 0.01, 0.1, 1.0]

    ssq, nnz = JWAS.bayesr_sigma_sufficient_statistics(α, δ, gamma)

    expected_ssq = α[2]^2 / gamma[2] + α[3]^2 / gamma[4] + α[4]^2 / gamma[3]
    @test ssq ≈ expected_ssq
    @test nnz == 3
end

@testset "BayesR variance update matches direct ssq formula" begin
    geno = DummyBayesRState(
        α=[Float64[0.0, 0.4, -0.3, 0.1, 0.0]],
        δ=[Int[1, 2, 4, 3, 1]],
        G=DummyVarianceState(df=4.0, scale=0.2),
    )
    gamma = Float64[0.0, 0.01, 0.1, 1.0]

    ssq, nnz = JWAS.bayesr_sigma_sufficient_statistics(geno.α[1], geno.δ[1], gamma)
    Random.seed!(1234)
    expected = (ssq + geno.G.df * geno.G.scale) / rand(Chisq(nnz + geno.G.df))

    Random.seed!(1234)
    JWAS.sample_marker_effect_variance(geno)

    @test geno.G.val ≈ expected
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
    @test occursin("expected class counts", lowercase(printed))
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

@testset "BayesR fast_blocks dispatch" begin
    global geno = get_genotypes(genofile, 1.0, separator=',',
                                method="BayesR",
                                Pi=Float64[0.95, 0.03, 0.015, 0.005],
                                estimatePi=false,
                                estimate_variance=true)
    install_bayesr_block_dispatch_probe(typeof(geno))
    BAYESR_BLOCK_DISPATCH_COUNTER[] = 0
    model = build_model("y1 = intercept + geno", 1.0)
    outdir = tempname()

    _ = runMCMC(model, phenotypes,
                chain_length=20,
                burnin=5,
                output_samples_frequency=5,
                output_folder=outdir,
                seed=777,
                printout_model_info=false,
                outputEBV=false,
                output_heritability=false,
                fast_blocks=true)

    @test BAYESR_BLOCK_DISPATCH_COUNTER[] > 0
    isdir(outdir) && rm(outdir, recursive=true)
end
