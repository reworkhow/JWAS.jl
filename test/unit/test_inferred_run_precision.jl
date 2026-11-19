using Test, JWAS, DataFrames, CSV, Random
if !isdefined(@__MODULE__, :bayesc_memory_case)
    include("bayesc_memory_fixture.jl")
end

function inferred_precision_multiterm(types)
    rng = MersenneTwister(20261118)
    X = types[1].(rand(rng, 0:2, 24, 7))
    Y = types[2].(rand(rng, 0:2, 24, 7))
    global inferred_first = get_genotypes(X, 1.0; method="BayesC", Pi=0.65,
        quality_control=false, double_precision=(types[1] === Float64))
    global inferred_second = get_genotypes(Y, 1.0; method="BayesC", Pi=0.65,
        quality_control=false, double_precision=(types[2] === Float64))
    model = build_model("y = intercept + inferred_first + inferred_second", 1.0)
    phenotypes = DataFrame(ID=string.(1:24), y=randn(rng, 24))
    (; model, phenotypes, parents=(X, Y))
end

function inferred_precision_small_run(model, phenotypes, output_folder; kwargs...)
    runMCMC(model, phenotypes; chain_length=20, burnin=5, output_samples_frequency=5,
        seed=73, printout_model_info=false, outputEBV=false, output_heritability=false,
        output_folder, kwargs...)
end

@testset "Genotype-inferred run precision" begin
    @testset "Automatic and explicitly matched production precision" begin
        for traits in (1, 2), sampler in (traits == 1 ? (:I,) : (:I, :II)),
            T in (Float32, Float64), seed in (73, 123)
            @testset "traits=$traits sampler=$sampler type=$T seed=$seed" begin
                automatic = bayesc_memory_case(; traits, sampler, precision=T, seed,
                    mode="sequential", outer=80, burnin=20, inherit_run_precision=true)
                explicit = bayesc_memory_case(; traits, sampler, precision=T, seed,
                    mode="sequential", outer=80, burnin=20)
                @test automatic.model.MCMCinfo.double_precision === (T === Float64)
                @test eltype(automatic.model.sol) === T
                @test eltype(automatic.model.invweights) === T
                @test all(a -> eltype(a) === T, automatic.g.α)
                @test eltype(automatic.g.genotypes) === T
                @test JWAS.dense_genotype_parent(automatic.g.genotypes) === automatic.input
                @test JWAS.dense_genotype_parent(automatic.g.output_genotypes) === automatic.input
                @test all(block -> block isa SubArray && Base.mightalias(block, automatic.input), automatic.g.MArray)
                aligned = compare_bayesc_memory(automatic, explicit)
                @test aligned.animals.EBV == aligned.animals.EBV_1
                @test aligned.animals.PEV == aligned.animals.PEV_1
                @test aligned.markers.Estimate == aligned.markers.Estimate_1
                @test aligned.markers.SD == aligned.markers.SD_1
                @test aligned.markers.Model_Frequency == aligned.markers.Model_Frequency_1
                if traits == 1 && seed == 73
                    println("Inspected auto/explicit $T animal join: ", first(aligned.animals[:, [:Trait, :ID, :EBV, :EBV_1]], 2))
                    println("Inspected auto/explicit $T marker join: ", first(aligned.markers[:, [:Trait, :Marker_ID, :Estimate, :Estimate_1]], 2))
                end
            end
        end
    end

    @testset "Consistent multiple terms and conventional defaults" begin
        for T in (Float32, Float64)
            fixture = inferred_precision_multiterm((T, T))
            @test JWAS.infer_run_double_precision(fixture.model) === (T === Float64)
            @test all(i -> fixture.model.M[i].genotypes === fixture.parents[i], 1:2)
        end
        model = build_model("y = intercept", 1.0)
        @test JWAS.infer_run_double_precision(model) === false
        ph = DataFrame(ID=string.(1:12), y=sin.(1:12))
        @test inferred_precision_small_run(model, ph, joinpath(mktempdir(), "conventional")) isa AbstractDict
        @test model.MCMCinfo.double_precision === false
        @test eltype(model.sol) === Float32
        empty_model = build_model("y = intercept", 1.0)
        empty_model.M = JWAS.Genotypes[]
        @test JWAS.infer_run_double_precision(empty_model) === false
    end

    @testset "Mixed and unsupported types fail before run setup" begin
        fixture = inferred_precision_multiterm((Float32, Float64))
        initial_mcmcinfo = fixture.model.MCMCinfo
        output_folder = joinpath(mktempdir(), "mixed")
        err = try
            inferred_precision_small_run(fixture.model, fixture.phenotypes, output_folder)
            nothing
        catch error
            error
        end
        @test err isa ErrorException
        @test occursin("mixed genotype precision", sprint(showerror, err))
        @test occursin("double_precision", sprint(showerror, err))
        @test !ispath(output_folder)
        @test fixture.model.MCMCinfo === initial_mcmcinfo
        @test all(i -> fixture.model.M[i].genotypes === fixture.parents[i], 1:2)
        @test eltype(fixture.parents[1]) === Float32
        @test eltype(fixture.parents[2]) === Float64

        # Explicit true bypasses inference and retains guarded legacy promotion.
        @test inferred_precision_small_run(fixture.model, fixture.phenotypes, output_folder;
            double_precision=true) isa AbstractDict
        @test fixture.model.MCMCinfo.double_precision === true
        @test all(g -> eltype(g.genotypes) === Float64, fixture.model.M)
        @test fixture.model.M[1].genotypes !== fixture.parents[1]
        @test fixture.model.M[2].genotypes === fixture.parents[2]

        # Explicit false also bypasses inference; strict sharing still performs
        # its existing one-term validation before any setup or allocation.
        explicit_false = inferred_precision_multiterm((Float32, Float64))
        @test_throws r"shared_genotypes=true currently requires exactly one genotype term" inferred_precision_small_run(
            explicit_false.model, explicit_false.phenotypes, joinpath(mktempdir(), "explicit_false");
            double_precision=false, shared_genotypes=true)

        unsupported = inferred_precision_multiterm((Float32, Float32))
        unsupported.model.M[1].genotypes = Float16.(unsupported.parents[1])
        unsupported_parent = unsupported.model.M[1].genotypes
        unsupported_folder = joinpath(mktempdir(), "unsupported")
        @test_throws r"unsupported genotype precision.*Float16" inferred_precision_small_run(
            unsupported.model, unsupported.phenotypes, unsupported_folder)
        @test !ispath(unsupported_folder)
        @test unsupported.model.M[1].genotypes === unsupported_parent
    end

    @testset "Streaming input keeps its Float32 backend" begin
        mktempdir() do tmpdir
            path = joinpath(tmpdir, "genotypes.csv")
            data = DataFrame(ID=string.(1:12), m1=repeat([0, 1, 2], 4),
                m2=repeat([2, 0, 1], 4), m3=repeat([1, 2, 0], 4))
            CSV.write(path, data)
            prefix = JWAS.prepare_streaming_genotypes(path; quality_control=false)
            global inferred_stream = get_genotypes(prefix, 1.0; method="BayesC", storage=:stream)
            model = build_model("y = intercept + inferred_stream", 1.0)
            parent, backend = inferred_stream.genotypes, inferred_stream.stream_backend
            @test JWAS.infer_run_double_precision(model) === false
            ph = DataFrame(ID=string.(1:12), y=sin.(1:12))
            @test inferred_precision_small_run(model, ph, joinpath(tmpdir, "stream")) isa AbstractDict
            @test model.MCMCinfo.double_precision === false
            @test eltype(model.sol) === Float32
            @test all(a -> eltype(a) === Float32, inferred_stream.α)
            @test inferred_stream.genotypes === parent
            @test inferred_stream.stream_backend === backend
        end
    end
end
