using Test, JWAS, DataFrames

function bayesc_revision_fixture(; storage_type=Float32)
    X = storage_type.(repeat([0 0 2 1 0 2; 1 1 1 0 2 1; 2 2 0 2 1 0], 4, 1))
    global bayesc_revision_genotypes = get_genotypes(X, 1.0; method="BayesC", Pi=0.5,
        quality_control=false, double_precision=(storage_type === Float64))
    model = build_model("y = intercept + bayesc_revision_genotypes", 1.0)
    phenotypes = DataFrame(ID=string.(bayesc_revision_genotypes.obsID), y=sin.(1:size(X, 1)))
    (; X, genotypes=bayesc_revision_genotypes, model, phenotypes)
end

function bayesc_revision_run(fixture; kwargs...)
    runMCMC(fixture.model, fixture.phenotypes; chain_length=12, burnin=2,
        output_samples_frequency=2, seed=123, printout_model_info=false,
        output_heritability=false, output_folder=joinpath(mktempdir(), "mcmc"), kwargs...)
end

@testset "BayesC no-copy regressions" begin
    @testset "Requested precision promotion uses the copying fallback" begin
        fixture = bayesc_revision_fixture()
        output = bayesc_revision_run(fixture; double_precision=true)
        @test output isa AbstractDict
        @test eltype(fixture.genotypes.genotypes) === Float64
        @test fixture.genotypes.genotypes !== fixture.X

        explicit = bayesc_revision_fixture()
        @test_throws r"shared_genotypes=true requires an initial Matrix\{Float64\}" bayesc_revision_run(
            explicit; double_precision=true, shared_genotypes=true)
        @test explicit.genotypes.genotypes === explicit.X
        @test eltype(explicit.X) === Float32
    end

    @testset "Precision conversion checks source and destination before allocating" begin
        fixture = bayesc_revision_fixture()
        peak_bytes = sizeof(fixture.X) + length(fixture.X) * sizeof(Float64)
        ratio = 0.5 * peak_bytes / Sys.total_memory()
        err = try
            bayesc_revision_run(fixture; double_precision=true, shared_genotypes=false,
                memory_guard_ratio=ratio)
            nothing
        catch error
            error
        end
        @test err isa ErrorException
        @test occursin("genotype precision conversion", sprint(showerror, err))
        @test occursin("retained_dense_parents=", sprint(showerror, err))
        @test occursin("destination=", sprint(showerror, err))
        @test fixture.genotypes.genotypes === fixture.X
        @test eltype(fixture.genotypes.genotypes) === Float32

        off = bayesc_revision_fixture()
        @test bayesc_revision_run(off; double_precision=true, shared_genotypes=false,
            memory_guard=:off, memory_guard_ratio=ratio) isa AbstractDict
        @test eltype(off.genotypes.genotypes) === Float64

        warned = bayesc_revision_fixture()
        @test_logs match_mode=:any (:warn, r"genotype precision conversion") begin
            bayesc_revision_run(warned; double_precision=true, shared_genotypes=false,
                memory_guard=:warn, memory_guard_ratio=ratio)
        end
        @test eltype(warned.genotypes.genotypes) === Float64
    end

    @testset "Conversion estimate counts unique parents and each destination" begin
        X = reshape(Float32.(1:72), 12, 6)
        a = get_genotypes(X; method="BayesC", quality_control=false, center=false)
        b = get_genotypes(X; method="BayesC", quality_control=false, center=false)
        b.genotypes = view(X, 1:8, :)
        b.output_genotypes = view(X, 9:12, :)
        b.MArray = [view(X, :, 1:3), view(X, :, 4:6)]
        Y = Float64.(X)
        c = get_genotypes(Y; method="BayesC", quality_control=false, center=false,
            double_precision=true)
        estimate = JWAS.estimate_genotype_precision_conversion([a, b, c])
        @test estimate.retained_bytes == sizeof(X) + sizeof(Y)
        @test estimate.destination_bytes == (length(a.genotypes) + length(b.genotypes)) * sizeof(Float64)
        @test estimate.peak_bytes == estimate.retained_bytes + estimate.destination_bytes
        @test JWAS.estimate_genotype_precision_conversion([c]).destination_bytes == 0
    end

    @testset "Dense memory guard uses actual stored precision" begin
        fixture = bayesc_revision_fixture(; storage_type=Float64)
        values = length(fixture.X) + fixture.genotypes.nMarkers
        ratio = values * 6 / Sys.total_memory() # Between Float32 and Float64 payloads.
        err = try
            bayesc_revision_run(fixture; double_precision=false, shared_genotypes=false,
                memory_guard_ratio=ratio)
            nothing
        catch error
            error
        end
        @test err isa ErrorException
        @test occursin("precision=Float64", sprint(showerror, err))
        @test eltype(fixture.genotypes.genotypes) === Float64
    end

    @testset "All-pass QC retains the adopted genotype parent" begin
        for T in (Float32, Float64), center in (false, true)
            X = T.(repeat([0 2 1; 1 1 0; 2 0 2], 4, 1))
            original = copy(X)
            annotations = reshape(Float64.(1:3), 3, 1)
            g = get_genotypes(X; method="BayesC", Pi=0.5,
                double_precision=(T === Float64), center, annotations)
            @test g.genotypes === X
            @test g.nMarkers == 3
            @test g.markerID == ["1", "2", "3"]
            @test g.alleleFreq == fill(T(0.5), 1, 3)
            @test g.genotypes == (center ? original .- one(T) : original)
            @test g.annotations.design_matrix == hcat(ones(3), annotations)
        end
    end

    @testset "Actual QC filtering keeps the existing data and metadata selection" begin
        X = Float32.(repeat([0 0 2 2; 1 0 1 2; 2 0 0 2], 4, 1))
        expected = copy(X[:, [1, 3]]) .- 1f0
        annotations = reshape(Float64[0, 0, 1, 1], 4, 1)
        g = get_genotypes(X; method="BayesC", Pi=0.5, annotations)
        @test g.genotypes !== X
        @test g.genotypes == expected
        @test g.nMarkers == 2
        @test g.markerID == ["1", "3"]
        @test g.alleleFreq == fill(0.5f0, 1, 2)
        @test g.annotations.design_matrix == hcat(ones(2), annotations[[1, 3], :])
    end

    @testset "Matching precision retains default adoption and row sharing" begin
        fixture = bayesc_revision_fixture()
        @test bayesc_revision_run(fixture) isa AbstractDict
        @test fixture.genotypes.genotypes === fixture.X
        @test fixture.genotypes.output_genotypes === fixture.X
        @test JWAS.genotype_storage_report(fixture.genotypes).backing_matrix_count == 1
    end
end
