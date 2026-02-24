using Test
using JWAS
using DataFrames

function write_stream_test_geno(path::AbstractString)
    open(path, "w") do io
        println(io, "ID,m1,m2,m3")
        println(io, "a1,0,1,2")
        println(io, "a2,1,0,1")
        println(io, "a3,2,1,0")
        println(io, "a4,0,2,1")
    end
end

@testset "Streaming mode constraints" begin
    mktempdir() do tmpdir
        cd(tmpdir) do
            write_stream_test_geno("geno.csv")
            prefix = JWAS.prepare_streaming_genotypes("geno.csv";
                                                      separator=',',
                                                      header=true,
                                                      quality_control=false,
                                                      center=true)

            phenotypes = DataFrame(ID=["a1", "a2", "a3", "a4"],
                                   y1=Float32[1.0, -0.2, 0.5, -0.8])

            @testset "fast_blocks is rejected" begin
                global geno = get_genotypes(prefix, 1.0; method="BayesC", storage=:stream)
                model = build_model("y1 = intercept + geno", 1.0)
                @test_throws ErrorException runMCMC(model, phenotypes;
                                                    chain_length=20,
                                                    burnin=5,
                                                    output_samples_frequency=5,
                                                    output_folder="stream_fast_blocks",
                                                    outputEBV=false,
                                                    output_heritability=false,
                                                    seed=11,
                                                    fast_blocks=true,
                                                    memory_guard=:off)
            end

            @testset "non-unit residual weights are rejected" begin
                phenotypes_weighted = DataFrame(ID=["a1", "a2", "a3", "a4"],
                                                y1=Float32[1.0, -0.2, 0.5, -0.8],
                                                weights=Float32[1.0, 2.0, 1.0, 1.0])
                global geno = get_genotypes(prefix, 1.0; method="BayesC", storage=:stream)
                model = build_model("y1 = intercept + geno", 1.0)
                @test_throws ErrorException runMCMC(model, phenotypes_weighted;
                                                    chain_length=20,
                                                    burnin=5,
                                                    output_samples_frequency=5,
                                                    output_folder="stream_weighted",
                                                    outputEBV=false,
                                                    output_heritability=false,
                                                    heterogeneous_residuals=true,
                                                    seed=12,
                                                    memory_guard=:off)
            end

            @testset "multi-trait is rejected" begin
                pheno_mt = DataFrame(ID=["a1", "a2", "a3", "a4"],
                                     y1=Float32[1.0, -0.2, 0.5, -0.8],
                                     y2=Float32[0.1, 0.2, -0.1, 0.3])
                global geno = get_genotypes(prefix, [1.0 0.0; 0.0 1.0]; method="BayesC", storage=:stream)
                model_mt = build_model("y1 = intercept + geno\ny2 = intercept + geno", [1.0 0.0; 0.0 1.0])
                @test_throws ErrorException runMCMC(model_mt, pheno_mt;
                                                    chain_length=20,
                                                    burnin=5,
                                                    output_samples_frequency=5,
                                                    output_folder="stream_mt",
                                                    outputEBV=false,
                                                    output_heritability=false,
                                                    seed=13,
                                                    memory_guard=:off)
            end

            @testset "double_precision is rejected" begin
                global geno = get_genotypes(prefix, 1.0; method="BayesC", storage=:stream)
                model = build_model("y1 = intercept + geno", 1.0)
                @test_throws ErrorException runMCMC(model, phenotypes;
                                                    chain_length=20,
                                                    burnin=5,
                                                    output_samples_frequency=5,
                                                    output_folder="stream_double",
                                                    outputEBV=false,
                                                    output_heritability=false,
                                                    seed=14,
                                                    double_precision=true,
                                                    memory_guard=:off)
            end
        end
    end
end
