using Test
using JWAS
using DataFrames
using Random

function write_simple_genotype_file(path::AbstractString; with_missing::Bool=false)
    open(path, "w") do io
        println(io, "ID,m1,m2,m3,m4")
        println(io, "a1,0,1,2,0")
        println(io, "a2,1,0,1,2")
        println(io, with_missing ? "a3,2,9,0,1" : "a3,2,1,0,1")
        println(io, "a4,0,2,1,0")
        println(io, "a5,1,1,2,2")
        println(io, "a6,2,0,0,1")
    end
end

@testset "Streaming genotype codec and BayesC equivalence" begin
    mktempdir() do tmpdir
        cd(tmpdir) do
            @testset "prepare + load + decode" begin
                write_simple_genotype_file("geno_missing.csv"; with_missing=true)
                prefix = JWAS.prepare_streaming_genotypes("geno_missing.csv";
                                                          separator=',',
                                                          header=true,
                                                          quality_control=true,
                                                          center=true,
                                                          missing_value=9.0)
                @test isfile(prefix * ".jgb2")
                @test isfile(prefix * ".meta")

                dense = get_genotypes("geno_missing.csv", 1.0;
                                      separator=',',
                                      header=true,
                                      method="BayesC",
                                      quality_control=true,
                                      center=true)
                stream = get_genotypes(prefix, 1.0;
                                       method="BayesC",
                                       storage=:stream)

                @test stream.storage_mode == :stream
                @test stream.nObs == dense.nObs
                @test stream.nMarkers == dense.nMarkers

                marker_buf = zeros(Float32, stream.nObs)
                for j in 1:stream.nMarkers
                    JWAS.decode_marker!(marker_buf, stream.stream_backend, j)
                    @test marker_buf ≈ vec(Float32.(dense.genotypes[:, j])) atol=1e-5
                end
            end

            @testset "runMCMC dense vs stream" begin
                write_simple_genotype_file("geno_nomissing.csv"; with_missing=false)
                prefix = JWAS.prepare_streaming_genotypes("geno_nomissing.csv";
                                                          separator=',',
                                                          header=true,
                                                          quality_control=false,
                                                          center=true,
                                                          missing_value=9.0)

                ids = ["a1", "a2", "a3", "a4", "a5", "a6"]
                phenotypes = DataFrame(ID=ids,
                                       y1=Float32[1.1, -0.3, 0.8, -0.9, 0.5, -0.1])

                global geno = get_genotypes("geno_nomissing.csv", 1.0;
                                            separator=',',
                                            header=true,
                                            method="BayesC",
                                            quality_control=false,
                                            center=true)
                model_dense = build_model("y1 = intercept + geno", 1.0)
                out_dense = runMCMC(model_dense, phenotypes;
                                    chain_length=40,
                                    burnin=10,
                                    output_samples_frequency=10,
                                    output_folder="dense_results",
                                    outputEBV=false,
                                    output_heritability=false,
                                    seed=2026,
                                    memory_guard=:off)

                global geno = get_genotypes(prefix, 1.0;
                                            method="BayesC",
                                            storage=:stream)
                model_stream = build_model("y1 = intercept + geno", 1.0)
                out_stream = runMCMC(model_stream, phenotypes;
                                     chain_length=40,
                                     burnin=10,
                                     output_samples_frequency=10,
                                     output_folder="stream_results",
                                     outputEBV=false,
                                     output_heritability=false,
                                     seed=2026,
                                     memory_guard=:off)

                dense_eff = Vector{Float64}(out_dense["marker effects geno"][!, :Estimate])
                stream_eff = Vector{Float64}(out_stream["marker effects geno"][!, :Estimate])
                @test length(dense_eff) == length(stream_eff)
                @test stream_eff ≈ dense_eff atol=1e-4

                dense_resid = Float64(out_dense["residual variance"][1, :Estimate])
                stream_resid = Float64(out_stream["residual variance"][1, :Estimate])
                @test stream_resid ≈ dense_resid atol=1e-4
            end
        end
    end
end
