using Test
using JWAS
using LinearAlgebra

function write_lowmem_genotype_file(path::AbstractString)
    open(path, "w") do io
        println(io, "ID,m1,m2,m3,m4,m5")
        println(io, "a1,0,1,2,0,1")
        println(io, "a2,1,0,1,2,0")
        println(io, "a3,2,9,0,1,2")
        println(io, "a4,0,2,1,0,1")
        println(io, "a5,1,1,2,2,0")
        println(io, "a6,2,0,0,1,2")
    end
end

@testset "Streaming low-memory prepare behavior" begin
    mktempdir() do tmpdir
        cd(tmpdir) do
            write_lowmem_genotype_file("geno.csv")

            @testset "QC and decode parity vs dense" begin
                prefix = JWAS.prepare_streaming_genotypes("geno.csv";
                                                          separator=',',
                                                          header=true,
                                                          quality_control=true,
                                                          center=true,
                                                          missing_value=9.0)
                dense = get_genotypes("geno.csv", 1.0;
                                      separator=',',
                                      header=true,
                                      method="BayesC",
                                      quality_control=true,
                                      center=true)
                stream = get_genotypes(prefix, 1.0;
                                       method="BayesC",
                                       storage=:stream)

                @test stream.markerID == dense.markerID
                @test stream.nMarkers == dense.nMarkers
                @test stream.nObs == dense.nObs

                marker_buf = zeros(Float32, stream.nObs)
                for j in 1:stream.nMarkers
                    JWAS.decode_marker!(marker_buf, stream.stream_backend, j)
                    @test marker_buf ≈ vec(Float32.(dense.genotypes[:, j])) atol=1e-5
                    @test dot(marker_buf, marker_buf) ≈ stream.stream_backend.xpRinvx[j] atol=1e-5
                end
            end

            @testset "xpRinvx parity when center=false" begin
                prefix = JWAS.prepare_streaming_genotypes("geno.csv";
                                                          output_prefix="geno_uncentered_stream",
                                                          separator=',',
                                                          header=true,
                                                          quality_control=false,
                                                          center=false,
                                                          missing_value=9.0)
                stream = get_genotypes(prefix, 1.0;
                                       method="BayesC",
                                       storage=:stream)
                marker_buf = zeros(Float32, stream.nObs)
                for j in 1:stream.nMarkers
                    JWAS.decode_marker!(marker_buf, stream.stream_backend, j)
                    @test dot(marker_buf, marker_buf) ≈ stream.stream_backend.xpRinvx[j] atol=1e-5
                end
            end

            @testset "cleanup_temp toggles stage artifact removal" begin
                tmproot = joinpath(tmpdir, "stream_stage")
                mkpath(tmproot)

                before = Set(readdir(tmproot))
                _ = JWAS.prepare_streaming_genotypes("geno.csv";
                                                     output_prefix="cleanup_true_stream",
                                                     separator=',',
                                                     header=true,
                                                     quality_control=false,
                                                     center=true,
                                                     tmpdir=tmproot,
                                                     cleanup_temp=true)
                after_cleanup = Set(readdir(tmproot))
                @test after_cleanup == before

                _ = JWAS.prepare_streaming_genotypes("geno.csv";
                                                     output_prefix="cleanup_false_stream",
                                                     separator=',',
                                                     header=true,
                                                     quality_control=false,
                                                     center=true,
                                                     tmpdir=tmproot,
                                                     cleanup_temp=false)
                after_keep = Set(readdir(tmproot))
                @test !isempty(setdiff(after_keep, before))
            end

            @testset "disk guard failure path" begin
                err = nothing
                try
                    JWAS.prepare_streaming_genotypes("geno.csv";
                                                     output_prefix="disk_guard_fail_stream",
                                                     separator=',',
                                                     header=true,
                                                     quality_control=false,
                                                     center=true,
                                                     disk_guard_ratio=1e-12)
                catch e
                    err = e
                end
                @test err !== nothing
                @test occursin("Insufficient disk", sprint(showerror, err))
            end
        end
    end
end
