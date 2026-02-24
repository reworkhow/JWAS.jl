using JWAS
using DelimitedFiles
using LinearAlgebra
using Random
using DataFrames
using CSV
using Dates

function get_arg_value(args::Vector{String}, key::String, default::String)
    prefix = "--" * key * "="
    for arg in args
        if startswith(arg, prefix)
            return arg[length(prefix)+1:end]
        end
    end
    return default
end

function get_arg_bool(args::Vector{String}, key::String, default::Bool)
    raw = lowercase(get_arg_value(args, key, default ? "true" : "false"))
    if raw in ("1","true","yes","y","on")
        return true
    elseif raw in ("0","false","no","n","off")
        return false
    end
    error("Invalid boolean value for --$key: $raw")
end

function get_arg_int(args::Vector{String}, key::String, default::Int)
    return parse(Int, get_arg_value(args, key, string(default)))
end

function get_arg_float(args::Vector{String}, key::String, default::Float64)
    return parse(Float64, get_arg_value(args, key, string(default)))
end

function get_arg_int_list(args::Vector{String}, key::String, default::Vector{Int})
    raw = get_arg_value(args, key, join(default, ","))
    if isempty(strip(raw))
        return Int[]
    end
    return [parse(Int, strip(x)) for x in split(raw, ",")]
end

function rss_kb()
    pid = getpid()
    out = read(`ps -o rss= -p $pid`, String)
    raw = strip(out)
    if isempty(raw)
        return 0
    end
    return parse(Int, split(raw)[1])
end

function start_rss_sampler(interval_sec::Float64)
    max_rss = Ref(rss_kb())
    keep_running = Ref(true)
    task = @async begin
        while keep_running[]
            sleep(interval_sec)
            current = rss_kb()
            if current > max_rss[]
                max_rss[] = current
            end
        end
    end
    return keep_running, max_rss, task
end

function stop_rss_sampler!(keep_running::Base.RefValue{Bool}, task::Task)
    keep_running[] = false
    wait(task)
    return nothing
end

function write_meta(path::String, kvs::Vector{Pair{String,String}})
    open(path, "w") do io
        for (k, v) in kvs
            println(io, k, "\t", v)
        end
    end
    return nothing
end

function ensure_sparse_stream_backend(prefix::String, nObs::Int, nMarkers::Int; centered::Bool=true)
    prefix = abspath(prefix)
    data_path = prefix * ".jgb2"
    meta_path = prefix * ".meta"
    obs_path = prefix * ".obsid.txt"
    marker_path = prefix * ".markerid.txt"
    mean_path = prefix * ".mean.f32"
    xp_path = prefix * ".xpRinvx.f32"
    afreq_path = prefix * ".afreq.f32"
    stride_bytes = cld(nObs, 4)
    total_bytes = Int64(nMarkers) * Int64(stride_bytes)

    println("Creating sparse streaming backend at prefix: ", prefix)
    println("nObs=", nObs, ", nMarkers=", nMarkers, ", stride_bytes=", stride_bytes)
    println("packed logical bytes=", total_bytes)

    # Create sparse packed payload.
    open(data_path, "w") do io
        seek(io, total_bytes - 1)
        write(io, UInt8(0))
    end

    open(obs_path, "w") do io
        for i in 1:nObs
            println(io, "id", i)
        end
    end

    open(marker_path, "w") do io
        for j in 1:nMarkers
            println(io, "m", j)
        end
    end

    means = fill(Float32(1.0), nMarkers)
    xp = fill(Float32(nObs), nMarkers)
    afreq = fill(Float32(0.5), nMarkers)
    open(mean_path, "w") do io
        write(io, means)
    end
    open(xp_path, "w") do io
        write(io, xp)
    end
    open(afreq_path, "w") do io
        write(io, afreq)
    end

    write_meta(meta_path, [
        "version" => "1",
        "data_path" => data_path,
        "obs_path" => obs_path,
        "marker_path" => marker_path,
        "mean_path" => mean_path,
        "xp_path" => xp_path,
        "afreq_path" => afreq_path,
        "nObs" => string(nObs),
        "nMarkers" => string(nMarkers),
        "stride_bytes" => string(stride_bytes),
        "centered" => centered ? "1" : "0",
        "sum2pq" => "0.0",
    ])

    return prefix
end

function benchmark_streaming_sweeps(prefix::String;
                                    nSweeps::Int,
                                    sampleMarkersPerSweep::Int,
                                    samplerIntervalSec::Float64,
                                    seed::Int)
    Random.seed!(seed)
    geno = get_genotypes(prefix, 1.0; method="BayesC", storage=:stream)
    backend = geno.stream_backend

    nObs = backend.nObs
    nMarkers = backend.nMarkers
    markers_this_sweep = sampleMarkersPerSweep <= 0 ? nMarkers : min(sampleMarkersPerSweep, nMarkers)
    total_markers_processed = Int64(markers_this_sweep) * Int64(nSweeps)

    ycorr = randn(Float32, nObs)
    xbuf = Vector{Float32}(undef, nObs)

    keep_running, max_rss_ref, sampler_task = start_rss_sampler(samplerIntervalSec)
    timed = @timed begin
        sink = 0.0f0
        for _ in 1:nSweeps
            for j in 1:markers_this_sweep
                JWAS.decode_marker!(xbuf, backend, j)
                sink += dot(xbuf, ycorr)
                BLAS.axpy!(0.01f0, xbuf, ycorr)
            end
        end
        sink
    end
    stop_rss_sampler!(keep_running, sampler_task)

    elapsed = timed.time
    markers_per_sec = total_markers_processed / elapsed
    projected_full_sweep_sec = nMarkers / markers_per_sec
    projected_full_sweep_min = projected_full_sweep_sec / 60
    projected_full_sweep_hour = projected_full_sweep_sec / 3600

    est_dense = JWAS.estimate_marker_memory(nObs, nMarkers;
                                       element_bytes=4,
                                       has_nonunit_weights=false,
                                       block_starts=false,
                                       storage_mode=:dense)
    est_stream = JWAS.estimate_marker_memory(nObs, nMarkers;
                                        element_bytes=4,
                                        has_nonunit_weights=false,
                                        block_starts=false,
                                        storage_mode=:stream)

    return (
        timestamp = string(now()),
        prefix = prefix,
        nObs = nObs,
        nMarkers = nMarkers,
        nSweeps = nSweeps,
        markers_this_sweep = markers_this_sweep,
        total_markers_processed = total_markers_processed,
        elapsed_sec = elapsed,
        allocated_bytes = timed.bytes,
        peak_rss_kb = max_rss_ref[],
        markers_per_sec = markers_per_sec,
        projected_full_sweep_sec = projected_full_sweep_sec,
        projected_full_sweep_min = projected_full_sweep_min,
        projected_full_sweep_hour = projected_full_sweep_hour,
        system_ram_bytes = Sys.total_memory(),
        dense_est_bytes = est_dense.bytes_total,
        stream_est_bytes = est_stream.bytes_total,
        dense_est_human = JWAS.format_bytes_human(est_dense.bytes_total),
        stream_est_human = JWAS.format_bytes_human(est_stream.bytes_total),
    )
end

function main(args::Vector{String})
    nObs = get_arg_int(args, "n-obs", 500_000)
    nMarkers = get_arg_int(args, "n-markers", 2_000_000)
    prefix = get_arg_value(args, "prefix", "/tmp/jwas_stream_500k_2m")
    create_backend = get_arg_bool(args, "create-backend", false)
    force_backend = get_arg_bool(args, "force-create", false)
    sweeps = get_arg_int_list(args, "sweeps", [1, 5, 10])
    sample_markers = get_arg_int(args, "sample-markers-per-sweep", 0)
    sampler_interval = get_arg_float(args, "rss-sample-interval-sec", 0.2)
    seed = get_arg_int(args, "seed", 2026)
    csv_out = get_arg_value(args, "csv-out", "benchmarks/streaming_large_benchmark_results.csv")

    prefix = abspath(prefix)
    csv_out = abspath(csv_out)
    mkpath(dirname(csv_out))

    if create_backend
        if force_backend || !isfile(prefix * ".meta")
            ensure_sparse_stream_backend(prefix, nObs, nMarkers; centered=true)
        else
            println("Backend metadata already exists; reuse: ", prefix * ".meta")
        end
    elseif !isfile(prefix * ".meta")
        error("Backend metadata not found: " * (prefix * ".meta") * ". Use --create-backend=true.")
    end

    println("Benchmark configuration:")
    println("  prefix=", prefix)
    println("  sweeps=", sweeps)
    println("  sample-markers-per-sweep=", sample_markers <= 0 ? "FULL" : sample_markers)
    println("  csv-out=", csv_out)

    results = DataFrame()
    for s in sweeps
        println("\nRunning sweep benchmark: ", s, " sweep(s)")
        r = benchmark_streaming_sweeps(prefix;
                                       nSweeps=s,
                                       sampleMarkersPerSweep=sample_markers,
                                       samplerIntervalSec=sampler_interval,
                                       seed=seed)
        push!(results, (; r...))

        println("  elapsed_sec=", round(r.elapsed_sec, digits=3))
        println("  markers_per_sec=", round(r.markers_per_sec, digits=2))
        println("  projected_full_sweep_hour=", round(r.projected_full_sweep_hour, digits=3))
        println("  peak_rss_kb=", r.peak_rss_kb)
        println("  dense_est_human=", r.dense_est_human)
        println("  stream_est_human=", r.stream_est_human)
    end

    CSV.write(csv_out, results)
    println("\nBenchmark results written to: ", csv_out)
    return nothing
end

main(ARGS)
