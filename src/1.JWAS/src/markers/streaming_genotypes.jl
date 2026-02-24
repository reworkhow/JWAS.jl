"""
    Packed2BitBackend

Streaming backend for marker-major 2-bit packed genotypes used by
`storage=:stream` in `get_genotypes`.
"""
mutable struct Packed2BitBackend
    prefix::String
    data_path::String
    io::IOStream
    nObs::Int
    nMarkers::Int
    stride_bytes::Int
    centered::Bool
    marker_means::Vector{Float32}
    xpRinvx::Vector{Float32}
    allele_freq::Vector{Float32}
    sum2pq::Float64
    obsID::Vector{String}
    markerID::Vector{String}
    packed_buffer::Vector{UInt8}
end

function close_streaming_backend!(backend::Packed2BitBackend)
    if isopen(backend.io)
        close(backend.io)
    end
    return nothing
end

function _write_f32_vector(path::AbstractString, values::AbstractVector{Float32})
    open(path, "w") do io
        write(io, values)
    end
    return nothing
end

function _read_f32_vector(path::AbstractString, n::Int)
    vals = Vector{Float32}(undef, n)
    open(path, "r") do io
        read!(io, vals)
    end
    return vals
end

function _write_string_lines(path::AbstractString, values::Vector{String})
    open(path, "w") do io
        for value in values
            println(io, value)
        end
    end
    return nothing
end

function _read_string_lines(path::AbstractString)
    return filter(!isempty, readlines(path))
end

function _write_streaming_manifest(path::AbstractString, entries::Vector{Pair{String,String}})
    open(path, "w") do io
        for (k, v) in entries
            println(io, k, "\t", v)
        end
    end
    return nothing
end

function _read_streaming_manifest(path::AbstractString)
    entries = Dict{String,String}()
    for line in eachline(path)
        parts = split(line, '\t', limit=2)
        if length(parts) == 2
            entries[parts[1]] = parts[2]
        end
    end
    return entries
end

function _resolve_streaming_prefix(path::AbstractString)
    if endswith(path, ".meta")
        return path[1:end-5]
    elseif endswith(path, ".jgb2")
        return path[1:end-5]
    else
        return path
    end
end

"""
    prepare_streaming_genotypes(file::AbstractString;
                                output_prefix=nothing,
                                separator=',',
                                header=true,
                                missing_value=9.0,
                                quality_control=true,
                                MAF=0.01,
                                center=true)

Convert a dense text genotype file to a marker-major 2-bit packed backend.
"""
function prepare_streaming_genotypes(file::AbstractString;
                                     output_prefix=nothing,
                                     separator=',',
                                     header=true,
                                     missing_value=9.0,
                                     quality_control=true,
                                     MAF=0.01,
                                     center=true)
    output_prefix = output_prefix === nothing ? splitext(file)[1] * "_stream" : string(output_prefix)
    prefix_abs = abspath(output_prefix)
    data_path = prefix_abs * ".jgb2"
    meta_path = prefix_abs * ".meta"
    obs_path = prefix_abs * ".obsid.txt"
    marker_path = prefix_abs * ".markerid.txt"
    mean_path = prefix_abs * ".mean.f32"
    xp_path = prefix_abs * ".xpRinvx.f32"
    afreq_path = prefix_abs * ".afreq.f32"

    myfile = open(file)
    row1 = split(readline(myfile), [separator, '\n'], keepempty=false)
    close(myfile)

    marker_id_all = header ? string.(row1[2:end]) : string.(1:length(row1[2:end]))

    ncol = length(row1)
    etv = Array{DataType}(undef, ncol)
    fill!(etv, Float32)
    etv[1] = String
    data = CSV.read(file, DataFrame, types=etv, delim=separator, header=false, skipto=(header ? 2 : 1))

    obsID = map(string, data[!, 1])
    X = map(Float32, Matrix(data[!, 2:end]))
    nObs, nMarkersAll = size(X)

    if nObs == 0 || nMarkersAll == 0
        error("Genotype data is empty.")
    end

    missing_value32 = Float32(missing_value)
    maf32 = Float32(MAF)
    keep = trues(nMarkersAll)
    marker_means = Vector{Float32}(undef, nMarkersAll)
    allele_freq = Vector{Float32}(undef, nMarkersAll)
    xpRinvx_all = Vector{Float32}(undef, nMarkersAll)
    missing_masks = Vector{BitVector}(undef, nMarkersAll)

    @showprogress "preparing packed genotypes ..." for j in 1:nMarkersAll
        col = view(X, :, j)
        miss = BitVector(undef, nObs)
        n_nonmissing = 0
        sum_nonmissing = 0.0f0

        @inbounds for i in 1:nObs
            v = col[i]
            is_missing = (v == missing_value32)
            miss[i] = is_missing
            if !is_missing
                if !(v == 0.0f0 || v == 1.0f0 || v == 2.0f0)
                    error("Only 0/1/2 genotypes (and missing_value=$missing_value) are supported in storage=:stream MVP.")
                end
                n_nonmissing += 1
                sum_nonmissing += v
            end
        end

        if n_nonmissing == 0
            error("Marker $(marker_id_all[j]) has only missing values.")
        end

        μ = sum_nonmissing / n_nonmissing
        marker_means[j] = μ
        allele_freq[j] = μ / 2.0f0
        missing_masks[j] = miss

        ss_centered = 0.0f0
        ss_raw = 0.0f0
        @inbounds for i in 1:nObs
            v = miss[i] ? μ : col[i]
            d = v - μ
            ss_centered += d * d
            ss_raw += v * v
        end
        xpRinvx_all[j] = center ? ss_centered : ss_raw

        if quality_control
            maf_ok = (maf32 < allele_freq[j] < (1.0f0 - maf32))
            var_ok = (ss_centered != 0.0f0)
            keep[j] = maf_ok && var_ok
        end
    end

    selected = quality_control ? findall(keep) : collect(1:nMarkersAll)
    if isempty(selected)
        error("No markers remain after streaming genotype quality control.")
    end

    markerID = marker_id_all[selected]
    means_kept = marker_means[selected]
    afreq_kept = allele_freq[selected]
    xp_kept = xpRinvx_all[selected]
    nMarkers = length(selected)
    stride_bytes = cld(nObs, 4)
    sum2pq = Float64(sum(2.0f0 .* afreq_kept .* (1.0f0 .- afreq_kept)))

    if quality_control
        printstyled("$(nMarkersAll - nMarkers) loci which are fixed or have minor allele frequency < $MAF are removed.\n", bold=true)
    end

    packed_row = zeros(UInt8, stride_bytes)
    open(data_path, "w") do io
        @showprogress "writing packed genotypes ..." for idx in 1:nMarkers
            j = selected[idx]
            miss = missing_masks[j]
            col = view(X, :, j)
            fill!(packed_row, 0x00)

            @inbounds for i in 1:nObs
                code = miss[i] ? UInt8(3) : UInt8(round(Int, col[i]))
                byte_idx = ((i - 1) >>> 2) + 1
                shift = ((i - 1) & 0x03) << 1
                packed_row[byte_idx] |= code << shift
            end
            write(io, packed_row)
        end
    end

    _write_string_lines(obs_path, obsID)
    _write_string_lines(marker_path, markerID)
    _write_f32_vector(mean_path, means_kept)
    _write_f32_vector(xp_path, xp_kept)
    _write_f32_vector(afreq_path, afreq_kept)

    _write_streaming_manifest(meta_path, [
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
        "centered" => string(center ? 1 : 0),
        "sum2pq" => string(sum2pq)
    ])

    printstyled("Streaming genotype files are created with prefix $(prefix_abs).\n", bold=false, color=:green)
    return prefix_abs
end

"""
    load_streaming_backend(path::AbstractString)

Load a packed genotype backend from a prefix (or `.jgb2` / `.meta` path).
"""
function load_streaming_backend(path::AbstractString)
    prefix = _resolve_streaming_prefix(abspath(path))
    meta_path = prefix * ".meta"
    if !isfile(meta_path)
        error("Streaming manifest is not found: $meta_path")
    end

    meta = _read_streaming_manifest(meta_path)
    nObs = parse(Int, meta["nObs"])
    nMarkers = parse(Int, meta["nMarkers"])
    stride_bytes = parse(Int, meta["stride_bytes"])
    centered = parse(Int, meta["centered"]) == 1
    sum2pq = parse(Float64, meta["sum2pq"])

    data_path = meta["data_path"]
    obs_path = meta["obs_path"]
    marker_path = meta["marker_path"]
    mean_path = meta["mean_path"]
    xp_path = meta["xp_path"]
    afreq_path = meta["afreq_path"]

    expected_bytes = Int64(nMarkers) * Int64(stride_bytes)
    if filesize(data_path) != expected_bytes
        error("Packed genotype file size does not match metadata for $data_path")
    end

    marker_means = _read_f32_vector(mean_path, nMarkers)
    xpRinvx = _read_f32_vector(xp_path, nMarkers)
    allele_freq = _read_f32_vector(afreq_path, nMarkers)
    obsID = _read_string_lines(obs_path)
    markerID = _read_string_lines(marker_path)

    if length(obsID) != nObs
        error("Number of IDs in $obs_path does not match nObs in manifest.")
    end
    if length(markerID) != nMarkers
        error("Number of markers in $marker_path does not match nMarkers in manifest.")
    end

    io = open(data_path, "r")
    backend = Packed2BitBackend(
        prefix,
        data_path,
        io,
        nObs,
        nMarkers,
        stride_bytes,
        centered,
        marker_means,
        xpRinvx,
        allele_freq,
        sum2pq,
        obsID,
        markerID,
        zeros(UInt8, stride_bytes)
    )

    finalizer(backend) do b
        close_streaming_backend!(b)
    end

    return backend
end

"""
    decode_marker!(dest, backend, marker_index)

Decode one marker into `dest` with centering and missing-value mean-imputation.
"""
function decode_marker!(dest::AbstractVector, backend::Packed2BitBackend, marker_index::Integer)
    if marker_index < 1 || marker_index > backend.nMarkers
        error("marker_index out of bounds.")
    end
    if length(dest) != backend.nObs
        error("Destination vector length must equal nObs.")
    end

    offset = (marker_index - 1) * backend.stride_bytes
    seek(backend.io, offset)
    read!(backend.io, backend.packed_buffer)

    μ = backend.marker_means[marker_index]
    centered = backend.centered

    @inbounds for i in 1:backend.nObs
        byte_idx = ((i - 1) >>> 2) + 1
        shift = ((i - 1) & 0x03) << 1
        code = (backend.packed_buffer[byte_idx] >> shift) & 0x03
        v = code == 0x03 ? μ : Float32(code)
        dest[i] = centered ? (v - μ) : v
    end

    return dest
end

"""
    streaming_mul_alpha!(out, backend, α)

Compute `out = X*α` from a streaming backend without materializing dense `X`.
"""
function streaming_mul_alpha!(out::AbstractVector, backend::Packed2BitBackend, α::AbstractVector)
    if length(out) != backend.nObs
        error("Output vector length must equal nObs.")
    end
    if length(α) != backend.nMarkers
        error("length(α) must equal nMarkers.")
    end

    fill!(out, zero(eltype(out)))
    marker_buffer = Vector{eltype(out)}(undef, backend.nObs)
    for j in 1:backend.nMarkers
        αj = α[j]
        if αj != zero(αj)
            decode_marker!(marker_buffer, backend, j)
            BLAS.axpy!(αj, marker_buffer, out)
        end
    end
    return out
end
