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

@inline function _parse_genotype_token(token::AbstractString, missing_value32::Float32, missing_value::Real)
    token_stripped = strip(token)
    if isempty(token_stripped)
        error("Empty genotype value is not supported in storage=:stream MVP.")
    end

    parsed = tryparse(Float32, token_stripped)
    if parsed === nothing
        error("Non-numeric genotype value '$token_stripped' is not supported in storage=:stream MVP.")
    end
    value = parsed::Float32

    if value == missing_value32
        return (value, true)
    end
    if value == 0.0f0 || value == 1.0f0 || value == 2.0f0
        return (value, false)
    end
    error("Only 0/1/2 genotypes (and missing_value=$missing_value) are supported in storage=:stream MVP.")
end

function _read_first_line(path::AbstractString)
    open(path, "r") do io
        if eof(io)
            error("Genotype data is empty.")
        end
        return readline(io)
    end
end

function _header_marker_ids(file::AbstractString, separator, header::Bool)
    first_line = _read_first_line(file)
    fields = collect(eachsplit(chomp(first_line), separator; keepempty=true))
    if length(fields) < 2
        error("Genotype file must contain an ID column and at least one marker column.")
    end
    nmarkers = length(fields) - 1
    marker_ids = header ? map(x -> String(x), fields[2:end]) : string.(1:nmarkers)
    return marker_ids, nmarkers
end

function _df_available_bytes(path::AbstractString)
    parent = isdir(path) ? path : dirname(path)
    if isempty(parent)
        parent = "."
    end
    out = read(`df -kP $parent`, String)
    lines = split(chomp(out), '\n')
    if length(lines) < 2
        error("Unable to read free-space information for path '$parent'.")
    end
    fields = split(strip(lines[end]))
    if length(fields) < 4
        error("Unable to parse free-space information for path '$parent'.")
    end
    available_kb = parse(Int128, fields[4])
    filesystem = fields[1]
    return available_kb * Int128(1024), filesystem
end

function _check_streaming_disk_guard!(output_dir::AbstractString,
                                      temp_dir::AbstractString,
                                      required_output_bytes::Int128,
                                      required_temp_bytes::Int128,
                                      disk_guard_ratio::Float64)
    if !(0.0 <= disk_guard_ratio <= 1.0)
        error("disk_guard_ratio must be between 0.0 and 1.0.")
    end

    # Skip guard when ratio is exactly 1.0 and requirements are tiny, but keep deterministic behavior.
    avail_out, fs_out = _df_available_bytes(output_dir)
    avail_tmp, fs_tmp = _df_available_bytes(temp_dir)

    if fs_out == fs_tmp
        required_total = required_output_bytes + required_temp_bytes
        threshold = floor(Int128, avail_out * disk_guard_ratio)
        if required_total > threshold
            error("Insufficient disk for streaming conversion (same filesystem).\n" *
                  "required (output + temp): $(required_total) bytes\n" *
                  "available: $(avail_out) bytes\n" *
                  "threshold via disk_guard_ratio=$(disk_guard_ratio): $(threshold) bytes\n" *
                  "Try a larger filesystem via tmpdir, reduce data size, or lower disk_guard_ratio.")
        end
        return nothing
    end

    threshold_out = floor(Int128, avail_out * disk_guard_ratio)
    threshold_tmp = floor(Int128, avail_tmp * disk_guard_ratio)
    if required_output_bytes > threshold_out || required_temp_bytes > threshold_tmp
        error("Insufficient disk for streaming conversion (separate filesystems).\n" *
              "required output: $(required_output_bytes) bytes; available output: $(avail_out) bytes; threshold: $(threshold_out) bytes\n" *
              "required temp: $(required_temp_bytes) bytes; available temp: $(avail_tmp) bytes; threshold: $(threshold_tmp) bytes\n" *
              "Try a larger tmpdir/output path, reduce data size, or lower disk_guard_ratio.")
    end
    return nothing
end

function _scan_streaming_stats!(file::AbstractString,
                                separator,
                                header::Bool,
                                nmarkers::Int,
                                marker_ids::Vector{String},
                                missing_value::Real,
                                quality_control::Bool,
                                maf::Real,
                                center::Bool,
                                obsid_stage_path::AbstractString)
    missing_value32 = Float32(missing_value)
    maf32 = Float32(maf)

    nonmissing_counts = zeros(Int64, nmarkers)
    missing_counts = zeros(Int64, nmarkers)
    sums = zeros(Float32, nmarkers)
    sumsq = zeros(Float32, nmarkers)
    keep = trues(nmarkers)
    marker_means = zeros(Float32, nmarkers)
    allele_freq = zeros(Float32, nmarkers)
    xp_all = zeros(Float32, nmarkers)

    nobs = 0
    open(obsid_stage_path, "w") do obs_io
        open(file, "r") do io
            if header
                if eof(io)
                    error("Genotype data is empty.")
                end
                readline(io) # skip header line
            end

            for line in eachline(io)
                nobs += 1
                token_iter = eachsplit(chomp(line), separator; keepempty=true)
                first_state = iterate(token_iter)
                first_state === nothing && error("Encountered an empty genotype row.")
                obs_token, state = first_state
                println(obs_io, obs_token)

                j = 0
                while true
                    next_state = iterate(token_iter, state)
                    next_state === nothing && break
                    token, state = next_state
                    j += 1
                    if j > nmarkers
                        error("Row $nobs has more marker columns than expected ($nmarkers).")
                    end
                    value, is_missing = _parse_genotype_token(token, missing_value32, missing_value)
                    if is_missing
                        missing_counts[j] += 1
                    else
                        nonmissing_counts[j] += 1
                        sums[j] += value
                        sumsq[j] += value * value
                    end
                end
                if j != nmarkers
                    error("Row $nobs has $j marker columns but expected $nmarkers.")
                end
            end
        end
    end

    if nobs == 0
        error("Genotype data is empty.")
    end

    @showprogress "preparing packed genotypes ..." for j in 1:nmarkers
        nn = nonmissing_counts[j]
        if nn == 0
            marker_name = marker_ids[j]
            error("Marker $(marker_name) has only missing values.")
        end

        μ = sums[j] / nn
        marker_means[j] = μ
        allele_freq[j] = μ / 2.0f0
        ss_centered = sumsq[j] - μ * sums[j]
        ss_raw = sumsq[j] + missing_counts[j] * (μ * μ)
        xp_all[j] = center ? ss_centered : ss_raw

        if quality_control
            maf_ok = (maf32 < allele_freq[j] < (1.0f0 - maf32))
            var_ok = (ss_centered != 0.0f0)
            keep[j] = maf_ok && var_ok
        end
    end

    selected = quality_control ? findall(keep) : collect(1:nmarkers)
    if isempty(selected)
        error("No markers remain after streaming genotype quality control.")
    end

    marker_kept = marker_ids[selected]
    means_kept = marker_means[selected]
    afreq_kept = allele_freq[selected]
    xp_kept = xp_all[selected]
    sum2pq = Float64(sum(2.0f0 .* afreq_kept .* (1.0f0 .- afreq_kept)))

    return (
        nObs = nobs,
        nMarkersAll = nmarkers,
        selected = selected,
        markerID = marker_kept,
        means = means_kept,
        afreq = afreq_kept,
        xp = xp_kept,
        sum2pq = sum2pq
    )
end

function _write_row_major_spool!(file::AbstractString,
                                 separator,
                                 header::Bool,
                                 nmarkers_all::Int,
                                 selected_map::Vector{Int32},
                                 nmarkers_kept::Int,
                                 nobs::Int,
                                 missing_value::Real,
                                 rowmajor_path::AbstractString)
    missing_value32 = Float32(missing_value)
    row_stride = cld(nmarkers_kept, 4)
    packed_row = zeros(UInt8, row_stride)
    rows_seen = 0

    open(rowmajor_path, "w") do out_io
        open(file, "r") do io
            if header
                if eof(io)
                    error("Genotype data is empty.")
                end
                readline(io) # skip header line
            end

            prog = Progress(nobs; desc="writing packed genotypes ...")
            for line in eachline(io)
                rows_seen += 1
                next!(prog)
                fill!(packed_row, 0x00)

                token_iter = eachsplit(chomp(line), separator; keepempty=true)
                first_state = iterate(token_iter)
                first_state === nothing && error("Encountered an empty genotype row.")
                _, state = first_state # skip ID token

                j = 0
                while true
                    next_state = iterate(token_iter, state)
                    next_state === nothing && break
                    token, state = next_state
                    j += 1
                    if j > nmarkers_all
                        error("Row $rows_seen has more marker columns than expected ($nmarkers_all).")
                    end

                    pos = selected_map[j]
                    if pos > 0
                        value, is_missing = _parse_genotype_token(token, missing_value32, missing_value)
                        code = is_missing ? UInt8(3) : UInt8(round(Int, value))
                        byte_idx = ((pos - 1) >>> 2) + 1
                        shift = ((pos - 1) & 0x03) << 1
                        packed_row[byte_idx] |= code << shift
                    end
                end
                if j != nmarkers_all
                    error("Row $rows_seen has $j marker columns but expected $nmarkers_all.")
                end
                write(out_io, packed_row)
            end
            finish!(prog)
        end
    end

    if rows_seen != nobs
        error("Number of rows changed between conversion passes ($rows_seen vs expected $nobs).")
    end
    return row_stride
end

function _transpose_row_major_to_marker_major!(rowmajor_path::AbstractString,
                                               markermajor_path::AbstractString,
                                               nobs::Int,
                                               nmarkers::Int)
    row_stride = cld(nmarkers, 4)
    marker_stride = cld(nobs, 4)
    total_bytes = Int64(nmarkers) * Int64(marker_stride)
    row_tile_target_bytes = 32 * 1024 * 1024
    marker_tile = 4096

    row_tile = max(4, Int(floor(row_tile_target_bytes / max(row_stride, 1))))
    row_tile = min(row_tile, 4096)
    row_tile -= row_tile % 4
    row_tile = max(row_tile, 4)
    max_row_tile_bytes = row_tile * row_stride
    max_row_tile_packed_bytes = cld(row_tile, 4)

    raw_rows = Vector{UInt8}(undef, max_row_tile_bytes)
    tile_packed = Matrix{UInt8}(undef, marker_tile, max_row_tile_packed_bytes)
    src_byte_idx = Vector{Int}(undef, marker_tile)
    src_shift = Vector{UInt8}(undef, marker_tile)

    open(markermajor_path, "w+") do out_io
        seek(out_io, total_bytes - 1)
        write(out_io, UInt8(0))
        flush(out_io)
        seekstart(out_io)
        open(rowmajor_path, "r") do in_io
            @showprogress "transposing packed genotypes ..." for row_start in 1:row_tile:nobs
                rows_this_tile = min(row_tile, nobs - row_start + 1)
                row_bytes_this_tile = rows_this_tile * row_stride
                read!(in_io, view(raw_rows, 1:row_bytes_this_tile))

                packed_row_bytes_this_tile = cld(rows_this_tile, 4)
                byte_start = ((row_start - 1) >>> 2) + 1

                for marker_start in 1:marker_tile:nmarkers
                    markers_this_tile = min(marker_tile, nmarkers - marker_start + 1)
                    fill!(view(tile_packed, 1:markers_this_tile, 1:packed_row_bytes_this_tile), 0x00)

                    @inbounds for mj in 1:markers_this_tile
                        g = marker_start + mj - 1
                        src_byte_idx[mj] = ((g - 1) >>> 2) + 1
                        src_shift[mj] = UInt8(((g - 1) & 0x03) << 1)
                    end

                    @inbounds for r in 1:rows_this_tile
                        row_offset = (r - 1) * row_stride
                        packed_byte = ((r - 1) >>> 2) + 1
                        packed_shift = UInt8(((r - 1) & 0x03) << 1)
                        for mj in 1:markers_this_tile
                            code = (raw_rows[row_offset + src_byte_idx[mj]] >> src_shift[mj]) & 0x03
                            tile_packed[mj, packed_byte] |= code << packed_shift
                        end
                    end

                    @inbounds for mj in 1:markers_this_tile
                        g = marker_start + mj - 1
                        dest_start = (g - 1) * marker_stride + byte_start
                        seek(out_io, dest_start - 1)
                        write(out_io, view(tile_packed, mj, 1:packed_row_bytes_this_tile))
                    end
                end
            end
        end
    end
    return marker_stride
end

function _safe_rm(path::AbstractString)
    if isfile(path)
        rm(path; force=true)
    end
    return nothing
end

function _count_genotype_rows(file::AbstractString, header::Bool)
    nobs = 0
    open(file, "r") do io
        if header
            if eof(io)
                error("Genotype data is empty.")
            end
            readline(io)
        end
        for _ in eachline(io)
            nobs += 1
        end
    end
    if nobs == 0
        error("Genotype data is empty.")
    end
    return nobs
end

function _normalize_conversion_mode(conversion_mode)
    mode = conversion_mode isa Symbol ? conversion_mode : Symbol(conversion_mode)
    if mode in (:auto, :dense, :lowmem)
        return mode
    end
    error("conversion_mode must be one of :auto, :dense, :lowmem.")
end

function _choose_conversion_mode(requested_mode::Symbol,
                                 nobs::Int,
                                 nmarkers::Int,
                                 auto_dense_max_bytes::Int128)
    if requested_mode != :auto
        return requested_mode
    end
    dense_bytes = Int128(nobs) * Int128(nmarkers) * Int128(sizeof(Float32))
    return dense_bytes <= auto_dense_max_bytes ? :dense : :lowmem
end

function _prepare_streaming_genotypes_dense!(file::AbstractString,
                                             prefix_abs::AbstractString;
                                             separator,
                                             header::Bool,
                                             missing_value::Real,
                                             quality_control::Bool,
                                             MAF::Real,
                                             center::Bool,
                                             disk_guard_ratio::Float64,
                                             marker_ids_all::Vector{String},
                                             nmarkers_all::Int)
    output_dir = dirname(prefix_abs)
    mkpath(output_dir)

    data_path = prefix_abs * ".jgb2"
    meta_path = prefix_abs * ".meta"
    obs_path = prefix_abs * ".obsid.txt"
    marker_path = prefix_abs * ".markerid.txt"
    mean_path = prefix_abs * ".mean.f32"
    xp_path = prefix_abs * ".xpRinvx.f32"
    afreq_path = prefix_abs * ".afreq.f32"

    ncol = nmarkers_all + 1
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
    if nMarkersAll != nmarkers_all
        error("Genotype marker count changed while reading file.")
    end

    missing_value32 = Float32(missing_value)
    maf32 = Float32(MAF)
    keep = trues(nMarkersAll)
    marker_means = Vector{Float32}(undef, nMarkersAll)
    allele_freq = Vector{Float32}(undef, nMarkersAll)
    xpRinvx_all = Vector{Float32}(undef, nMarkersAll)
    missing_masks = Vector{BitVector}(undef, nMarkersAll)

    @showprogress "preparing packed genotypes (dense) ..." for j in 1:nMarkersAll
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
            error("Marker $(marker_ids_all[j]) has only missing values.")
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

    markerID = marker_ids_all[selected]
    means_kept = marker_means[selected]
    afreq_kept = allele_freq[selected]
    xp_kept = xpRinvx_all[selected]
    nMarkers = length(selected)
    stride_bytes = cld(nObs, 4)
    sum2pq = Float64(sum(2.0f0 .* afreq_kept .* (1.0f0 .- afreq_kept)))

    if quality_control
        printstyled("$(nMarkersAll - nMarkers) loci which are fixed or have minor allele frequency < $MAF are removed.\n", bold=true)
    end

    marker_text_bytes = sum(length(id) + 1 for id in markerID)
    obs_text_bytes = sum(length(id) + 1 for id in obsID)
    sidecar_bytes = Int128(marker_text_bytes + obs_text_bytes + 256) + Int128(12) * Int128(nMarkers)
    required_output_bytes = Int128(nMarkers) * Int128(stride_bytes) + sidecar_bytes
    _check_streaming_disk_guard!(output_dir, output_dir, required_output_bytes, Int128(0), disk_guard_ratio)

    packed_row = zeros(UInt8, stride_bytes)
    open(data_path, "w") do io
        @showprogress "writing packed genotypes (dense) ..." for idx in 1:nMarkers
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
    return nothing
end

function _prepare_streaming_genotypes_lowmem!(file::AbstractString,
                                              prefix_abs::AbstractString;
                                              separator,
                                              header::Bool,
                                              missing_value::Real,
                                              quality_control::Bool,
                                              MAF::Real,
                                              center::Bool,
                                              tmpdir,
                                              cleanup_temp::Bool,
                                              disk_guard_ratio::Float64,
                                              marker_ids_all::Vector{String},
                                              nmarkers_all::Int)
    output_dir = dirname(prefix_abs)
    mkpath(output_dir)

    data_path = prefix_abs * ".jgb2"
    meta_path = prefix_abs * ".meta"
    obs_path = prefix_abs * ".obsid.txt"
    marker_path = prefix_abs * ".markerid.txt"
    mean_path = prefix_abs * ".mean.f32"
    xp_path = prefix_abs * ".xpRinvx.f32"
    afreq_path = prefix_abs * ".afreq.f32"
    temp_root = tmpdir === nothing ? output_dir : abspath(string(tmpdir))
    mkpath(temp_root)
    stage_dir = mktempdir(temp_root)
    stage_obs_path = joinpath(stage_dir, "obsid.stage.txt")
    stage_rowmajor_path = joinpath(stage_dir, "rowmajor.stage.jgb2")

    nonce = string(getpid(), "_", string(time_ns()))
    data_tmp = data_path * ".tmp." * nonce
    meta_tmp = meta_path * ".tmp." * nonce
    obs_tmp = obs_path * ".tmp." * nonce
    marker_tmp = marker_path * ".tmp." * nonce
    mean_tmp = mean_path * ".tmp." * nonce
    xp_tmp = xp_path * ".tmp." * nonce
    afreq_tmp = afreq_path * ".tmp." * nonce

    published = false
    try
        # Stage A: scan file, collect per-marker stats, and stream obs IDs.
        stats = _scan_streaming_stats!(file, separator, header, nmarkers_all, marker_ids_all,
                                       missing_value, quality_control, MAF, center,
                                       stage_obs_path)
        nObs = stats.nObs
        selected = stats.selected
        markerID = stats.markerID
        means_kept = stats.means
        afreq_kept = stats.afreq
        xp_kept = stats.xp
        sum2pq = stats.sum2pq
        nMarkers = length(selected)
        stride_bytes = cld(nObs, 4)

        if quality_control
            printstyled("$(nmarkers_all - nMarkers) loci which are fixed or have minor allele frequency < $MAF are removed.\n", bold=true)
        end

        marker_text_bytes = sum(length(id) + 1 for id in markerID)
        obs_text_bytes = filesize(stage_obs_path)
        sidecar_bytes = Int128(marker_text_bytes + obs_text_bytes + 256) + Int128(12) * Int128(nMarkers)
        required_output_bytes = Int128(nMarkers) * Int128(stride_bytes) + sidecar_bytes
        required_temp_bytes = Int128(nObs) * Int128(cld(nMarkers, 4))
        _check_streaming_disk_guard!(output_dir, stage_dir, required_output_bytes, required_temp_bytes, disk_guard_ratio)

        selected_map = zeros(Int32, nmarkers_all)
        @inbounds for (k, j) in enumerate(selected)
            selected_map[j] = Int32(k)
        end

        # Stage B: second pass writes row-major packed spool for retained markers.
        _write_row_major_spool!(file, separator, header, nmarkers_all, selected_map, nMarkers, nObs,
                                missing_value, stage_rowmajor_path)

        # Stage C: out-of-core transpose from row-major spool to marker-major backend payload.
        _transpose_row_major_to_marker_major!(stage_rowmajor_path, data_tmp, nObs, nMarkers)

        # Stage D: write sidecar files and atomically publish all outputs.
        cp(stage_obs_path, obs_tmp; force=true)
        _write_string_lines(marker_tmp, markerID)
        _write_f32_vector(mean_tmp, means_kept)
        _write_f32_vector(xp_tmp, xp_kept)
        _write_f32_vector(afreq_tmp, afreq_kept)

        _write_streaming_manifest(meta_tmp, [
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

        mv(data_tmp, data_path; force=true)
        mv(meta_tmp, meta_path; force=true)
        mv(obs_tmp, obs_path; force=true)
        mv(marker_tmp, marker_path; force=true)
        mv(mean_tmp, mean_path; force=true)
        mv(xp_tmp, xp_path; force=true)
        mv(afreq_tmp, afreq_path; force=true)
        published = true
    finally
        if !published
            _safe_rm(data_tmp)
            _safe_rm(meta_tmp)
            _safe_rm(obs_tmp)
            _safe_rm(marker_tmp)
            _safe_rm(mean_tmp)
            _safe_rm(xp_tmp)
            _safe_rm(afreq_tmp)
        end
        if cleanup_temp
            if isdir(stage_dir)
                rm(stage_dir; recursive=true, force=true)
            end
        end
    end
    return nothing
end

"""
    prepare_streaming_genotypes(file::AbstractString;
                                output_prefix=nothing,
                                separator=',',
                                header=true,
                                missing_value=9.0,
                                quality_control=true,
                                MAF=0.01,
                                center=true,
                                conversion_mode=:lowmem,
                                auto_dense_max_bytes=2^30,
                                tmpdir=nothing,
                                cleanup_temp=true,
                                disk_guard_ratio=0.9)

Convert a dense text genotype file to a marker-major 2-bit packed backend.

Conversion backend selection:
- `conversion_mode=:lowmem` uses out-of-core staged conversion (disk-backed).
- `conversion_mode=:dense` uses in-memory conversion.
- `conversion_mode=:auto` chooses between the two based on `auto_dense_max_bytes`.

Low-memory conversion options:
- `tmpdir`: optional location for temporary conversion files.
- `cleanup_temp`: remove temporary files after successful conversion.
- `disk_guard_ratio`: fail fast when estimated required bytes exceed
  `disk_guard_ratio * available_bytes`.
"""
function prepare_streaming_genotypes(file::AbstractString;
                                     output_prefix=nothing,
                                     separator=',',
                                     header=true,
                                     missing_value=9.0,
                                     quality_control=true,
                                     MAF=0.01,
                                     center=true,
                                     conversion_mode=:lowmem,
                                     auto_dense_max_bytes=2^30,
                                     tmpdir=nothing,
                                     cleanup_temp=true,
                                     disk_guard_ratio=0.9)
    output_prefix = output_prefix === nothing ? splitext(file)[1] * "_stream" : string(output_prefix)
    prefix_abs = abspath(output_prefix)

    marker_ids_all, nmarkers_all = _header_marker_ids(file, separator, header)
    if auto_dense_max_bytes < 0
        error("auto_dense_max_bytes must be non-negative.")
    end
    requested_mode = _normalize_conversion_mode(conversion_mode)
    nobs_for_decision = requested_mode == :auto ? _count_genotype_rows(file, header) : 0
    selected_mode = _choose_conversion_mode(requested_mode, nobs_for_decision, nmarkers_all, Int128(auto_dense_max_bytes))

    if requested_mode == :auto
        dense_bytes_est = Int128(nobs_for_decision) * Int128(nmarkers_all) * Int128(sizeof(Float32))
        printstyled("Auto conversion mode selected :$selected_mode (estimated dense bytes=$dense_bytes_est, auto_dense_max_bytes=$(Int128(auto_dense_max_bytes))).\n",
                    bold=false, color=:light_black)
    end

    if selected_mode == :dense
        _prepare_streaming_genotypes_dense!(file, prefix_abs;
                                            separator=separator,
                                            header=header,
                                            missing_value=missing_value,
                                            quality_control=quality_control,
                                            MAF=MAF,
                                            center=center,
                                            disk_guard_ratio=disk_guard_ratio,
                                            marker_ids_all=marker_ids_all,
                                            nmarkers_all=nmarkers_all)
    else
        _prepare_streaming_genotypes_lowmem!(file, prefix_abs;
                                             separator=separator,
                                             header=header,
                                             missing_value=missing_value,
                                             quality_control=quality_control,
                                             MAF=MAF,
                                             center=center,
                                             tmpdir=tmpdir,
                                             cleanup_temp=cleanup_temp,
                                             disk_guard_ratio=disk_guard_ratio,
                                             marker_ids_all=marker_ids_all,
                                             nmarkers_all=nmarkers_all)
    end

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
