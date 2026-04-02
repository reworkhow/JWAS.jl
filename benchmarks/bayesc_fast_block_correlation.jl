using CSV
using DataFrames
using Distributions
using Random
using Statistics
using JWAS

function env_int(name, default)
    return parse(Int, get(ENV, name, string(default)))
end

function env_float(name, default)
    return parse(Float64, get(ENV, name, string(default)))
end

function env_string(name, default)
    return String(get(ENV, name, default))
end

function parse_seed_list(raw::AbstractString)
    isempty(strip(raw)) && error("Seed list must not be empty.")
    return parse.(Int, strip.(split(raw, ",")))
end

function safe_cor(x::AbstractVector, y::AbstractVector)
    length(x) == length(y) || error("Vectors must have the same length.")
    if isempty(x)
        return NaN
    end
    xf = Float64.(x)
    yf = Float64.(y)
    if all(isapprox(v, xf[1]; atol=1e-12) for v in xf) || all(isapprox(v, yf[1]; atol=1e-12) for v in yf)
        return all(isapprox.(xf, yf; atol=1e-12)) ? 1.0 : NaN
    end
    return cor(xf, yf)
end

function joint_probs_from_conditionals(p1::Real, p2::Real, p3::Real)
    probs = Float64[
        1 - p1,
        p1 * (1 - p2),
        p1 * p2 * (1 - p3),
        p1 * p2 * p3,
    ]
    isapprox(sum(probs), 1.0; atol=1e-12) || error("Conditional probabilities must map to a valid class distribution.")
    return probs
end

function build_annotated_dataset(; seed=20260401, n_obs=200, n_markers=1000, h2=0.45)
    rng = MersenneTwister(seed)
    baseline_probs = joint_probs_from_conditionals(0.10, 0.20, 0.20)
    enriched_probs = joint_probs_from_conditionals(0.30, 0.60, 0.60)
    gamma = Float64[0.0, 0.01, 0.1, 1.0]

    ids = ["id_$(i)" for i in 1:n_obs]
    marker_ids = ["m$(j)" for j in 1:n_markers]

    X = Matrix{Float64}(undef, n_obs, n_markers)
    allele_freq = rand(rng, n_markers) .* 0.40 .+ 0.05
    for j in 1:n_markers
        p = allele_freq[j]
        for i in 1:n_obs
            X[i, j] = (rand(rng) < p) + (rand(rng) < p)
        end
    end

    n_enriched = max(1, round(Int, 0.25 * n_markers))
    annotation_1 = zeros(Float64, n_markers)
    annotation_1[randperm(rng, n_markers)[1:n_enriched]] .= 1.0

    n_nuisance = max(1, round(Int, 0.50 * n_markers))
    annotation_2 = zeros(Float64, n_markers)
    annotation_2[randperm(rng, n_markers)[1:n_nuisance]] .= 1.0
    annotations = hcat(annotation_1, annotation_2)

    true_class = Vector{Int}(undef, n_markers)
    for j in 1:n_markers
        probs = annotation_1[j] == 1.0 ? enriched_probs : baseline_probs
        true_class[j] = rand(rng, Categorical(probs))
    end
    if !any(true_class .> 1)
        true_class[1] = 4
    end

    raw_beta = zeros(Float64, n_markers)
    for j in 1:n_markers
        cj = true_class[j]
        if cj > 1
            raw_beta[j] = randn(rng) * sqrt(gamma[cj])
        end
    end

    genetic_raw = X * raw_beta
    scale = sqrt(var(genetic_raw))
    beta_true = raw_beta / scale
    genetic_value = X * beta_true
    residual_sd = sqrt((1 - h2) / h2)
    y = 1.0 .+ genetic_value .+ randn(rng, n_obs) .* residual_sd

    return (
        ids=ids,
        marker_ids=marker_ids,
        X=X,
        y=y,
        annotations=annotations,
        target_h2=h2,
    )
end

function write_dataset(outdir, data)
    mkpath(outdir)
    geno_df = DataFrame(ID=data.ids)
    for (j, marker_id) in enumerate(data.marker_ids)
        geno_df[!, marker_id] = data.X[:, j]
    end
    CSV.write(joinpath(outdir, "genotypes.csv"), geno_df)
    CSV.write(joinpath(outdir, "phenotypes.csv"), DataFrame(ID=data.ids, y1=data.y))
end

function run_bayesc_case(case_outdir;
                         data,
                         seed,
                         chain_length,
                         burnin,
                         output_samples_frequency,
                         fast_blocks,
                         annotated)
    write_dataset(case_outdir, data)
    geno_path = joinpath(case_outdir, "genotypes.csv")
    pheno_path = joinpath(case_outdir, "phenotypes.csv")
    phenotypes = CSV.read(pheno_path, DataFrame)
    vary = var(data.y)
    start_h2 = data.target_h2
    start_vare = vary * (1 - start_h2)
    start_genetic_variance = vary * start_h2

    geno_kwargs = (
        separator=',',
        method="BayesC",
        Pi=0.0,
        estimatePi=true,
        estimate_variance=true,
        estimate_scale=false,
        quality_control=false,
        center=false,
    )

    global geno = if annotated
        get_genotypes(geno_path, start_genetic_variance; annotations=data.annotations, geno_kwargs...)
    else
        get_genotypes(geno_path, start_genetic_variance; geno_kwargs...)
    end

    model = build_model("y1 = intercept + geno", start_vare)
    outputEBV(model, geno.obsID)
    output_dir = joinpath(case_outdir, "run")
    isdir(output_dir) && rm(output_dir; recursive=true, force=true)
    output = runMCMC(
        model,
        phenotypes;
        chain_length=chain_length,
        burnin=burnin,
        output_samples_frequency=output_samples_frequency,
        output_folder=output_dir,
        seed=seed,
        printout_model_info=false,
        outputEBV=true,
        output_heritability=false,
        fast_blocks=fast_blocks,
    )
    return output, model, phenotypes
end

function extract_run_summary(label, output, model, phenotypes)
    marker_effects = output["marker effects geno"][:, [:Marker_ID, :Estimate, :Model_Frequency]]
    ebv = haskey(output, "EBV_y1") ? output["EBV_y1"][:, [:ID, :EBV]] : DataFrame(ID=String[], EBV=Float64[])
    annotation_coefficients = haskey(output, "annotation coefficients geno") ? output["annotation coefficients geno"] : DataFrame()

    return (
        label=label,
        marker_effects=marker_effects,
        ebv=ebv,
        residual_variance=Float64(output["residual variance"][1, :Estimate]),
        pi_state=copy(model.M[1].π),
        annotation_coefficients=annotation_coefficients,
    )
end

function compare_pair(method_label, seed, dense_summary, block_summary)
    dense_markers = sort(dense_summary.marker_effects, :Marker_ID)
    block_markers = sort(block_summary.marker_effects, :Marker_ID)
    dense_ebv = sort(dense_summary.ebv, :ID)
    block_ebv = sort(block_summary.ebv, :ID)

    row = (
        method=method_label,
        seed=seed,
        marker_estimate_correlation=safe_cor(dense_markers.Estimate, block_markers.Estimate),
        model_frequency_correlation=safe_cor(dense_markers.Model_Frequency, block_markers.Model_Frequency),
        ebv_correlation=safe_cor(dense_ebv.EBV, block_ebv.EBV),
        residual_variance_absdiff=abs(dense_summary.residual_variance - block_summary.residual_variance),
        pi_absdiff=if dense_summary.pi_state isa AbstractVector && block_summary.pi_state isa AbstractVector &&
                      length(dense_summary.pi_state) == 1 && length(block_summary.pi_state) == 1
            abs(Float64(dense_summary.pi_state[1]) - Float64(block_summary.pi_state[1]))
        elseif !(dense_summary.pi_state isa AbstractVector) && !(block_summary.pi_state isa AbstractVector)
            abs(Float64(dense_summary.pi_state) - Float64(block_summary.pi_state))
        else
            NaN
        end,
        pi_vector_correlation=if dense_summary.pi_state isa AbstractVector && block_summary.pi_state isa AbstractVector &&
                                 length(dense_summary.pi_state) > 1
            safe_cor(dense_summary.pi_state, block_summary.pi_state)
        else
            NaN
        end,
        annotation_coefficient_correlation=if !isempty(dense_summary.annotation_coefficients) && !isempty(block_summary.annotation_coefficients)
            safe_cor(dense_summary.annotation_coefficients.Estimate, block_summary.annotation_coefficients.Estimate)
        else
            NaN
        end,
        annotation_coefficient_max_absdiff=if !isempty(dense_summary.annotation_coefficients) && !isempty(block_summary.annotation_coefficients)
            maximum(abs.(Float64.(dense_summary.annotation_coefficients.Estimate) .- Float64.(block_summary.annotation_coefficients.Estimate)))
        else
            NaN
        end,
    )
    return row
end

function summarize_pairs(rows::DataFrame)
    out = NamedTuple[]
    for method in unique(rows.method)
        sub = rows[rows.method .== method, :]
        for metric in names(sub)
            metric in ("method", "seed") && continue
            vals = Float64.(sub[!, metric])
            if all(isnan, vals)
                continue
            end
            push!(out, (
                method=String(method),
                metric=String(metric),
                mean_value=mean(skipmissing(vals[.!isnan.(vals)])),
                min_value=minimum(vals[.!isnan.(vals)]),
                max_value=maximum(vals[.!isnan.(vals)]),
            ))
        end
    end
    return DataFrame(out)
end

function main(argv=ARGS)
    length(argv) > 1 && error("Usage: julia --project=. --startup-file=no benchmarks/bayesc_fast_block_correlation.jl [outdir]")

    outdir = length(argv) == 1 ? abspath(argv[1]) : joinpath(@__DIR__, "out", "bayesc_fast_block_correlation")
    data_seed = env_int("JWAS_BAYESC_CORR_DATA_SEED", 20260401)
    chain_length = env_int("JWAS_BAYESC_CORR_CHAIN_LENGTH", 4000)
    burnin = env_int("JWAS_BAYESC_CORR_BURNIN", 800)
    output_samples_frequency = env_int("JWAS_BAYESC_CORR_OUTPUT_FREQ", 50)
    n_obs = env_int("JWAS_BAYESC_CORR_N_OBS", 200)
    n_markers = env_int("JWAS_BAYESC_CORR_N_MARKERS", 1000)
    h2 = env_float("JWAS_BAYESC_CORR_H2", 0.45)
    seeds = parse_seed_list(env_string("JWAS_BAYESC_CORR_SEEDS", "2026,2027,2028,2029,2030"))

    mkpath(outdir)
    data = build_annotated_dataset(seed=data_seed, n_obs=n_obs, n_markers=n_markers, h2=h2)
    pair_rows = NamedTuple[]

    for seed in seeds
        dense_out, dense_model, phenotypes = run_bayesc_case(joinpath(outdir, "BayesC_dense_seed_$(seed)");
                                                             data=data, seed=seed, chain_length=chain_length, burnin=burnin,
                                                             output_samples_frequency=output_samples_frequency,
                                                             fast_blocks=false, annotated=false)
        block_out, block_model, _ = run_bayesc_case(joinpath(outdir, "BayesC_fast_blocks_1_seed_$(seed)");
                                                    data=data, seed=seed, chain_length=chain_length, burnin=burnin,
                                                    output_samples_frequency=output_samples_frequency,
                                                    fast_blocks=1, annotated=false)
        push!(pair_rows, compare_pair("BayesC", seed,
                                      extract_run_summary("BayesC_dense", dense_out, dense_model, phenotypes),
                                      extract_run_summary("BayesC_fast_blocks_1", block_out, block_model, phenotypes)))

        adense_out, adense_model, aphenotypes = run_bayesc_case(joinpath(outdir, "Annotated_BayesC_dense_seed_$(seed)");
                                                                data=data, seed=seed, chain_length=chain_length, burnin=burnin,
                                                                output_samples_frequency=output_samples_frequency,
                                                                fast_blocks=false, annotated=true)
        ablock_out, ablock_model, _ = run_bayesc_case(joinpath(outdir, "Annotated_BayesC_fast_blocks_1_seed_$(seed)");
                                                      data=data, seed=seed, chain_length=chain_length, burnin=burnin,
                                                      output_samples_frequency=output_samples_frequency,
                                                      fast_blocks=1, annotated=true)
        push!(pair_rows, compare_pair("Annotated_BayesC", seed,
                                      extract_run_summary("Annotated_BayesC_dense", adense_out, adense_model, aphenotypes),
                                      extract_run_summary("Annotated_BayesC_fast_blocks_1", ablock_out, ablock_model, aphenotypes)))
    end

    pair_df = sort!(DataFrame(pair_rows), [:method, :seed])
    summary_df = sort!(summarize_pairs(pair_df), [:method, :metric])
    CSV.write(joinpath(outdir, "pairwise_correlations.csv"), pair_df)
    CSV.write(joinpath(outdir, "summary.csv"), summary_df)
    println("WROTE ", joinpath(outdir, "pairwise_correlations.csv"))
    println("WROTE ", joinpath(outdir, "summary.csv"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
