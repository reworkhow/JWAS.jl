using CSV
using DataFrames
using Statistics
using JWAS
using JWAS.Datasets

const SIMULATED_ANNOTATIONS_START_PI = Float64[0.99, 0.006, 0.003, 0.001]

function parse_seed_list(raw::AbstractString)
    isempty(strip(raw)) && error("Seed list must not be empty.")
    return parse.(Int, strip.(split(raw, ",")))
end

env_int(name, default) = parse(Int, get(ENV, name, string(default)))
env_float(name, default) = parse(Float64, get(ENV, name, string(default)))
env_string(name, default) = String(get(ENV, name, default))

safe_cor(x, y) = (length(x) == 0 || length(y) == 0 || std(x) == 0 || std(y) == 0) ? NaN : cor(Float64.(x), Float64.(y))
to_float(x) = x isa Number ? Float64(x) : parse(Float64, strip(replace(String(x), "\"" => "")))

function annotation_vector(output_dir::AbstractString)
    path = joinpath(output_dir, "annotation_coefficients_geno.txt")
    isfile(path) || return Float64[]
    coeff_df = CSV.read(path, DataFrame)
    steps = if "Step" in names(coeff_df)
        String.(coeff_df[!, :Step])
    else
        fill("step1_zero_vs_nonzero", nrow(coeff_df))
    end
    order = sortperm(collect(zip(String.(coeff_df[!, :Annotation]), steps)))
    return to_float.(coeff_df[order, :Estimate])
end

function pi_values(output_dir::AbstractString)
    path = joinpath(output_dir, "pi_geno.txt")
    isfile(path) || return (Float64[], NaN)
    pi_df = CSV.read(path, DataFrame)
    estimates = to_float.(pi_df[!, :Estimate])
    if nrow(pi_df) == 1
        return (Float64[], estimates[1])
    end
    order = sortperm(to_float.(pi_df[!, :π]))
    return (estimates[order], NaN)
end

function ebv_table(output_dir::AbstractString)
    path = joinpath(output_dir, "EBV_y1.txt")
    isfile(path) || return DataFrame()
    ebv_df = CSV.read(path, DataFrame)[:, [:ID, :EBV]]
    sort!(ebv_df, :ID)
    return ebv_df
end

function marker_table(output_dir::AbstractString)
    marker_df = CSV.read(joinpath(output_dir, "marker_effects_geno.txt"), DataFrame)[:, [:Marker_ID, :Estimate, :Model_Frequency]]
    sort!(marker_df, :Marker_ID)
    return marker_df
end

function residual_variance_value(output_dir::AbstractString)
    resid_df = CSV.read(joinpath(output_dir, "residual_variance.txt"), DataFrame)
    return to_float(resid_df[1, :Estimate])
end

function dataset_paths()
    return (
        genotypes=Datasets.dataset("genotypes.csv", dataset_name="simulated_annotations"),
        phenotypes=Datasets.dataset("phenotypes.csv", dataset_name="simulated_annotations"),
        annotations=Datasets.dataset("annotations.csv", dataset_name="simulated_annotations"),
        truth=Datasets.dataset("truth.csv", dataset_name="simulated_annotations"),
    )
end

function run_one(case::NamedTuple, seed::Integer, paths;
                 chain_length::Integer, burnin::Integer, output_samples_frequency::Integer, start_h2::Real,
                 output_root::AbstractString)
    pheno_df = CSV.read(paths.phenotypes, DataFrame)
    geno_ids = CSV.read(paths.genotypes, DataFrame; select=["ID"])
    annotations_df = CSV.read(paths.annotations, DataFrame)
    annotations = Matrix{Float64}(annotations_df[:, [:functional, :random_anno]])

    vary = var(Float64.(pheno_df.y1))
    start_genetic_variance = vary * start_h2
    start_vare = vary * (1 - start_h2)

    kwargs = (
        separator=',',
        method=case.method,
        estimatePi=true,
        estimate_variance=true,
        estimate_scale=false,
        quality_control=false,
        center=false,
    )

    if case.method == "BayesC"
        global geno = if case.annotated
            get_genotypes(
                paths.genotypes, start_genetic_variance;
                kwargs...,
                Pi=0.0,
                annotations=annotations,
            )
        else
            get_genotypes(
                paths.genotypes, start_genetic_variance;
                kwargs...,
                Pi=0.0,
            )
        end
    elseif case.method == "BayesR"
        global geno = if case.annotated
            get_genotypes(
                paths.genotypes, start_genetic_variance;
                kwargs...,
                Pi=copy(SIMULATED_ANNOTATIONS_START_PI),
                G_is_marker_variance=false,
                annotations=annotations,
            )
        else
            get_genotypes(
                paths.genotypes, start_genetic_variance;
                kwargs...,
                Pi=copy(SIMULATED_ANNOTATIONS_START_PI),
                G_is_marker_variance=false,
            )
        end
    else
        error("Unsupported method $(case.method)")
    end

    model = build_model("y1 = intercept + geno", start_vare)
    outputEBV(model, geno_ids.ID)

    seed_dir = joinpath(output_root, case.variant, "seed_$(seed)")
    mkpath(seed_dir)
    output_dir = joinpath(seed_dir, "run")
    isdir(output_dir) && rm(output_dir; recursive=true, force=true)
    runMCMC(
        model,
        pheno_df;
        chain_length=chain_length,
        burnin=burnin,
        output_samples_frequency=output_samples_frequency,
        output_folder=output_dir,
        seed=seed,
        outputEBV=true,
        output_heritability=false,
        printout_model_info=false,
        printout_frequency=chain_length + 1,
        fast_blocks=case.fast_blocks,
    )

    pi_vec, pi_scalar_val = pi_values(output_dir)

    return (
        marker_df=marker_table(output_dir),
        annotation_vec=annotation_vector(output_dir),
        pi_vec=pi_vec,
        pi_scalar=pi_scalar_val,
        residual_variance=residual_variance_value(output_dir),
        marker_variance=Float64(model.M[1].meanVara),
        ebv_df=ebv_table(output_dir),
    )
end

function summarize_pair(case::NamedTuple, run_a, run_b)
    marker_compare = innerjoin(
        rename(copy(run_a.marker_df), :Estimate => :Estimate_seed_a, :Model_Frequency => :Model_Frequency_seed_a),
        rename(copy(run_b.marker_df), :Estimate => :Estimate_seed_b, :Model_Frequency => :Model_Frequency_seed_b),
        on=:Marker_ID,
    )

    ebv_corr = NaN
    if nrow(run_a.ebv_df) > 0 && nrow(run_b.ebv_df) > 0
        ebv_compare = innerjoin(
            rename(copy(run_a.ebv_df), :EBV => :EBV_seed_a),
            rename(copy(run_b.ebv_df), :EBV => :EBV_seed_b),
            on=:ID,
        )
        ebv_corr = safe_cor(ebv_compare.EBV_seed_a, ebv_compare.EBV_seed_b)
    end

    return (
        method=case.method,
        variant=case.variant,
        annotated=case.annotated,
        fast_blocks=string(case.fast_blocks),
        marker_corr=safe_cor(marker_compare.Estimate_seed_a, marker_compare.Estimate_seed_b),
        pip_corr=safe_cor(marker_compare.Model_Frequency_seed_a, marker_compare.Model_Frequency_seed_b),
        ebv_corr=ebv_corr,
        annotation_coeff_corr=safe_cor(run_a.annotation_vec, run_b.annotation_vec),
        pi_vector_corr=safe_cor(run_a.pi_vec, run_b.pi_vec),
        pi_scalar_abs_diff=abs(run_a.pi_scalar - run_b.pi_scalar),
        residual_variance_abs_diff=abs(run_a.residual_variance - run_b.residual_variance),
        marker_variance_abs_diff=abs(run_a.marker_variance - run_b.marker_variance),
    )
end

function default_output_dir()
    return joinpath(@__DIR__, "out", "simulated_annotations_method_matrix")
end

function main(args=ARGS)
    output_dir = length(args) == 1 ? abspath(args[1]) :
                 length(args) == 0 ? default_output_dir() :
                 error("Usage: julia --project=. --startup-file=no benchmarks/simulated_annotations_method_matrix.jl [OUTPUT_DIR]")
    mkpath(output_dir)

    seeds = parse_seed_list(env_string("JWAS_SIMULATED_ANNOTATIONS_SEEDS", "100,110"))
    length(seeds) == 2 || error("JWAS_SIMULATED_ANNOTATIONS_SEEDS must contain exactly two seeds.")

    chain_length = env_int("JWAS_SIMULATED_ANNOTATIONS_CHAIN_LENGTH", 5000)
    burnin = env_int("JWAS_SIMULATED_ANNOTATIONS_BURNIN", 1000)
    output_freq = env_int("JWAS_SIMULATED_ANNOTATIONS_OUTPUT_FREQ", 10)
    start_h2 = env_float("JWAS_SIMULATED_ANNOTATIONS_START_H2", 0.5)

    cases = [
        (method="BayesC", variant="BayesC_dense", annotated=false, fast_blocks=false),
        (method="BayesC", variant="BayesC_fast_blocks_1", annotated=false, fast_blocks=1),
        (method="BayesC", variant="Annotated_BayesC_dense", annotated=true, fast_blocks=false),
        (method="BayesC", variant="Annotated_BayesC_fast_blocks_1", annotated=true, fast_blocks=1),
        (method="BayesR", variant="BayesR_dense", annotated=false, fast_blocks=false),
        (method="BayesR", variant="BayesR_fast_blocks_1", annotated=false, fast_blocks=1),
        (method="BayesR", variant="Annotated_BayesR_dense", annotated=true, fast_blocks=false),
        (method="BayesR", variant="Annotated_BayesR_fast_blocks_1", annotated=true, fast_blocks=1),
    ]

    paths = dataset_paths()
    summary_rows = NamedTuple[]
    for case in cases
        run_a = run_one(case, seeds[1], paths;
                        chain_length=chain_length,
                        burnin=burnin,
                        output_samples_frequency=output_freq,
                        start_h2=start_h2,
                        output_root=output_dir)
        run_b = run_one(case, seeds[2], paths;
                        chain_length=chain_length,
                        burnin=burnin,
                        output_samples_frequency=output_freq,
                        start_h2=start_h2,
                        output_root=output_dir)
        push!(summary_rows, summarize_pair(case, run_a, run_b))
    end

    summary_df = DataFrame(summary_rows)
    CSV.write(joinpath(output_dir, "cross_seed_summary.csv"), summary_df)
    println("Wrote simulated_annotations seed matrix summary to ", joinpath(output_dir, "cross_seed_summary.csv"))
end

main()
