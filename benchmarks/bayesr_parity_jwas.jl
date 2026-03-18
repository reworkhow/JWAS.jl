using CSV
using DataFrames
using Statistics
using JWAS

include(joinpath(@__DIR__, "bayesr_parity_common.jl"))

const DEFAULT_GAMMA = Float64[0.0, 0.01, 0.1, 1.0]
const DEFAULT_START_PI = Float64[0.95, 0.03, 0.015, 0.005]

function env_int(name, default)
    return parse(Int, get(ENV, name, string(default)))
end

function env_float(name, default)
    return parse(Float64, get(ENV, name, string(default)))
end

outdir = get(ENV, "JWAS_PARITY_OUTDIR", joinpath(@__DIR__, "out", "bayesr_parity"))
mode = get(ENV, "JWAS_PARITY_MODE", "fixed_pi")
mode in ("fixed_pi", "estimate_pi") || error("JWAS_PARITY_MODE must be fixed_pi or estimate_pi.")

seed = env_int("JWAS_PARITY_SEED", 2026)
chain_length = env_int("JWAS_PARITY_CHAIN_LENGTH", 80)
burnin = env_int("JWAS_PARITY_BURNIN", 20)
n_obs = env_int("JWAS_PARITY_N_OBS", 40)
n_markers = env_int("JWAS_PARITY_N_MARKERS", 12)
start_h2 = env_float("JWAS_PARITY_START_H2", 0.5)
estimate_pi = mode == "estimate_pi"

dataset = build_bayesr_parity_dataset(seed=seed, n_obs=n_obs, n_markers=n_markers)
vary = var(dataset.y)
start_vare = vary * (1 - start_h2)
start_sigma_sq = vary * start_h2 / (n_markers * sum(DEFAULT_GAMMA .* DEFAULT_START_PI))

datadir = joinpath(outdir, "data")
write_parity_dataset(datadir;
                     ids=dataset.ids,
                     marker_ids=dataset.marker_ids,
                     X=dataset.X,
                     y=dataset.y,
                     gamma=DEFAULT_GAMMA,
                     start_pi=DEFAULT_START_PI,
                     estimate_pi=estimate_pi,
                     chain_length=chain_length,
                     burnin=burnin,
                     start_h2=start_h2,
                     start_sigma_sq=start_sigma_sq,
                     start_vare=start_vare,
                     seed=seed)

geno_path = joinpath(datadir, "genotypes.csv")
pheno_path = joinpath(datadir, "phenotypes.csv")
phenotypes = CSV.read(pheno_path, DataFrame)

global geno = get_genotypes(geno_path, start_sigma_sq,
                            separator=',',
                            method="BayesR",
                            Pi=DEFAULT_START_PI,
                            estimatePi=estimate_pi,
                            G_is_marker_variance=true,
                            estimate_variance=true,
                            estimate_scale=false,
                            quality_control=false,
                            center=false)
model = build_model("y1 = intercept + geno", start_vare)

tmp_output = joinpath(outdir, "tmp_jwas_run")
isdir(tmp_output) && rm(tmp_output; recursive=true, force=true)
output = runMCMC(model, phenotypes,
                 chain_length=chain_length,
                 burnin=burnin,
                 output_samples_frequency=1,
                 output_folder=tmp_output,
                 seed=seed,
                 printout_model_info=false,
                 outputEBV=false,
                 output_heritability=false,
                 fast_blocks=false)

summary_dir = joinpath(outdir, estimate_pi ? "jwas_estimate_pi" : "jwas_fixed_pi")
isdir(summary_dir) && rm(summary_dir; recursive=true, force=true)
write_jwas_parity_summary(output, summary_dir;
                          sigma_sq=model.M[1].meanVara,
                          pi_values=model.M[1].π)

isdir(tmp_output) && rm(tmp_output; recursive=true, force=true)

println("WROTE ", summary_dir)
