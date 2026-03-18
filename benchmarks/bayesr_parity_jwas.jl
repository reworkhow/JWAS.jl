using CSV
using DataFrames
using Distributions
using LinearAlgebra
using Random
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

function env_bool(name, default)
    raw = lowercase(get(ENV, name, default ? "true" : "false"))
    return raw in ("1", "true", "yes", "on")
end

function run_jwas_bayesr_fixed_pi_trace(datadir, geno_path;
                                        chain_length,
                                        gamma,
                                        start_pi)
    initial_state = read_parity_initial_state(datadir)
    geno = get_genotypes(geno_path, initial_state.sigmaSq0,
                         separator=',',
                         method="BayesR",
                         Pi=Float64.(start_pi),
                         estimatePi=false,
                         G_is_marker_variance=true,
                         estimate_variance=true,
                         estimate_scale=false,
                         quality_control=false,
                         center=false)
    T = eltype(geno.genotypes)
    mGibbs = JWAS.GibbsMats(geno.genotypes, ones(eltype(geno.genotypes), geno.nObs))
    geno.mArray = mGibbs.xArray
    geno.mRinvArray = mGibbs.xRinvArray
    geno.mpRinvm = mGibbs.xpRinvx
    geno.α = [T.(initial_state.beta0)]
    geno.δ = [copy(initial_state.delta0)]
    geno.G.val = T(initial_state.sigmaSq0)

    y = T.(CSV.read(joinpath(datadir, "phenotypes.csv"), DataFrame)[!, :y1])
    mu = T(initial_state.mu0)
    vare = T(initial_state.vare0)
    ycorr = T.(y .- mu .- geno.genotypes * geno.α[1])
    nub = 4.0
    nue = 4.0
    scaleb = (nub - 2.0) / nub * initial_state.sigmaSq0
    scalee = (nue - 2.0) / nue * initial_state.vare0

    trace = DataFrame(iter=Int[], sigmaSq=Float64[], ssq=Float64[], nnz=Int[], vare=Float64[])

    for iter in 1:chain_length
        ycorr .+= mu
        rhs_mu = sum(ycorr)
        inv_lhs_mu = 1 / length(y)
        mu_hat = inv_lhs_mu * rhs_mu
        mu = mu_hat + randn() * sqrt(inv_lhs_mu * vare)
        ycorr .-= mu

        JWAS.BayesR!(geno, ycorr, vare)

        ssq, nnz = JWAS.bayesr_sigma_sufficient_statistics(geno.α[1], geno.δ[1], gamma)
        geno.G.val = T((ssq + nub * scaleb) / rand(Chisq(nnz + nub)))
        vare = T((dot(ycorr, ycorr) + nue * scalee) / rand(Chisq(length(y) + nue)))

        push!(trace, (iter=iter, sigmaSq=Float64(geno.G.val), ssq=Float64(ssq), nnz=nnz, vare=Float64(vare)))
    end

    return trace
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
trace_mode = env_bool("JWAS_PARITY_TRACE", false)

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
summary_dir = joinpath(outdir, estimate_pi ? "jwas_estimate_pi" : "jwas_fixed_pi")
isdir(summary_dir) && rm(summary_dir; recursive=true, force=true)

if trace_mode
    mode == "fixed_pi" || error("JWAS parity trace mode currently supports fixed_pi only.")
    burnin == 0 || error("JWAS parity trace mode requires burnin=0.")
    Random.seed!(seed)
    trace = run_jwas_bayesr_fixed_pi_trace(datadir, geno_path;
                                           chain_length=chain_length,
                                           gamma=DEFAULT_GAMMA,
                                           start_pi=DEFAULT_START_PI)
    write_parity_trace(joinpath(summary_dir, "trace_fixed_pi.csv"), trace)
else
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

    write_jwas_parity_summary(output, summary_dir;
                              sigma_sq=model.M[1].meanVara,
                              pi_values=model.M[1].π)

    isdir(tmp_output) && rm(tmp_output; recursive=true, force=true)
end

println("WROTE ", summary_dir)
