using CSV
using DataFrames
using Distributions
using LinearAlgebra
using Statistics
using JWAS

include(joinpath(@__DIR__, "bayesr_parity_common.jl"))

const DEFAULT_GAMMA = Float64[0.0, 0.01, 0.1, 1.0]
const DEFAULT_START_PI = Float64[0.95, 0.03, 0.015, 0.005]
const DEFAULT_REPLAY_SEED = 20260318

function env_int(name, default)
    parse(Int, get(ENV, name, string(default)))
end

function env_float(name, default)
    parse(Float64, get(ENV, name, string(default)))
end

function parse_config(datadir)
    cfg = CSV.read(joinpath(datadir, "config.csv"), DataFrame)
    Dict(String(row.key) => String(row.value) for row in eachrow(cfg))
end

function parse_numeric_vector(raw)
    parse.(Float64, split(raw, ","))
end

function draw_value(draws, kind::AbstractString, idx::Integer=1)
    matches = draws[(draws.kind .== kind) .& (draws.index .== idx), :]
    nrow(matches) == 1 || error("Expected exactly one replay draw for $(kind)[$idx].")
    Float64(matches.value[1])
end

function choose_class(probs, u_class)
    cumprob = cumsum(probs)
    idx = findfirst(>=(u_class), cumprob)
    return idx === nothing ? length(probs) : idx
end

function bayesr_class_probabilities(rhs, xpx, vare, sigmaSq, pi, gamma)
    nclasses = length(gamma)
    log_probs = zeros(Float64, nclasses)
    log_probs[1] = log(pi[1])
    invVarRes = 1 / vare
    for k in 2:nclasses
        varEffect = gamma[k] * sigmaSq
        invVarEffect = 1 / varEffect
        lhs = xpx * invVarRes + invVarEffect
        invLhs = 1 / lhs
        betaHat = invLhs * rhs
        log_probs[k] = 0.5 * (log(invLhs) - log(varEffect) + betaHat * rhs) + log(pi[k])
    end
    max_log = maximum(log_probs)
    probs = exp.(log_probs .- max_log)
    probs ./ sum(probs)
end

function ensure_replay_dataset(outdir)
    datadir = joinpath(outdir, "data")
    if !isfile(joinpath(datadir, "config.csv"))
        seed = env_int("JWAS_PARITY_SEED", 2026)
        n_obs = env_int("JWAS_PARITY_N_OBS", 40)
        n_markers = env_int("JWAS_PARITY_N_MARKERS", 12)
        start_h2 = env_float("JWAS_PARITY_START_H2", 0.5)
        dataset = build_bayesr_parity_dataset(seed=seed, n_obs=n_obs, n_markers=n_markers)
        vary = var(dataset.y)
        start_vare = vary * (1 - start_h2)
        start_sigma_sq = vary * start_h2 / (n_markers * sum(DEFAULT_GAMMA .* DEFAULT_START_PI))
        write_parity_dataset(datadir;
                             ids=dataset.ids,
                             marker_ids=dataset.marker_ids,
                             X=dataset.X,
                             y=dataset.y,
                             gamma=DEFAULT_GAMMA,
                             start_pi=DEFAULT_START_PI,
                             estimate_pi=false,
                             chain_length=1,
                             burnin=0,
                             start_h2=start_h2,
                             start_sigma_sq=start_sigma_sq,
                             start_vare=start_vare,
                             seed=seed)
    end

    if !isfile(joinpath(datadir, "replay_draws_iteration1.csv"))
        cfg = parse_config(datadir)
        gamma = parse_numeric_vector(cfg["gamma"])
        draws = build_bayesr_replay_draws(length(gamma) == 0 ? 0 : ncol(CSV.read(joinpath(datadir, "genotypes.csv"), DataFrame)) - 1;
                                          seed=DEFAULT_REPLAY_SEED)
        write_bayesr_replay_draws(joinpath(datadir, "replay_draws_iteration1.csv"), draws)
    end

    datadir
end

function run_replay(outdir)
    datadir = ensure_replay_dataset(outdir)
    cfg = parse_config(datadir)
    gamma = parse_numeric_vector(cfg["gamma"])
    start_pi = parse_numeric_vector(cfg["start_pi"])
    initial_state = read_parity_initial_state(datadir)
    draws = read_bayesr_replay_draws(joinpath(datadir, "replay_draws_iteration1.csv"))

    geno_path = joinpath(datadir, "genotypes.csv")
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
    mGibbs = JWAS.GibbsMats(geno.genotypes, ones(T, geno.nObs))
    geno.mArray = mGibbs.xArray
    geno.mRinvArray = mGibbs.xRinvArray
    geno.mpRinvm = mGibbs.xpRinvx
    geno.α = [T.(initial_state.beta0)]
    geno.δ = [copy(initial_state.delta0)]
    geno.G.val = T(initial_state.sigmaSq0)

    y = T.(CSV.read(joinpath(datadir, "phenotypes.csv"), DataFrame)[!, :y1])
    mu_old = T(initial_state.mu0)
    vare_old = T(initial_state.vare0)
    ycorr = T.(y .- mu_old .- geno.genotypes * geno.α[1])
    marker_ids = String.(CSV.read(joinpath(datadir, "initial_state.csv"), DataFrame)[!, :marker_id])

    ycorr .+= mu_old
    rhs_mu = sum(ycorr)
    inv_lhs_mu = one(T) / length(y)
    mu_hat = inv_lhs_mu * rhs_mu
    z_mu = T(draw_value(draws, "mu_normal", 1))
    mu_new = mu_hat + z_mu * sqrt(inv_lhs_mu * vare_old)
    ycorr .-= mu_new

    marker_rows = NamedTuple[]
    invVarRes = one(T) / vare_old

    for j in eachindex(geno.α[1])
        x = geno.mArray[j]
        xRinv = geno.mRinvArray[j]
        oldAlpha = geno.α[1][j]
        rhs = (dot(xRinv, ycorr) + geno.mpRinvm[j] * oldAlpha) * invVarRes
        probs = bayesr_class_probabilities(Float64(rhs), Float64(geno.mpRinvm[j]), Float64(vare_old), Float64(geno.G.val), start_pi, gamma)
        u_class = draw_value(draws, "marker_class_uniform", j)
        chosen_class = choose_class(probs, u_class)
        geno.δ[1][j] = chosen_class

        beta_hat_chosen = 0.0
        inv_lhs_chosen = 0.0
        z_beta = draw_value(draws, "marker_beta_normal", j)

        if chosen_class == 1
            if oldAlpha != 0
                BLAS.axpy!(oldAlpha, x, ycorr)
            end
            geno.α[1][j] = zero(T)
        else
            varEffect = gamma[chosen_class] * Float64(geno.G.val)
            invVarEffect = 1 / varEffect
            lhs = Float64(geno.mpRinvm[j]) * Float64(invVarRes) + invVarEffect
            inv_lhs_chosen = 1 / lhs
            beta_hat_chosen = inv_lhs_chosen * Float64(rhs)
            geno.α[1][j] = T(beta_hat_chosen + z_beta * sqrt(inv_lhs_chosen))
            BLAS.axpy!(oldAlpha - geno.α[1][j], x, ycorr)
        end

        push!(marker_rows, (
            marker_id=marker_ids[j],
            rhs=Float64(rhs),
            old_alpha=Float64(oldAlpha),
            p_class1=probs[1],
            p_class2=probs[2],
            p_class3=probs[3],
            p_class4=probs[4],
            u_class=u_class,
            chosen_class=chosen_class,
            beta_hat_chosen=beta_hat_chosen,
            inv_lhs_chosen=inv_lhs_chosen,
            z_beta=z_beta,
            new_alpha=Float64(geno.α[1][j]),
            ycorr_norm_after=Float64(norm(ycorr)),
        ))
    end

    ssq, nnz = JWAS.bayesr_sigma_sufficient_statistics(geno.α[1], geno.δ[1], gamma)
    nub = 4.0
    nue = 4.0
    scaleb = (nub - 2.0) / nub * initial_state.sigmaSq0
    scalee = (nue - 2.0) / nue * initial_state.vare0
    chisq_sigma = draw_value(draws, "sigma_chisq", 1)
    chisq_vare = draw_value(draws, "vare_chisq", 1)
    sigmaSq_new = (Float64(ssq) + nub * scaleb) / chisq_sigma
    vare_new = (Float64(dot(ycorr, ycorr)) + nue * scalee) / chisq_vare

    marker_df = DataFrame(marker_rows)
    scalar_df = DataFrame(
        field=[
            "mu_old",
            "mu_hat",
            "z_mu",
            "mu_new",
            "sigmaSq_old",
            "ssq",
            "nnz",
            "chisq_sigma",
            "sigmaSq_new",
            "vare_old",
            "chisq_vare",
            "vare_new",
        ],
        value=[
            Float64(mu_old),
            Float64(mu_hat),
            Float64(z_mu),
            Float64(mu_new),
            Float64(geno.G.val),
            Float64(ssq),
            Float64(nnz),
            chisq_sigma,
            sigmaSq_new,
            Float64(vare_old),
            chisq_vare,
            vare_new,
        ],
    )

    summary_dir = joinpath(outdir, "jwas_fixed_pi")
    mkpath(summary_dir)
    CSV.write(joinpath(summary_dir, "replay_marker_iteration1.csv"), marker_df)
    CSV.write(joinpath(summary_dir, "replay_scalar_iteration1.csv"), scalar_df)
    println("WROTE ", summary_dir)
end

length(ARGS) == 1 || error("Usage: julia --project=. benchmarks/bayesr_parity_replay_jwas.jl <outdir>")
run_replay(ARGS[1])
