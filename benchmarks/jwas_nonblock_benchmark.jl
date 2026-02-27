using LinearAlgebra
using Random
using DataFrames
using Printf
using Dates
using Statistics

using JWAS

const N = parse(Int, get(ENV, "BENCH_N", "50000"))
const P_LIST = parse.(Int, split(get(ENV, "BENCH_P_LIST", "100000,200000"), ","))
const L_LIST = parse.(Int, split(get(ENV, "BENCH_L_LIST", "2,4"), ","))
const TARGET_P = parse(Int, get(ENV, "TARGET_P", "2000000"))
const TARGET_L = parse(Int, get(ENV, "TARGET_L", "2000"))
const REPS = parse(Int, get(ENV, "BENCH_REPS", "2"))
const WARMUP_N = parse(Int, get(ENV, "WARMUP_N", "2000"))
const WARMUP_P = parse(Int, get(ENV, "WARMUP_P", "2000"))
const WARMUP_L = parse(Int, get(ENV, "WARMUP_L", "5"))
const OUT_TXT = get(ENV, "BENCH_TXT", "/home/qtlcheng/jwas_nonblock_benchmark_results.txt")

cpus = try parse(Int, get(ENV, "SLURM_CPUS_PER_TASK", "1")) catch; 1 end
BLAS.set_num_threads(max(cpus, 1))

function logline(msg)
    println("[", Dates.format(now(), "HH:MM:SS"), "] ", msg)
    flush(stdout)
end

function run_case(n::Int, p::Int, l::Int; seed::Int=2026, do_log::Bool=true, cleanup::Bool=true)
    if do_log
        logline(@sprintf("CASE_START N=%d P=%d L=%d seed=%d", n, p, l, seed))
    end

    Random.seed!(seed)
    pheno = DataFrame(ID = string.(1:n), y1 = randn(Float32, n))

    t_x = time()
    X = rand(Float32, n, p)
    t_x = time() - t_x
    if do_log
        logline(@sprintf("X_GENERATED N=%d P=%d t=%.3f", n, p, t_x))
    end

    t_geno = time()
    global geno = get_genotypes(X, 1.0;
                                method="BayesC",
                                quality_control=false,
                                center=false,
                                estimatePi=false,
                                estimate_variance=false,
                                estimate_scale=false)
    t_geno = time() - t_geno
    if do_log
        logline(@sprintf("GENO_BUILT P=%d t=%.3f", p, t_geno))
    end

    X = nothing
    GC.gc()

    model = build_model("y1 = intercept + geno", 1.0)
    outdir = "/home/qtlcheng/jwas_nonblock_bench_N$(n)_P$(p)_L$(l)_$(Dates.format(now(), "yyyymmdd_HHMMSS"))"

    burnin = max(1, Int(floor(0.1 * l)))

    t_mcmc = time()
    runMCMC(model, pheno;
            chain_length=l,
            burnin=burnin,
            output_samples_frequency=l + 1,
            output_folder=outdir,
            seed=seed,
            printout_model_info=false,
            printout_frequency=l + 1,
            outputEBV=false,
            output_heritability=false,
            fast_blocks=false,
            memory_guard=:warn)
    t_mcmc = time() - t_mcmc

    outer = model.MCMCinfo.chain_length

    if cleanup
        rm(outdir; recursive=true, force=true)
    end
    model = nothing
    pheno = nothing
    global geno = nothing
    GC.gc()

    if do_log
        logline(@sprintf("CASE_DONE N=%d P=%d L=%d outer=%d t_mcmc=%.3f", n, p, l, outer, t_mcmc))
    end

    return Dict(
        "N" => n,
        "P" => p,
        "L_input" => l,
        "outer" => outer,
        "t_x" => t_x,
        "t_geno" => t_geno,
        "t_mcmc" => t_mcmc,
        "seed" => seed,
    )
end

function fit_line(xs::Vector{Float64}, ys::Vector{Float64})
    xbar = mean(xs)
    ybar = mean(ys)
    denom = sum((x - xbar)^2 for x in xs)
    if denom == 0.0
        return (ybar, 0.0)
    end
    b = sum((x - xbar) * (y - ybar) for (x, y) in zip(xs, ys)) / denom
    a = ybar - b * xbar
    return (a, b)
end

logline("WARMUP_START")
_ = run_case(WARMUP_N, WARMUP_P, WARMUP_L; seed=123456, do_log=false, cleanup=true)
logline("WARMUP_DONE")

cases = Any[]
for rep in 1:REPS
    for p in P_LIST
        for l in L_LIST
            seed = 3000000 + rep * 10000 + p + l
            c = run_case(N, p, l; seed=seed, do_log=true, cleanup=true)
            c["rep"] = rep
            push!(cases, c)
        end
    end
end

per_p = Dict{Int,Any}()
for p in P_LIST
    subset = [c for c in cases if c["P"] == p]
    if length(subset) < 2
        continue
    end

    xs = Float64[c["outer"] for c in subset]
    ys = Float64[c["t_mcmc"] for c in subset]
    a, b = fit_line(xs, ys)

    target_t_mcmc = a + b * TARGET_L

    per_p[p] = Dict(
        "n_points" => length(subset),
        "intercept" => a,
        "per_iter" => b,
        "target_t_mcmc" => target_t_mcmc,
    )

    logline(@sprintf("DECOMP P=%d points=%d intercept=%.3f per_iter=%.3f target_L=%d target_t=%.3f",
                     p, length(subset), a, b, TARGET_L, target_t_mcmc))
end

sorted_p = sort(collect(keys(per_p)))
est_seconds = NaN
est_hours = NaN
if length(sorted_p) >= 2
    xs = Float64[float(p) for p in sorted_p]
    ys = Float64[per_p[p]["target_t_mcmc"] for p in sorted_p]
    a, b = fit_line(xs, ys)
    est_seconds = a + b * TARGET_P
    est_hours = est_seconds / 3600

    logline(@sprintf("EXTRAPOLATED using P=%s => target P=%d L=%d t_mcmc=%.3f sec (%.3f h)",
                     join(sorted_p, ","), TARGET_P, TARGET_L, est_seconds, est_hours))
end

open(OUT_TXT, "w") do io
    println(io, "timestamp=", now())
    println(io, "N=", N)
    println(io, "P_LIST=", join(P_LIST, ","))
    println(io, "L_LIST=", join(L_LIST, ","))
    println(io, "REPS=", REPS)
    println(io, "WARMUP_N=", WARMUP_N)
    println(io, "WARMUP_P=", WARMUP_P)
    println(io, "WARMUP_L=", WARMUP_L)
    println(io, "TARGET_P=", TARGET_P)
    println(io, "TARGET_L=", TARGET_L)
    println(io, "BLAS_THREADS=", BLAS.get_num_threads())

    for c in cases
        @printf(io, "CASE rep=%d N=%d P=%d L=%d outer=%d t_x=%.6f t_geno=%.6f t_mcmc=%.6f seed=%d\n",
                c["rep"], c["N"], c["P"], c["L_input"], c["outer"], c["t_x"], c["t_geno"], c["t_mcmc"], c["seed"])
    end

    for p in sorted_p
        d = per_p[p]
        @printf(io, "DECOMP P=%d points=%d intercept=%.6f per_iter=%.6f target_t=%.6f\n",
                p, d["n_points"], d["intercept"], d["per_iter"], d["target_t_mcmc"])
    end

    if isfinite(est_seconds)
        @printf(io, "SUMMARY_EST target_seconds=%.6f target_hours=%.6f\n", est_seconds, est_hours)
    end
end

logline("RESULTS_WRITTEN $(OUT_TXT)")
if isfinite(est_seconds)
    println(@sprintf("SUMMARY_EST target_seconds=%.3f target_hours=%.3f", est_seconds, est_hours))
end
