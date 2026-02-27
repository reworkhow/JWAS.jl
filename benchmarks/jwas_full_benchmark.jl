using LinearAlgebra
using Random
using DataFrames
using Printf
using Dates

# Use project environment already activated via --project in launcher.
using JWAS

const N = parse(Int, get(ENV, "BENCH_N", "50000"))
const P_LIST = parse.(Int, split(get(ENV, "BENCH_P_LIST", "100000,200000"), ","))
const L_LIST = parse.(Int, split(get(ENV, "BENCH_L_LIST", "223,446"), ","))
const TARGET_P = parse(Int, get(ENV, "TARGET_P", "2000000"))
const TARGET_L = parse(Int, get(ENV, "TARGET_L", "2000"))
const OUT_TXT = get(ENV, "BENCH_TXT", "/home/qtlcheng/jwas_full_benchmark_results.txt")

cpus = try parse(Int, get(ENV, "SLURM_CPUS_PER_TASK", "1")) catch; 1 end
BLAS.set_num_threads(max(cpus, 1))

function logline(msg)
    println("[", Dates.format(now(), "HH:MM:SS"), "] ", msg)
    flush(stdout)
end

function run_case(n::Int, p::Int, l::Int; seed::Int=2026)
    logline(@sprintf("CASE_START N=%d P=%d L=%d", n, p, l))
    Random.seed!(seed)

    # Synthetic phenotype for timed MCMC run.
    pheno = DataFrame(ID = string.(1:n), y1 = randn(Float32, n))

    t_x = time()
    X = rand(Float32, n, p)
    t_x = time() - t_x
    logline(@sprintf("X_GENERATED N=%d P=%d t=%.3f", n, p, t_x))

    t_geno = time()
    global geno = get_genotypes(X, 1.0;
                                method="BayesC",
                                quality_control=false,
                                center=false,
                                estimatePi=false,
                                estimate_variance=false,
                                estimate_scale=false)
    t_geno = time() - t_geno
    logline(@sprintf("GENO_BUILT P=%d t=%.3f", p, t_geno))

    # Free source matrix after genotype object is created.
    X = nothing
    GC.gc()

    model = build_model("y1 = intercept + geno", 1.0)
    outdir = "/home/qtlcheng/jwas_full_bench_N$(n)_P$(p)_L$(l)_$(Dates.format(now(), "yyyymmdd_HHMMSS"))"

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
            fast_blocks=true,
            memory_guard=:warn)
    t_mcmc = time() - t_mcmc

    outer = model.MCMCinfo.chain_length
    block_size = Int(floor(sqrt(n)))

    # Cleanup case artifacts and large objects.
    rm(outdir; recursive=true, force=true)
    model = nothing
    pheno = nothing
    global geno = nothing
    GC.gc()

    logline(@sprintf("CASE_DONE N=%d P=%d L=%d outer=%d block_size=%d t_mcmc=%.3f", n, p, l, outer, block_size, t_mcmc))
    return Dict(
        "N" => n,
        "P" => p,
        "L_input" => l,
        "outer" => outer,
        "block_size" => block_size,
        "t_x" => t_x,
        "t_geno" => t_geno,
        "t_mcmc" => t_mcmc,
    )
end

function fit_linear_two_points(x1, y1, x2, y2)
    b = (y2 - y1) / (x2 - x1)
    a = y1 - b * x1
    return a, b
end

cases = Any[]
for p in P_LIST
    for l in L_LIST
        push!(cases, run_case(N, p, l; seed=2026 + p + l))
    end
end

# Organize by P and infer per-outer decomposition from two chain lengths.
per_p = Dict{Int,Any}()
for p in P_LIST
    subset = [c for c in cases if c["P"] == p]
    sort!(subset, by = c -> c["outer"])

    if length(subset) < 2
        continue
    end

    c1, c2 = subset[1], subset[2]
    o1, o2 = c1["outer"], c2["outer"]
    t1, t2 = c1["t_mcmc"], c2["t_mcmc"]

    per_outer = (t2 - t1) / (o2 - o1)
    precompute_plus_overhead = t1 - per_outer * o1

    target_outer = Int(floor(TARGET_L / c1["block_size"]))
    target_t_mcmc = precompute_plus_overhead + per_outer * target_outer

    per_p[p] = Dict(
        "outer_points" => [o1, o2],
        "t_points" => [t1, t2],
        "per_outer" => per_outer,
        "precompute_plus_overhead" => precompute_plus_overhead,
        "target_outer" => target_outer,
        "target_t_mcmc" => target_t_mcmc,
    )

    logline(@sprintf("P=%d decomposition pre=%.3f per_outer=%.3f target_outer=%d target_t=%.3f", p, precompute_plus_overhead, per_outer, target_outer, target_t_mcmc))
end

# Extrapolate target_t_mcmc across P to TARGET_P using two largest available P values.
sorted_p = sort(collect(keys(per_p)))
est_seconds = NaN
est_hours = NaN
if length(sorted_p) >= 2
    p1, p2 = sorted_p[end-1], sorted_p[end]
    y1, y2 = per_p[p1]["target_t_mcmc"], per_p[p2]["target_t_mcmc"]
    a, b = fit_linear_two_points(p1, y1, p2, y2)
    est_seconds = a + b * TARGET_P
    est_hours = est_seconds / 3600

    logline(@sprintf("EXTRAPOLATED from P=%d,%d => target P=%d L=%d t_mcmc=%.3f sec (%.3f h)", p1, p2, TARGET_P, TARGET_L, est_seconds, est_hours))
end

open(OUT_TXT, "w") do io
    println(io, "timestamp=", now())
    println(io, "N=", N)
    println(io, "P_LIST=", join(P_LIST, ","))
    println(io, "L_LIST=", join(L_LIST, ","))
    println(io, "TARGET_P=", TARGET_P)
    println(io, "TARGET_L=", TARGET_L)
    println(io, "BLAS_THREADS=", BLAS.get_num_threads())
    for c in cases
        @printf(io, "CASE N=%d P=%d L=%d outer=%d block_size=%d t_x=%.6f t_geno=%.6f t_mcmc=%.6f\n",
                c["N"], c["P"], c["L_input"], c["outer"], c["block_size"], c["t_x"], c["t_geno"], c["t_mcmc"])
    end
    for p in sorted_p
        d = per_p[p]
        @printf(io, "DECOMP P=%d pre=%.6f per_outer=%.6f target_outer=%d target_t=%.6f\n",
                p, d["precompute_plus_overhead"], d["per_outer"], d["target_outer"], d["target_t_mcmc"])
    end
    if isfinite(est_seconds)
        @printf(io, "SUMMARY_EST target_seconds=%.6f target_hours=%.6f\n", est_seconds, est_hours)
    end
end

logline("RESULTS_WRITTEN $(OUT_TXT)")
if isfinite(est_seconds)
    println(@sprintf("SUMMARY_EST target_seconds=%.3f target_hours=%.3f", est_seconds, est_hours))
end
