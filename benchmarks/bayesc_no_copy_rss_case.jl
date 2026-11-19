# One fresh-process, production-path case for run_bayesc_no_copy_rss.py.
# Arguments: output-directory n p source quality-control seed outer burnin block-size.
using JWAS, DataFrames, CSV, Random, LinearAlgebra, Statistics

BLAS.set_num_threads(1)
length(ARGS) == 9 || error("Expected output-directory n p adopted|copied|file true|false seed outer burnin block-size")
case_dir = abspath(ARGS[1])
n, p = parse.(Int, ARGS[2:3])
source = ARGS[4]
qc = parse(Bool, ARGS[5])
seed, outer, burnin, block_size = parse.(Int, ARGS[6:9])
source in ("adopted", "copied", "file") || error("Unsupported source")
0 < burnin < outer || error("Require 0 < burnin < outer")
n >= 8 && p > block_size >= 1 || error("Require n >= 8 and p > block-size >= 1")
mkpath(case_dir)
source_path = realpath(pathof(JWAS))
println("JWAS source: ", source_path, "; Julia: ", VERSION)
println("Case: n=", n, " p=", p, " source=", source, " qc=", qc, " seed=", seed)
flush(stdout)

# Sys.maxrss documents bytes on every supported platform. Each snapshot is a
# process-lifetime high-water mark, including Julia, packages, JIT and fixtures.
rss_after_packages = Sys.maxrss()
rng = MersenneTwister(20261118)
X = rand(rng, Float32, n, p)
X .= floor.(3f0 .* X)
# Deterministically keep every SNP above the MAF threshold and polymorphic.
X[1, :] .= 0f0
X[2, :] .= 1f0
X[3, :] .= 2f0
ids = ["animal_$(lpad(i, 6, '0'))" for i in 1:n]
markers = ["rs_$(1000 + 7i)" for i in 1:p]
y = 0.3 .* Float64.(X[:, 1]) .- 0.2 .* Float64.(X[:, 2]) .+ randn(rng, n)
t = fld(3n, 4)
# Training-only input exercises JWAS's actual phenotype/genotype alignment.
phenotypes = DataFrame(ID=ids[1:t], y=y[1:t])

function write_fixture(path, matrix, ids, markers)
    frame = DataFrame([view(matrix, :, j) for j in axes(matrix, 2)], Symbol.(markers); copycols=false)
    insertcols!(frame, 1, :ID => ids; copycols=false)
    CSV.write(path, frame)
end

input = X
if source == "file"
    input = joinpath(case_dir, "fixture.csv")
    write_fixture(input, X, ids, markers)
    X = nothing
    GC.gc()
end
rss_before_import = Sys.maxrss()
started = time_ns()
global rss_genotypes = get_genotypes(input, 1.0;
    method="BayesC", Pi=0.65, quality_control=qc,
    copy_genotypes=(source == "copied"), double_precision=false)
import_seconds = (time_ns() - started) / 1e9
g = rss_genotypes
if source != "file"
    g.obsID .= ids
    g.markerID .= markers
end
@assert g.nMarkers == p "This benchmark requires all SNPs to pass QC"
@assert Set(g.markerID) == Set(markers)
@assert length(unique(g.markerID)) == p
import_aliases_caller = X !== nothing && g.genotypes === X
source == "adopted" && @assert import_aliases_caller "Adoption/all-pass QC must retain the input buffer"
import_storage = JWAS.genotype_storage_report(g)
rss_after_import = Sys.maxrss()

model = build_model("y = intercept + rss_genotypes", 1.0)
outputEBV(model, ids)
rss_before_mcmc = Sys.maxrss()
started = time_ns()
output = runMCMC(model, phenotypes;
    chain_length=outer * block_size, burnin=burnin,
    output_samples_frequency=10, seed=seed,
    fast_blocks=block_size, independent_blocks=false,
    shared_genotypes=(source != "copied"), double_precision=false,
    output_folder=joinpath(case_dir, "mcmc"), output_heritability=false,
    output_marker_effect_samples=false, printout_model_info=false,
    printout_frequency=outer + 1)
mcmc_seconds = (time_ns() - started) / 1e9
rss_after_mcmc = Sys.maxrss()
storage = JWAS.genotype_storage_report(g)
model_aliases_caller = X !== nothing &&
    (Base.mightalias(g.genotypes, X) || Base.mightalias(g.output_genotypes, X))
caller_bytes = X === nothing ? 0 : sizeof(X)
model_and_caller_bytes = storage.backing_matrix_bytes + (model_aliases_caller ? 0 : caller_bytes)

animals = CSV.read(joinpath(case_dir, "mcmc", "EBV_y.txt"), DataFrame; types=Dict(:ID => String))
animals.Trait = fill("y", nrow(animals))
marker = output["marker effects rss_genotypes"]
@assert nrow(animals) == n && length(unique(animals.ID)) == n
@assert Set(animals.ID) == Set(ids)
@assert nrow(marker) == p && length(unique(string.(marker.Marker_ID))) == p
@assert Set(string.(marker.Marker_ID)) == Set(markers)
@assert all(isfinite, animals.EBV) && all(isfinite, animals.PEV)
@assert all(isfinite, marker.Estimate) && all(isfinite, marker.SD)
@assert model.MCMCinfo.chain_length == outer
@assert g.nObs == t
CSV.write(joinpath(case_dir, "animals.csv"), animals)
CSV.write(joinpath(case_dir, "markers.csv"), marker)
metrics = (; n, p, training_rows=t, source, qc, seed, outer, burnin, block_size,
    saved_draws=fld(outer - burnin, 10), input_payload_bytes=4n*p,
    import_aliases_caller, import_backing_count=import_storage.backing_matrix_count,
    import_backing_bytes=import_storage.backing_matrix_bytes,
    model_backing_count=storage.backing_matrix_count,
    model_backing_bytes=storage.backing_matrix_bytes,
    caller_matrix_bytes=caller_bytes, model_aliases_caller,
    model_and_caller_unique_bytes=model_and_caller_bytes,
    block_gram_bytes=storage.block_gram_bytes,
    all_blocks_share_training=storage.all_blocks_share_training,
    output_shares_training_parent=storage.output_shares_training_parent,
    rss_after_packages, rss_before_import, rss_after_import, rss_before_mcmc, rss_after_mcmc,
    import_seconds, mcmc_seconds, julia_version=string(VERSION), source_path)
CSV.write(joinpath(case_dir, "metrics.csv"), DataFrame([metrics]))
println("Measured metrics: ", metrics)
flush(stdout)
# Keep compact exports; the launcher records logs, and raw draw files add no
# evidence to these keyed comparisons.
rm(joinpath(case_dir, "mcmc"); recursive=true)
source == "file" && rm(input)
