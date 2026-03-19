# BayesR Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a first BayesR implementation to JWAS for single-trait dense genotype analysis, with concise setup printouts, strict validation, and targeted tests.

**Architecture:** BayesR will live inside the existing Bayesian Alphabet path. The implementation adds a new dense sampler file and method-specific branches in validation, hyperparameter setup, MCMC dispatch, variance updates, and output handling while explicitly rejecting streaming, fast-blocks, annotations, RRM, and multi-trait configurations in v1.

**Tech Stack:** Julia, JWAS internal MCMC code, `Distributions.jl`, `Test`, `DataFrames.jl`, `CSV.jl`, shell commands (`rg`, `julia --project=.`), git.

---

**Assumption:** When `method="BayesR"` and `Pi=0.0`, use a fixed 4-class starting vector such as `Float64[0.95, 0.03, 0.015, 0.005]`. If product direction changes, update Task 1 and Task 3 consistently before implementation.

### Task 1: Register BayesR And Add Validation Tests

**Files:**
- Create: `test/unit/test_bayesr.jl`
- Modify: `test/runtests.jl`
- Modify: `src/1.JWAS/src/input_data_validation.jl`
- Modify: `src/1.JWAS/src/markers/readgenotypes.jl`

**Step 1: Write the failing BayesR validation test**

Create `test/unit/test_bayesr.jl` with an initial testset like:

```julia
using Test, JWAS, DataFrames, CSV, JWAS.Datasets

phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
genofile = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstring=["NA"])

@testset "BayesR validation" begin
    geno = get_genotypes(genofile, 1.0, separator=',',
                         method="BayesR",
                         Pi=Float64[0.95, 0.03, 0.015, 0.005],
                         estimatePi=true)
    @test geno.method == "BayesR"
    @test geno.π == Float64[0.95, 0.03, 0.015, 0.005]

    bad_len = get_genotypes(genofile, 1.0, separator=',',
                            method="BayesR",
                            Pi=Float64[0.95, 0.05, 0.0])
    model_bad_len = build_model("y1 = intercept + bad_len", 1.0)
    err = @test_throws Exception runMCMC(model_bad_len, phenotypes;
                                         chain_length=5,
                                         burnin=0,
                                         output_samples_frequency=1,
                                         output_folder="test_bayesr_bad_len",
                                         seed=123,
                                         printout_model_info=false,
                                         outputEBV=false,
                                         output_heritability=false)
    @test occursin("length 4", sprint(showerror, err))

    stream_err = @test_throws Exception get_genotypes(genofile, 1.0, separator=',',
                                                      method="BayesR",
                                                      storage=:stream)
    @test occursin("dense", sprint(showerror, stream_err)) ||
          occursin("stream", sprint(showerror, stream_err))
end
```

Add the new file to `test/runtests.jl`:

```julia
@testset "BayesR" begin
    include(joinpath(@__DIR__, "unit", "test_bayesr.jl"))
end
```

**Step 2: Run the test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected: FAIL because `BayesR` is not yet an accepted method and/or BayesR-specific validation is missing.

**Step 3: Write the minimal validation implementation**

Update `src/1.JWAS/src/input_data_validation.jl` to:

```julia
if !(Mi.method in ["BayesL","BayesC","BayesB","BayesA","BayesR","RR-BLUP","GBLUP"])
    error(Mi.method, " is not available in JWAS. Please read the documentation.")
end

if Mi.method == "BayesR"
    mme.nModels == 1 || error("BayesR v1 supports single-trait analysis only.")
    Mi.storage_mode == :dense || error("BayesR v1 supports storage=:dense only.")
    mme.MCMCinfo.fast_blocks == false || error("BayesR v1 does not support fast_blocks.")
    mme.MCMCinfo.RRM == false || error("BayesR v1 does not support random regression model (RRM).")
    Mi.annotations === false || error("BayesR v1 does not support annotations.")
    (Mi.π isa AbstractVector && length(Mi.π) == 4) || Mi.π == 0.0 ||
        error("BayesR Pi must be a length 4 vector or 0.0 for defaults.")
end
```

Update `src/1.JWAS/src/markers/readgenotypes.jl` so the method is documented and `storage=:stream` rejects BayesR early:

```julia
if storage == :stream && method == "BayesR"
    error("storage=:stream MVP does not support BayesR.")
end
```

Add BayesR `Pi` normalization helpers in `readgenotypes.jl` so:
- `Pi=0.0` is preserved as the BayesR default sentinel
- explicit vectors are converted to `Vector{Float64}`
- invalid vector length / sign / sum are rejected before sampling

**Step 4: Run the test to verify it passes**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected: PASS for the validation-only checks.

**Step 5: Commit**

```bash
git add test/unit/test_bayesr.jl test/runtests.jl src/1.JWAS/src/input_data_validation.jl src/1.JWAS/src/markers/readgenotypes.jl
git commit -m "test: add BayesR validation coverage"
```

### Task 2: Implement The Dense BayesR Sampler And MCMC Dispatch

**Files:**
- Modify: `test/unit/test_bayesr.jl`
- Create: `src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`
- Modify: `src/1.JWAS/src/JWAS.jl`
- Modify: `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`

**Step 1: Write the failing sampler tests**

Extend `test/unit/test_bayesr.jl` with:

```julia
@testset "BayesR dense sampler" begin
    x1 = Float64[0, 1, 2, 1]
    x2 = Float64[2, 1, 0, 1]
    xArray = [x1, x2]
    xRinvArray = [x1, x2]
    xpRinvx = Float64[dot(x1, x1), dot(x2, x2)]
    yCorr = Float64[0.8, -0.1, 0.3, 0.5]
    α = zeros(Float64, 2)
    δ = ones(Int, 2)
    π = Float64[0.95, 0.03, 0.015, 0.005]
    gamma = Float64[0.0, 0.01, 0.1, 1.0]

    JWAS.BayesR!(xArray, xRinvArray, xpRinvx, yCorr, α, δ, 1.0, 0.2, π, gamma)

    @test all(1 .<= δ .<= 4)
    @test length(α) == 2
end

@testset "BayesR single-trait run" begin
    geno = get_genotypes(genofile, 1.0, separator=',',
                         method="BayesR",
                         Pi=Float64[0.95, 0.03, 0.015, 0.005],
                         estimatePi=true,
                         estimate_variance=true)
    model = build_model("y1 = intercept + geno", 1.0)

    output = runMCMC(model, phenotypes;
                     chain_length=30,
                     burnin=5,
                     output_samples_frequency=10,
                     output_folder="test_bayesr_dense",
                     seed=123,
                     printout_model_info=false,
                     outputEBV=false,
                     output_heritability=false,
                     fast_blocks=false)

    @test haskey(output, "marker effects geno")
    @test haskey(output, "location parameters")
    rm("test_bayesr_dense", recursive=true)
end
```

**Step 2: Run the tests to verify they fail**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected: FAIL with `UndefVarError`/method errors because `BayesR!` and the MCMC dispatch branch do not exist yet.

**Step 3: Write the minimal implementation**

Create `src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl` with a BayesABC-shaped dense sampler:

```julia
const BAYESR_GAMMA = Float64[0.0, 0.01, 0.1, 1.0]

function BayesR!(genotypes, ycorr, vare)
    BayesR!(genotypes.mArray, genotypes.mRinvArray, genotypes.mpRinvm,
            ycorr, genotypes.α[1], genotypes.δ[1], vare, genotypes.G.val, genotypes.π, BAYESR_GAMMA)
end

function BayesR!(xArray, xRinvArray, xpRinvx, yCorr, α, δ, vare, sigmaSq, π, gamma)
    # dense single-trait Gibbs sweep
end
```

Inside the kernel, implement:
- explicit zero-class handling
- stable log-sum-exp class probability normalization
- class label updates `δ[j] ∈ 1:4`
- class-specific normal sampling for nonzero classes
- in-place `yCorr` update with `BLAS.axpy!`

Wire it into `src/1.JWAS/src/JWAS.jl`:

```julia
include("markers/BayesianAlphabet/BayesR.jl")
```

Add a single-trait BayesR branch in `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`:

```julia
elseif Mi.method == "BayesR"
    BayesR!(Mi, ycorr, mme.R.val)
```

Reject any multi-trait BayesR path in the same file with a direct error rather than falling through.

**Step 4: Run the tests to verify they pass**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected: PASS for the sampler smoke tests; later output-specific assertions may still fail.

**Step 5: Commit**

```bash
git add test/unit/test_bayesr.jl src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl src/1.JWAS/src/JWAS.jl src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl
git commit -m "feat: add dense single-trait BayesR sampler"
```

### Task 3: Add BayesR Hyperparameter Setup, Initialization, And Concise Printouts

**Files:**
- Modify: `test/unit/test_bayesr.jl`
- Modify: `src/1.JWAS/src/markers/tools4genotypes.jl`
- Modify: `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`
- Modify: `src/1.JWAS/src/variance_components.jl`

**Step 1: Write the failing initialization tests**

Extend `test/unit/test_bayesr.jl` with:

```julia
@testset "BayesR initialization and priors" begin
    geno = get_genotypes(genofile, 1.0, separator=',',
                         method="BayesR",
                         Pi=0.0,
                         estimatePi=true,
                         estimate_variance=true)
    model = build_model("y1 = intercept + geno", 1.0)

    io = IOBuffer()
    redirect_stdout(io) do
        runMCMC(model, phenotypes;
                chain_length=10,
                burnin=0,
                output_samples_frequency=5,
                output_folder="test_bayesr_init",
                seed=321,
                printout_model_info=true,
                outputEBV=false,
                output_heritability=false,
                fast_blocks=false)
    end
    printed = String(take!(io))

    @test occursin("BayesR", printed)
    @test occursin("starting pi", lowercase(printed))
    @test occursin("gamma", lowercase(printed)) || occursin("mixture", lowercase(printed))
    rm("test_bayesr_init", recursive=true)
end
```

**Step 2: Run the test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected: FAIL because BayesR default `Pi`, `sigmaSq` derivation, class initialization, or printouts are incomplete.

**Step 3: Write the minimal implementation**

Update `src/1.JWAS/src/markers/tools4genotypes.jl` so `set_marker_hyperparameters_variances_and_pi` handles BayesR:

```julia
if mme.nModels == 1 && Mi.method == "BayesR"
    if Mi.π == 0.0
        Mi.π = Float64[0.95, 0.03, 0.015, 0.005]
    end
    denom = Mi.sum2pq * sum(BAYESR_GAMMA .* Mi.π)
    denom > 0 || error("BayesR implied variance denominator must be positive.")
    Mi.G.val = Mi.genetic_variance.val / denom
    Mi.G.scale = Mi.G.val * (Mi.G.df - 2) / Mi.G.df
end
```

Update `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl` to initialize BayesR state:

```julia
if Mi.method == "BayesR"
    Mi.δ = [ones(Int, Mi.nMarkers)]
    Mi.meanDelta = [zeros(Float64, Mi.nMarkers)]
    initialize_bayesr_indicators!(Mi)
    printstyled("BayesR gamma: [0.0, 0.01, 0.1, 1.0]\n", bold=false, color=:green)
    printstyled("BayesR starting pi: $(Mi.π)\n", bold=false, color=:green)
end
```

Implement `initialize_bayesr_indicators!(Mi)` in `MCMC_BayesianAlphabet.jl` or `BayesR.jl` so it:
- samples class labels from `Mi.π`
- prevents the degenerate all-zero start
- reports concise initial class counts

Update `src/1.JWAS/src/variance_components.jl` so `sample_marker_effect_variance(Mi)` handles BayesR by transforming current nonzero effects to the shared `sigmaSq` scale before calling the existing scalar variance sampler.

**Step 4: Run the test to verify it passes**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected: PASS for BayesR initialization, prior derivation, and concise setup printout checks.

**Step 5: Commit**

```bash
git add test/unit/test_bayesr.jl src/1.JWAS/src/markers/tools4genotypes.jl src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl src/1.JWAS/src/variance_components.jl
git commit -m "feat: add BayesR priors and initialization"
```

### Task 4: Add BayesR Pi Updates And Output Handling

**Files:**
- Modify: `test/unit/test_bayesr.jl`
- Modify: `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`
- Modify: `src/1.JWAS/src/markers/Pi.jl`
- Modify: `src/1.JWAS/src/output.jl`

**Step 1: Write the failing output tests**

Extend `test/unit/test_bayesr.jl` with:

```julia
@testset "BayesR output" begin
    geno = get_genotypes(genofile, 1.0, separator=',',
                         method="BayesR",
                         Pi=Float64[0.95, 0.03, 0.015, 0.005],
                         estimatePi=true,
                         estimate_variance=true)
    model = build_model("y1 = intercept + geno", 1.0)

    output = runMCMC(model, phenotypes;
                     chain_length=40,
                     burnin=10,
                     output_samples_frequency=10,
                     output_folder="test_bayesr_output",
                     seed=1234,
                     printout_model_info=false,
                     outputEBV=false,
                     output_heritability=false,
                     fast_blocks=false)

    @test haskey(output, "pi_geno")
    @test nrow(output["pi_geno"]) == 4
    @test Set(output["pi_geno"][!, :π]) == Set(["Pi1", "Pi2", "Pi3", "Pi4"])
    @test all(0.0 .<= output["marker effects geno"][!, :Model_Frequency] .<= 1.0)

    pi_sample_lines = filter(!isempty, readlines(joinpath("test_bayesr_output", "MCMC_samples_pi_geno.txt")))
    @test length(pi_sample_lines) == 4

    rm("test_bayesr_output", recursive=true)
end
```

**Step 2: Run the test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected: FAIL because BayesR `pi` is still collapsed or BayesR-specific posterior bookkeeping is missing.

**Step 3: Write the minimal implementation**

Update `src/1.JWAS/src/markers/Pi.jl` with a BayesR mixture update helper:

```julia
function samplePiBayesR(delta::AbstractVector{<:Integer}, nclasses::Integer)
    counts = [count(==(k), delta) + 1 for k in 1:nclasses]
    return vec(rand(Dirichlet(counts)))
end
```

Update `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`:

```julia
if Mi.method == "BayesR" && Mi.estimatePi == true
    Mi.π = samplePiBayesR(Mi.δ[1], 4)
end
```

Update `src/1.JWAS/src/output.jl` so BayesR preserves all four `pi` entries:

```julia
function collapse_pi_for_output(Mi, pi_value)
    if Mi.method == "BayesR"
        return OrderedDict("Pi1" => pi_value[1], "Pi2" => pi_value[2], "Pi3" => pi_value[3], "Pi4" => pi_value[4])
    end
    if Mi.annotations === false && Mi.ntraits == 1 && pi_value isa AbstractVector
        return pi_value[1]
    end
    return pi_value
end
```

Update posterior accumulation so BayesR `meanDelta` tracks nonzero frequency:

```julia
if Mi.method == "BayesR"
    Mi.meanDelta[trait] += ((Mi.δ[trait] .> 1) .- Mi.meanDelta[trait]) / nsamples
else
    Mi.meanDelta[trait] += (Mi.δ[trait] - Mi.meanDelta[trait]) / nsamples
end
```

**Step 4: Run the test to verify it passes**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected: PASS with 4-row `pi` output and BayesR marker frequencies in `[0, 1]`.

**Step 5: Commit**

```bash
git add test/unit/test_bayesr.jl src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl src/1.JWAS/src/markers/Pi.jl src/1.JWAS/src/output.jl
git commit -m "feat: add BayesR output handling"
```

### Task 5: Final Verification And Regression Pass

**Files:**
- Verify only: `test/unit/test_bayesr.jl`
- Verify only: `test/runtests.jl`
- Verify only: BayesR implementation files modified in Tasks 1-4

**Step 1: Run the focused BayesR tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected: PASS.

**Step 2: Run the full JWAS test suite**

Run:

```bash
julia --project=. --startup-file=no test/runtests.jl
```

Expected: PASS for the full suite, including the new BayesR tests.

**Step 3: Inspect git diff before final commit**

Run:

```bash
git status --short
git diff --stat
```

Expected: only BayesR-related files remain modified.

**Step 4: Commit the final verification state**

```bash
git add test/unit/test_bayesr.jl test/runtests.jl src/1.JWAS/src/input_data_validation.jl src/1.JWAS/src/markers/readgenotypes.jl src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl src/1.JWAS/src/JWAS.jl src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl src/1.JWAS/src/markers/tools4genotypes.jl src/1.JWAS/src/variance_components.jl src/1.JWAS/src/markers/Pi.jl src/1.JWAS/src/output.jl
git commit -m "feat: add single-trait dense BayesR"
```
