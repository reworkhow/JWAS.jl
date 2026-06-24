# Multi-Trait Annotated BayesR Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build dense 2-trait annotated BayesR with a 16-state trait-by-magnitude prior, trait-specific annotation effects for magnitude classes, and full cross-trait marker-effect covariance.

**Architecture:** Reuse the existing annotated multi-trait BayesC active-pattern tree as the first annotation layer, then add trait-specific nested BayesR magnitude layers. Add a dense multi-trait BayesR sampler that samples the active pattern by aggregating 16 state weights, samples the within-pattern magnitude state, stores full latent `beta` effects for the `G` update, and updates realized marker effects through BayesR scale matrices.

**Tech Stack:** Julia, JWAS dense genotype sampler internals, Distributions.jl, LinearAlgebra, Test stdlib, existing `MarkerAnnotations` storage.

---

## Reference Design

Read first:

- `docs/plans/2026-06-24-multitrait-annotated-bayesr-design.md`
- `src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`
- `src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl`
- `src/1.JWAS/src/MCMC/annotation_updates.jl`
- `src/1.JWAS/src/markers/annotation_setup.jl`
- `src/1.JWAS/src/variance_components.jl`
- `test/unit/test_annotated_bayesr.jl`
- `test/unit/test_multitrait_mcmc.jl`

Use TDD. Do not implement a later task until that task's failing tests fail for the expected reason.

### Task 1: Add 16-State Helper Tests

**Files:**

- Modify: `test/unit/test_annotated_bayesr.jl`
- Modify later: `src/1.JWAS/src/markers/annotation_setup.jl`

**Step 1: Write failing tests**

Add tests under `@testset "Annotated BayesR API and validation"`:

```julia
@testset "multi-trait BayesR state helpers" begin
    states = JWAS.annotated_bayesr_mt_state_keys()

    @test length(states) == 16
    @test states[1] == [0, 0]
    @test states[end] == [3, 3]
    @test length(unique(Tuple.(states))) == 16

    for (idx, state) in enumerate(states)
        @test JWAS.annotated_bayesr_mt_state_index(state) == idx
    end

    @test JWAS.annotated_bayesr_mt_pattern([0, 0]) == [0.0, 0.0]
    @test JWAS.annotated_bayesr_mt_pattern([2, 0]) == [1.0, 0.0]
    @test JWAS.annotated_bayesr_mt_pattern([0, 3]) == [0.0, 1.0]
    @test JWAS.annotated_bayesr_mt_pattern([1, 2]) == [1.0, 1.0]
end
```

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'
```

Expected: FAIL with `UndefVarError: annotated_bayesr_mt_state_keys not defined`.

**Step 3: Implement state helpers**

In `src/1.JWAS/src/markers/annotation_setup.jl`, add near existing annotated BayesC state helpers:

```julia
const ANNOTATED_BAYESR_MT_STATES = let
    states = Vector{Vector{Int}}()
    push!(states, [0, 0])
    for k in 1:3
        push!(states, [k, 0])
    end
    for k in 1:3
        push!(states, [0, k])
    end
    for k in 1:3, l in 1:3
        push!(states, [k, l])
    end
    Tuple(Tuple(s) for s in states)
end

function annotated_bayesr_mt_state_keys()
    return [collect(state) for state in ANNOTATED_BAYESR_MT_STATES]
end

function annotated_bayesr_mt_state_index(state::AbstractVector{<:Integer})
    length(state) == 2 || error("Annotated multi-trait BayesR expects 2-trait class labels.")
    key = (Int(state[1]), Int(state[2]))
    for (idx, candidate) in enumerate(ANNOTATED_BAYESR_MT_STATES)
        key == candidate && return idx
    end
    error("Annotated multi-trait BayesR expects class labels in 0:3.")
end

function annotated_bayesr_mt_pattern(state::AbstractVector{<:Integer})
    length(state) == 2 || error("Annotated multi-trait BayesR expects 2-trait class labels.")
    return [state[1] > 0 ? 1.0 : 0.0, state[2] > 0 ? 1.0 : 0.0]
end
```

**Step 4: Run test to verify it passes**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'
```

Expected: existing tests plus new helper test PASS.

**Step 5: Commit**

```bash
git add src/1.JWAS/src/markers/annotation_setup.jl test/unit/test_annotated_bayesr.jl
git commit -m "test: cover multi-trait BayesR state helpers"
```

### Task 2: Initialize Dense 2-Trait Annotated BayesR State

**Files:**

- Modify: `test/unit/test_annotated_bayesr.jl`
- Modify: `src/1.JWAS/src/markers/readgenotypes.jl`
- Modify: `src/1.JWAS/src/markers/annotation_setup.jl`
- Modify: `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`

**Step 1: Write failing tests**

Replace the current "rejects unsupported BayesR annotation modes" expectation for 2-trait annotated BayesR with a setup acceptance test:

```julia
@testset "initializes dense 2-trait annotated BayesR state" begin
    annotations = rand(Float64, 5, 2)
    geno = get_genotypes(
        genofile,
        [1.0 0.25; 0.25 1.0];
        method="BayesR",
        annotations=annotations,
        separator=',',
        quality_control=false,
    )
    model = build_model(
        "y1 = intercept + geno\ny2 = intercept + geno",
        [1.0 0.2; 0.2 1.0],
    )
    ann = model.M[1].annotations

    @test model.M[1].ntraits == 2
    @test ann !== false
    @test ann.nsteps == 7
    @test ann.nclasses == 16
    @test size(ann.coefficients) == (size(annotations, 2) + 1, 7)
    @test size(ann.snp_pi) == (model.M[1].nMarkers, 16)
    @test all(abs.(sum(ann.snp_pi, dims=2) .- 1.0) .< 1e-10)
end
```

Keep rejection tests for streaming, RRM, non-dense storage, and trait counts other than two.

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'
```

Expected: FAIL because annotated multi-trait BayesR is currently rejected or initialized as 4-class single-trait BayesR.

**Step 3: Implement model-dependent initialization**

In `readgenotypes.jl`, change `build_marker_annotations` for `method == "BayesR"` so it keeps the existing single-trait path but allows multi-trait finalization later when needed. Use `annotation_start_pi` similarly to BayesC.

In `annotation_setup.jl`, add:

```julia
function initialize_bayesr_mt_annotations!(genotypei::Genotypes)
    if genotypei.annotations === false || genotypei.method != "BayesR" || genotypei.ntraits == 1
        return nothing
    end
    genotypei.ntraits == 2 || error("Annotated multi-trait BayesR currently supports exactly 2 traits.")

    design_matrix = genotypei.annotations.design_matrix
    coeffs = zeros(Float64, size(design_matrix, 2), 7)
    start_row = fill(1.0 / 16.0, 16)

    genotypei.annotations = MarkerAnnotations(
        design_matrix;
        nsteps=7,
        nclasses=16,
        coefficients=coeffs,
        snp_pi=repeat(reshape(start_row, 1, :), genotypei.nMarkers, 1),
    )
    genotypei.π = copy(start_row)
    return nothing
end
```

Then update:

```julia
function finalize_marker_annotation_setup!(genotypei::Genotypes)
    ...
    if genotypei.method == "BayesR"
        genotypei.ntraits == 1 ? return nothing : initialize_bayesr_mt_annotations!(genotypei)
    end
end
```

In `MCMC_BayesianAlphabet.jl`, relax the multi-trait annotation guard to allow:

```julia
Mi.method == "BayesR" && mme.nModels == 2 && Mi.storage_mode == :dense && Mi.G.constraint == false
```

**Step 4: Run test to verify it passes**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'
```

Expected: setup tests PASS; later end-to-end tests may still be absent.

**Step 5: Commit**

```bash
git add src/1.JWAS/src/markers/readgenotypes.jl src/1.JWAS/src/markers/annotation_setup.jl src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl test/unit/test_annotated_bayesr.jl
git commit -m "feat: initialize multi-trait annotated BayesR"
```

### Task 3: Build 16-State Annotation Prior Refresh

**Files:**

- Modify: `test/unit/test_annotated_bayesr.jl`
- Modify: `src/1.JWAS/src/MCMC/annotation_updates.jl`

**Step 1: Write failing tests**

Add deterministic prior rebuild tests:

```julia
@testset "multi-trait annotated BayesR prior rows" begin
    design = [1.0 0.0; 1.0 1.0]
    ann = JWAS.MarkerAnnotations(
        design;
        nsteps=7,
        nclasses=16,
        coefficients=zeros(Float64, 2, 7),
        snp_pi=zeros(Float64, 2, 16),
    )

    JWAS.rebuild_bayesr_mt_priors!(ann)

    @test all(abs.(sum(ann.snp_pi, dims=2) .- 1.0) .< 1e-10)
    @test ann.snp_pi[1, JWAS.annotated_bayesr_mt_state_index([0, 0])] ≈ 0.5

    # Make trait 1 large-class probability annotation-sensitive without changing trait 2.
    ann.coefficients[:, 5] .= [0.0, 3.0]
    JWAS.rebuild_bayesr_mt_priors!(ann)

    idx_30 = JWAS.annotated_bayesr_mt_state_index([3, 0])
    idx_03 = JWAS.annotated_bayesr_mt_state_index([0, 3])
    @test ann.snp_pi[2, idx_30] > ann.snp_pi[1, idx_30]
    @test ann.snp_pi[2, idx_03] ≈ ann.snp_pi[1, idx_03]
end
```

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'
```

Expected: FAIL with `UndefVarError: rebuild_bayesr_mt_priors! not defined`.

**Step 3: Implement prior rebuild**

In `annotation_updates.jl`, add:

```julia
function rebuild_bayesr_mt_priors!(ann)
    probs = clamp.(cdf.(Normal(), ann.mu), eps(Float64), 1 - eps(Float64))

    p1 = probs[:, 1]
    p2 = probs[:, 2]
    p3 = probs[:, 3]

    pattern00 = 1 .- p1
    pattern11 = p1 .* p2
    pattern10 = p1 .* (1 .- p2) .* p3
    pattern01 = p1 .* (1 .- p2) .* (1 .- p3)

    q11 = probs[:, 4]
    q12 = probs[:, 5]
    q21 = probs[:, 6]
    q22 = probs[:, 7]

    mag1 = hcat(1 .- q11, q11 .* (1 .- q12), q11 .* q12)
    mag2 = hcat(1 .- q21, q21 .* (1 .- q22), q21 .* q22)

    ann.snp_pi[:, :] .= 0.0
    ann.snp_pi[:, annotated_bayesr_mt_state_index([0, 0])] .= pattern00

    for k in 1:3
        ann.snp_pi[:, annotated_bayesr_mt_state_index([k, 0])] .= pattern10 .* mag1[:, k]
        ann.snp_pi[:, annotated_bayesr_mt_state_index([0, k])] .= pattern01 .* mag2[:, k]
    end
    for k in 1:3, l in 1:3
        ann.snp_pi[:, annotated_bayesr_mt_state_index([k, l])] .= pattern11 .* mag1[:, k] .* mag2[:, l]
    end

    row_sums = sum(ann.snp_pi, dims=2)
    ann.snp_pi[:, :] ./= row_sums
    return nothing
end
```

**Step 4: Run test to verify it passes**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'
```

Expected: prior rebuild tests PASS.

**Step 5: Commit**

```bash
git add src/1.JWAS/src/MCMC/annotation_updates.jl test/unit/test_annotated_bayesr.jl
git commit -m "feat: build multi-trait BayesR annotation priors"
```

### Task 4: Add Multi-Trait BayesR Annotation Response Updates

**Files:**

- Modify: `test/unit/test_annotated_bayesr.jl`
- Modify: `src/1.JWAS/src/MCMC/annotation_updates.jl`

**Step 1: Write failing tests**

Add:

```julia
@testset "multi-trait BayesR annotation responses" begin
    delta = [Int[0, 1, 2, 3, 0, 3], Int[0, 0, 0, 0, 2, 1]]
    responses, active_sets = JWAS.bayesr_mt_step_indicators(delta)

    @test responses[1] == Int[0, 1, 1, 1, 1, 1] # active
    @test responses[2] == Int[0, 0, 0, 0, 0, 1] # shared among all markers, active set filters
    @test responses[3] == Int[0, 1, 1, 1, 0, 0] # 10 among singleton markers
    @test responses[4] == Int[0, 0, 1, 1, 0, 1] # trait 1 >= medium
    @test responses[5] == Int[0, 0, 0, 1, 0, 1] # trait 1 large
    @test responses[6] == Int[0, 0, 0, 0, 1, 0] # trait 2 >= medium
    @test responses[7] == Int[0, 0, 0, 0, 0, 0] # trait 2 large

    @test active_sets[1] == collect(1:6)
    @test active_sets[2] == [2, 3, 4, 5, 6]
    @test active_sets[3] == [2, 3, 4, 5]
    @test active_sets[4] == [2, 3, 4, 6]
    @test active_sets[5] == [3, 4, 6]
    @test active_sets[6] == [5, 6]
    @test active_sets[7] == [5]
end
```

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'
```

Expected: FAIL with missing `bayesr_mt_step_indicators`.

**Step 3: Implement responses and update dispatch**

In `annotation_updates.jl`, add:

```julia
function bayesr_mt_step_indicators(deltaArray::AbstractVector)
    c1 = Int.(deltaArray[1])
    c2 = Int.(deltaArray[2])
    active1 = c1 .> 0
    active2 = c2 .> 0
    any_active = active1 .| active2
    shared = active1 .& active2
    singleton = xor.(active1, active2)

    z1 = Int.(any_active)
    z2 = Int.(shared)
    z3 = Int.(active1 .& singleton)
    z4 = Int.(c1 .>= 2)
    z5 = Int.(c1 .== 3)
    z6 = Int.(c2 .>= 2)
    z7 = Int.(c2 .== 3)

    active_sets = (
        collect(eachindex(c1)),
        findall(any_active),
        findall(singleton),
        findall(active1),
        findall(c1 .>= 2),
        findall(active2),
        findall(c2 .>= 2),
    )
    return (z1, z2, z3, z4, z5, z6, z7), active_sets
end
```

Update `update_marker_annotation_priors!`:

```julia
if Mi.method == "BayesR" && Mi.ntraits == 2
    responses, active_sets = bayesr_mt_step_indicators(Mi.δ)
    for step in 1:ann.nsteps
        sample_nested_annotation_probit_step!(ann, step, responses[step], active_sets[step])
    end
    rebuild_bayesr_mt_priors!(ann)
    Mi.π = vec(mean(ann.snp_pi, dims=1))
    return nothing
end
```

Keep the existing single-trait BayesR branch guarded by `Mi.ntraits == 1`.

**Step 4: Run test to verify it passes**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'
```

Expected: response tests PASS.

**Step 5: Commit**

```bash
git add src/1.JWAS/src/MCMC/annotation_updates.jl test/unit/test_annotated_bayesr.jl
git commit -m "feat: update multi-trait BayesR annotation steps"
```

### Task 5: Add Marker Variance Initialization For 16-State Priors

**Files:**

- Modify: `test/unit/test_annotated_bayesr.jl`
- Modify: `src/1.JWAS/src/markers/tools4genotypes.jl`

**Step 1: Write failing tests**

Add:

```julia
@testset "multi-trait BayesR genetic2marker uses class scales" begin
    geno = get_genotypes(
        genofile,
        [2.0 0.4; 0.4 1.5];
        method="BayesR",
        annotations=rand(Float64, 5, 2),
        separator=',',
        quality_control=false,
    )
    model = build_model("y1 = intercept + geno\ny2 = intercept + geno", [1.0 0.2; 0.2 1.0])
    Mi = model.M[1]
    Mi.annotations.snp_pi .= 0.0
    Mi.annotations.snp_pi[:, JWAS.annotated_bayesr_mt_state_index([3, 3])] .= 1.0

    JWAS.genetic2marker(Mi, Mi.annotations.snp_pi)

    @test Mi.G.val ≈ Mi.genetic_variance.val ./ Mi.sum2pq
end
```

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'
```

Expected: FAIL because matrix-valued `Pi` is supported only for annotated 2-trait BayesC.

**Step 3: Implement BayesR matrix `genetic2marker`**

Extend `genetic2marker(M::Genotypes, snp_pi::AbstractMatrix)`:

```julia
if M.method == "BayesR" && M.ntraits == 2
    size(snp_pi, 1) == M.nMarkers || error("Annotated multi-trait BayesR snp_pi row count must match nMarkers.")
    size(snp_pi, 2) == 16 || error("Annotated multi-trait BayesR snp_pi must have 16 columns.")

    allele_freq = vec(Float64.(M.alleleFreq))
    twopq = 2.0 .* allele_freq .* (1.0 .- allele_freq)
    states = annotated_bayesr_mt_state_keys()
    denom = zeros(Float64, 2, 2)

    for (state_idx, state) in enumerate(states)
        s1 = sqrt(BAYESR_GAMMA[state[1] + 1])
        s2 = sqrt(BAYESR_GAMMA[state[2] + 1])
        probs = Float64.(snp_pi[:, state_idx])
        denom[1, 1] += sum(twopq .* probs .* s1^2)
        denom[2, 2] += sum(twopq .* probs .* s2^2)
        denom[1, 2] += sum(twopq .* probs .* s1 .* s2)
    end
    denom[2, 1] = denom[1, 2]
    all(denom .> 0) || error("Annotated multi-trait BayesR implied covariance denominator must be positive.")
    M.G.val = M.genetic_variance.val ./ denom
    return nothing
end
```

Update `set_marker_hyperparameters_variances_and_pi` to route annotated 2-trait BayesR through this matrix method.

**Step 4: Run test to verify it passes**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'
```

Expected: variance initialization test PASS.

**Step 5: Commit**

```bash
git add src/1.JWAS/src/markers/tools4genotypes.jl test/unit/test_annotated_bayesr.jl
git commit -m "feat: initialize multi-trait BayesR marker covariance"
```

### Task 6: Implement Dense Multi-Trait BayesR Marker Sampler

**Files:**

- Create: `src/1.JWAS/src/markers/BayesianAlphabet/MTBayesR.jl`
- Modify: `src/1.JWAS/src/JWAS.jl`
- Modify: `test/unit/test_multitrait_mcmc.jl`

**Step 1: Write failing exact-probability test**

Add helper to `test/unit/test_multitrait_mcmc.jl`:

```julia
function exact_mt_bayesr_state_probs(x, y_by_trait, vare, G, priors, states)
    xp = dot(x, x)
    Rinv = inv(vare)
    Ginv = inv(G)
    w = [dot(x, y_trait) for y_trait in y_by_trait]
    log_delta = zeros(Float64, length(states))

    for (idx, state) in enumerate(states)
        S = Diagonal(sqrt.(JWAS.BAYESR_GAMMA[Int.(state) .+ 1]))
        C = Ginv + xp * Matrix(S) * Rinv * Matrix(S)
        b = Matrix(S) * Rinv * w
        ghat = C \ b
        log_delta[idx] = log(priors[1, idx]) - 0.5 * log(det(C)) + 0.5 * dot(b, ghat)
    end

    log_delta .-= maximum(log_delta)
    probs = exp.(log_delta)
    probs ./= sum(probs)
    return probs
end
```

Add an empirical sampler test:

```julia
@testset "Multi-trait BayesR dense sampler matches one-marker target" begin
    x = [1.0, -0.5, 0.25, 1.2]
    y_by_trait = [[0.5, -0.1, 0.2, 0.8], [-0.3, 0.4, 0.1, 0.7]]
    vare = [1.0 0.2; 0.2 1.2]
    G = [0.8 0.25; 0.25 0.9]
    states = JWAS.annotated_bayesr_mt_state_keys()
    priors = reshape(fill(1 / 16, 16), 1, 16)
    exact = exact_mt_bayesr_state_probs(x, y_by_trait, vare, G, priors, states)

    empirical = empirical_mt_bayesr_state_probs(x, y_by_trait, vare, G, priors, states)
    @test maximum(abs.(empirical .- exact)) < 0.03
end
```

Define `empirical_mt_bayesr_state_probs` analogously to existing BayesC empirical tests, calling the new sampler repeatedly and counting sampled states.

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: FAIL because `MTBayesR!` does not exist.

**Step 3: Implement sampler**

Create `MTBayesR.jl` with:

```julia
function MTBayesR!(genotypes, ycorr_array, vare)
    priors = genotypes.annotations === false ? genotypes.π : genotypes.annotations.snp_pi
    MTBayesR!(genotypes.mArray, genotypes.mRinvArray, genotypes.mpRinvm,
              ycorr_array, genotypes.β, genotypes.δ, genotypes.α,
              vare, genotypes.G.val, priors, BAYESR_GAMMA,
              annotated_bayesr_mt_state_keys())
end
```

The internal method should:

- compute all 16 log weights
- aggregate by `00/10/01/11` pattern using `logsumexp`
- sample pattern
- sample class state within pattern
- sample `β` from `N(mu, inv(C))`
- store full latent `β` for every marker, including `00`
- set realized marker effects with `α = S*β`
- update residual arrays by old minus new realized `α`
- store `δ[trait][marker] = class` in `0:3`

**Step 4: Include the file**

In `src/1.JWAS/src/JWAS.jl`, include `MTBayesR.jl` next to other Bayesian Alphabet files.

**Step 5: Run test to verify it passes**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: one-marker sampler test PASS.

**Step 6: Commit**

```bash
git add src/1.JWAS/src/markers/BayesianAlphabet/MTBayesR.jl src/1.JWAS/src/JWAS.jl test/unit/test_multitrait_mcmc.jl
git commit -m "feat: add dense multi-trait BayesR sampler"
```

### Task 7: Add `G` Update For Multi-Trait BayesR

**Files:**

- Modify: `test/unit/test_annotated_bayesr.jl`
- Modify: `src/1.JWAS/src/variance_components.jl`

**Step 1: Write failing tests**

Add unit tests for sufficient statistics:

```julia
@testset "multi-trait BayesR G sufficient statistics" begin
    beta = [Float64[0.2, 0.0, 0.4, 0.7], Float64[0.3, 0.5, 0.0, -0.2]]
    delta = [Int[3, 0, 3, 0], Int[1, 3, 0, 0]]
    ssq, nmarkers = JWAS.bayesr_mt_sigma_sufficient_statistics(beta, delta)

    @test nmarkers == 4
    @test size(ssq) == (2, 2)
    @test isapprox(ssq, beta_matrix_outer_sum(beta); atol=1e-12)
end
```

Use a local helper `beta_matrix_outer_sum` in the test to compute the expected
sum of `beta_j * beta_j'` over all markers. The fourth marker is `00`, so this
test protects the intended all-marker `G` update.

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'
```

Expected: FAIL because `bayesr_mt_sigma_sufficient_statistics` does not exist.

**Step 3: Implement sufficient statistics**

In `variance_components.jl`, add:

```julia
function bayesr_mt_sigma_sufficient_statistics(betaArray, deltaArray)
    ntraits = length(betaArray)
    nmarkers = length(betaArray[1])
    ssq = zeros(Float64, ntraits, ntraits)
    for j in 1:nmarkers
        beta_j = [Float64(betaArray[t][j]) for t in 1:ntraits]
        ssq .+= beta_j * beta_j'
    end
    return ssq, nmarkers
end
```

Update `sample_marker_effect_variance(Mi)`:

```julia
elseif Mi.method == "BayesR" && Mi.ntraits == 2
    ssq, nmarkers = bayesr_mt_sigma_sufficient_statistics(Mi.β, Mi.δ)
    Mi.G.val = rand(InverseWishart(Mi.G.df + nmarkers, convert(Array, Symmetric(Mi.G.scale + ssq))))
```

**Step 4: Run test to verify it passes**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'
```

Expected: sufficient-statistic test PASS.

**Step 5: Commit**

```bash
git add src/1.JWAS/src/variance_components.jl test/unit/test_annotated_bayesr.jl
git commit -m "feat: sample multi-trait BayesR covariance"
```

### Task 8: Dispatch Multi-Trait BayesR In The MCMC Loop

**Files:**

- Modify: `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`
- Modify: `test/unit/test_multitrait_mcmc.jl`

**Step 1: Write failing integration test**

Add a tiny run:

```julia
@testset "Multi-trait annotated BayesR dense run" begin
    annotations = rand(Float64, 5, 2)
    global annotated_mt_bayesr = get_genotypes(
        genofile,
        [1.0 0.3; 0.3 1.0];
        separator=',',
        method="BayesR",
        quality_control=false,
        annotations=annotations,
    )
    model = build_model("y1 = intercept + annotated_mt_bayesr\ny2 = intercept + annotated_mt_bayesr", R)
    outdir = "test_mt_annotated_bayesr"
    output = runMCMC(
        model,
        phenotypes;
        chain_length=20,
        burnin=5,
        output_samples_frequency=5,
        output_folder=outdir,
        seed=123,
        printout_model_info=false,
        outputEBV=false,
        output_heritability=false,
    )

    @test haskey(output, "annotation coefficients annotated_mt_bayesr")
    @test all(abs.(sum(model.M[1].annotations.snp_pi, dims=2) .- 1.0) .< 1e-8)

    isdir(outdir) && rm(outdir, recursive=true, force=true)
end
```

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: FAIL because MCMC dispatch still lacks multi-trait BayesR call.

**Step 3: Implement dispatch**

In `MCMC_BayesianAlphabet.jl`, locate marker sampling dispatch and add:

```julia
elseif Mi.method == "BayesR" && Mi.ntraits == 2
    MTBayesR!(Mi, wArray, mme.R.val)
```

Ensure `Mi.β` is initialized as unscaled base effects and `Mi.α` remains realized effects.

**Step 4: Run test to verify it passes**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: integration test PASS.

**Step 5: Commit**

```bash
git add src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl test/unit/test_multitrait_mcmc.jl
git commit -m "feat: dispatch multi-trait annotated BayesR"
```

### Task 9: Update Output For 7 Annotation Steps And BayesR Classes

**Files:**

- Modify: `src/1.JWAS/src/output.jl`
- Modify: `test/unit/test_annotated_bayesr.jl`
- Modify: `test/unit/test_multitrait_mcmc.jl`

**Step 1: Write failing tests**

In the integration test, assert step labels:

```julia
coef = output["annotation coefficients annotated_mt_bayesr"]
@test Set(coef[!, :Step]) == Set([
    "zero_vs_active",
    "11_vs_singleton",
    "10_vs_01",
    "trait1_medium_or_large_vs_small",
    "trait1_large_vs_medium",
    "trait2_medium_or_large_vs_small",
    "trait2_large_vs_medium",
])
```

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: FAIL because output labels are missing or generic.

**Step 3: Implement labels**

In `output.jl`, extend annotation step labeling:

```julia
if Mi.method == "BayesR" && Mi.ntraits == 2 && Mi.annotations.nsteps == 7
    return [
        "zero_vs_active",
        "11_vs_singleton",
        "10_vs_01",
        "trait1_medium_or_large_vs_small",
        "trait1_large_vs_medium",
        "trait2_medium_or_large_vs_small",
        "trait2_large_vs_medium",
    ]
end
```

Also ensure `pi_<geno>` output can summarize 16-state priors, or add a separate clear output table for state probabilities.

**Step 4: Run test to verify it passes**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: output tests PASS.

**Step 5: Commit**

```bash
git add src/1.JWAS/src/output.jl test/unit/test_annotated_bayesr.jl test/unit/test_multitrait_mcmc.jl
git commit -m "feat: label multi-trait BayesR annotation output"
```

### Task 10: Guardrails And Regression Tests

**Files:**

- Modify: `test/unit/test_annotated_bayesr.jl`
- Modify: `test/unit/test_annotated_bayesc.jl`
- Modify: `test/unit/test_multitrait_mcmc.jl`
- Modify implementation files as needed

**Step 1: Add guardrail tests**

Cover:

- annotated multi-trait BayesR rejects `storage=:stream`
- annotated multi-trait BayesR rejects `ntraits != 2`
- annotated multi-trait BayesR rejects `constraint=true`
- annotated multi-trait BayesR rejects RRM
- single-trait annotated BayesR still initializes with 3 steps and 4 classes
- annotated multi-trait BayesC still initializes with 3 steps and 4 states

**Step 2: Run focused tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: all focused tests PASS.

**Step 3: Commit**

```bash
git add test/unit/test_annotated_bayesr.jl test/unit/test_annotated_bayesc.jl test/unit/test_multitrait_mcmc.jl src/1.JWAS/src
git commit -m "test: cover annotated BayesR guardrails"
```

### Task 11: Documentation And Implementation Record

**Files:**

- Create: `docs/plans/2026-06-24-multitrait-annotated-bayesr-implementation.md`
- Modify: docs manual files only if this feature should be user-facing immediately

**Step 1: Write implementation record**

Create:

```markdown
# Multi-Trait Annotated BayesR Implementation

## Summary

Implemented dense 2-trait annotated BayesR with a 16-state state space,
trait-specific magnitude annotation models, and a full global marker-effect
covariance matrix.

## Files Changed

- ...

## Tests

- ...

## Known Limitations

- dense 2-trait only
- no RRM
- no streaming
- full latent `beta` covariance update may need benchmarking for G12 mixing
```

**Step 2: Run docs only if manual docs changed**

Run:

```bash
julia --project=docs --startup-file=no docs/make.jl
```

Expected: PASS if docs changed. Skip and state not run if only plan docs changed.

**Step 3: Commit**

```bash
git add docs/plans/2026-06-24-multitrait-annotated-bayesr-implementation.md docs/src
git commit -m "docs: record multi-trait annotated BayesR implementation"
```

### Task 12: Final Verification

**Files:**

- No edits unless verification finds a defect.

**Step 1: Run focused tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: all PASS.

**Step 2: Run full test suite if time allows**

Run:

```bash
julia --project=. --startup-file=no test/runtests.jl
```

Expected: PASS.

**Step 3: Check git status**

Run:

```bash
git status --short
```

Expected: only unrelated pre-existing benchmark artifacts remain untracked.

**Step 4: Commit any verification fixes**

If fixes were needed:

```bash
git add <changed files>
git commit -m "fix: stabilize multi-trait annotated BayesR"
```

## Execution Options

Plan complete and saved to `docs/plans/2026-06-24-multitrait-annotated-bayesr-plan.md`.

Two execution options:

1. **Subagent-Driven (this session)** - dispatch fresh subagent per task, review between tasks, fast iteration.
2. **Parallel Session (separate)** - open a new session with `superpowers:executing-plans`, batch execution with checkpoints.

Choose one before implementation starts.
