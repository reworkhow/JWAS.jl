# Independent-Blocks Fast-Block Parallelism Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add `independent_blocks=true` as an explicit approximate fast-block mode for single-trait and multi-trait BayesianAlphabet samplers, enabling block-parallel computation when users assume supplied marker blocks are independent.

**Architecture:** Keep the current exact sequential fast-block sweep as the default. Add an independent-block path that freezes the sweep-level corrected phenotype / residual, updates blocks independently, and reconciles block deltas at the end of the sweep. Implement the independent-block path sequentially first, then add thread-parallel execution.

**Tech Stack:** Julia, existing JWAS fast-block matrices, BayesianAlphabet samplers, `Threads.@threads`, Documenter, JWAS unit tests and production benchmarks.

---

### Task 1: Add `independent_blocks` to the public MCMC API

**Files:**
- Modify: `src/1.JWAS/src/JWAS.jl`
- Modify: `src/1.JWAS/src/types.jl`
- Test: `test/unit/test_misc_coverage.jl`

**Step 1: Write failing API tests**

Add tests that expect:

- `runMCMC(...; independent_blocks=false)` to be accepted
- `runMCMC(...; independent_blocks=true, fast_blocks=true)` to be accepted
- `runMCMC(...; independent_blocks=true, fast_blocks=false)` to error clearly

**Step 2: Run the targeted test**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: failure because `independent_blocks` is not a recognized keyword yet.

**Step 3: Add API plumbing**

Implement:

- `independent_blocks::Bool=false` in both `runMCMC` signatures
- an `independent_blocks` field on `MCMCinfo`
- validation that `independent_blocks=true` requires `fast_blocks != false`

**Step 4: Pass `independent_blocks` through `MCMCinfo`**

Update the `MCMCinfo(...)` constructor call in `src/1.JWAS/src/JWAS.jl`.

**Step 5: Run the targeted test again**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: API tests pass.

**Step 6: Commit**

```bash
git add src/1.JWAS/src/JWAS.jl src/1.JWAS/src/types.jl test/unit/test_misc_coverage.jl
git commit -m "feat: add independent_blocks MCMC option"
```

### Task 2: Add explicit user-provided block starts

**Files:**
- Modify: `src/1.JWAS/src/JWAS.jl`
- Modify: `src/1.JWAS/src/input_data_validation.jl`
- Test: `test/unit/test_misc_coverage.jl`

**Step 1: Write failing block-partition tests**

Add tests that expect:

- `fast_blocks=[1, 101, 201]` to be accepted
- unsorted starts to error
- duplicated starts to error
- starts not beginning with `1` to error
- starts beyond `nMarkers` to error

**Step 2: Run the targeted test**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: failure because `fast_blocks` currently does not accept user-supplied block starts.

**Step 3: Implement block-start parsing**

Extend current `fast_blocks` parsing:

- `false`: unchanged
- `true`: unchanged heuristic block size
- `Number`: unchanged fixed block size
- `AbstractVector{<:Integer}`: validate and store as block starts

**Step 4: Define chain-length behavior for explicit block starts**

For `fast_blocks::AbstractVector{<:Integer}`:

- do not divide `chain_length` by a block size
- interpret one MCMC iteration as one full sweep over supplied blocks
- use one within-block pass per sweep unless a later option adds block repeats

Keep legacy chain-length scaling for `fast_blocks=true` and `fast_blocks=<Number>`.

**Step 5: Run the targeted test again**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: block partition tests pass.

**Step 6: Commit**

```bash
git add src/1.JWAS/src/JWAS.jl src/1.JWAS/src/input_data_validation.jl test/unit/test_misc_coverage.jl
git commit -m "feat: accept explicit fast block starts"
```

### Task 3: Thread `independent_blocks` through block sampler dispatch

**Files:**
- Modify: `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/BayesABC.jl`
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl`
- Test: `test/unit/test_misc_coverage.jl`
- Test: `test/unit/test_multitrait_mcmc.jl`

**Step 1: Write failing dispatch tests**

Add tests that verify:

- single-trait BayesABC block dispatch receives `independent_blocks`
- single-trait BayesR block dispatch receives `independent_blocks`
- multi-trait BayesC block dispatch receives `independent_blocks`

Use small smoke tests that exercise the production `runMCMC` path.

**Step 2: Run targeted tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl"); include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: failure because block samplers do not accept the new argument.

**Step 3: Add keyword or positional arguments to block sampler wrappers**

Update wrappers so they accept the mode:

- `BayesABC_block!(...; independent_blocks=false)` or equivalent positional argument
- `BayesR_block!(...; independent_blocks=false)` or equivalent positional argument
- `MTBayesABC_block!(...; independent_blocks=false)` or equivalent positional argument

Update `MCMC_BayesianAlphabet.jl` call sites to pass `mme.MCMCinfo.independent_blocks`.

**Step 4: Preserve current behavior for `independent_blocks=false`**

The exact sequential block path must remain unchanged.

**Step 5: Run targeted tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl"); include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: dispatch tests pass.

**Step 6: Commit**

```bash
git add src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl src/1.JWAS/src/markers/BayesianAlphabet/BayesABC.jl src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl test/unit/test_misc_coverage.jl test/unit/test_multitrait_mcmc.jl
git commit -m "feat: dispatch independent block mode"
```

### Task 4: Implement independent blocks for single-trait BayesABC

**Files:**
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/BayesABC.jl`
- Test: `test/unit/test_misc_coverage.jl`

**Step 1: Write failing BayesABC independent-block tests**

Add a small production-path test covering:

- single-trait BayesC with `fast_blocks=[...]` and `independent_blocks=true`
- single-trait annotated BayesC with the same mode

**Step 2: Run the targeted test**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: failure because the independent block path is not implemented.

**Step 3: Extract exact sequential block helper**

Move the current block body into:

- `_BayesABC_block_sequential!(...)`

This helper is the existing exact behavior.

**Step 4: Add independent block helper**

Implement:

- `_BayesABC_block_independent!(...)`

Behavior:

- copy/freeze `yCorr` at the start of the marker sweep
- for each block, compute local `XpRinvycorr` from the frozen `yCorr`
- update only the block’s markers
- store `alpha_old - alpha_new` for the block
- after all blocks finish, apply all block deltas to the global `yCorr`

Start sequentially. Do not add threading in this task.

**Step 5: Run targeted tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: BayesABC independent-block tests pass.

**Step 6: Commit**

```bash
git add src/1.JWAS/src/markers/BayesianAlphabet/BayesABC.jl test/unit/test_misc_coverage.jl
git commit -m "feat: add independent blocks for BayesABC"
```

### Task 5: Implement independent blocks for single-trait BayesR

**Files:**
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`
- Test: `test/unit/test_misc_coverage.jl`
- Test: `test/unit/test_annotated_bayesr.jl`

**Step 1: Write failing BayesR independent-block tests**

Add tests covering:

- plain BayesR with `fast_blocks=[...]` and `independent_blocks=true`
- annotated BayesR with the same mode

**Step 2: Run targeted tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl"); include("test/unit/test_annotated_bayesr.jl")'
```

Expected: failure because BayesR independent blocks are not implemented.

**Step 3: Extract exact sequential block helper**

Move current block code into:

- `_BayesR_block_sequential!(...)`

**Step 4: Add independent block helper**

Implement:

- `_BayesR_block_independent!(...)`

Use the same sweep-level pattern as BayesABC:

- freeze `yCorr`
- update blocks independently
- collect deltas
- apply global correction once at the end

**Step 5: Run targeted tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl"); include("test/unit/test_annotated_bayesr.jl")'
```

Expected: BayesR tests pass.

**Step 6: Commit**

```bash
git add src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl test/unit/test_misc_coverage.jl test/unit/test_annotated_bayesr.jl
git commit -m "feat: add independent blocks for BayesR"
```

### Task 6: Implement independent blocks for multi-trait BayesC

**Files:**
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl`
- Test: `test/unit/test_multitrait_mcmc.jl`
- Test: `test/unit/test_annotated_bayesc.jl`

**Step 1: Write failing MT BayesC independent-block tests**

Add production-path tests covering:

- plain multi-trait BayesC sampler I with `fast_blocks=[...]`, `independent_blocks=true`
- plain multi-trait BayesC sampler II with the same mode
- annotated 2-trait BayesC sampler I with the same mode
- annotated 2-trait BayesC sampler II with the same mode

**Step 2: Run targeted tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl"); include("test/unit/test_annotated_bayesc.jl")'
```

Expected: failure because MT BayesC independent blocks are not implemented.

**Step 3: Extract exact sequential helpers**

Keep current exact behavior in helpers such as:

- `_MTBayesABC_block_samplerI_sequential!(...)`
- `_MTBayesABC_block_samplerII_sequential!(...)`

**Step 4: Add independent helpers**

Implement:

- `_MTBayesABC_block_samplerI_independent!(...)`
- `_MTBayesABC_block_samplerII_independent!(...)`

For each helper:

- freeze trait-wise `wArray` vectors at the start of the marker sweep
- compute each block’s local RHS from those frozen vectors
- update only that block’s markers
- collect per-trait block deltas
- apply all block deltas to the global `wArray` after all blocks finish

Reuse:

- `GlobalPiPrior`
- `MarkerSpecificPiPrior`
- `mt_bayesc_sampler_mode(...)`
- existing sampler I / sampler II state logic

**Step 5: Run targeted tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl"); include("test/unit/test_annotated_bayesc.jl")'
```

Expected: MT BayesC independent-block tests pass.

**Step 6: Commit**

```bash
git add src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl test/unit/test_multitrait_mcmc.jl test/unit/test_annotated_bayesc.jl
git commit -m "feat: add independent blocks for multitrait BayesC"
```

### Task 7: Add thread-parallel independent-block execution

**Files:**
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/BayesABC.jl`
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl`
- Test: `test/unit/test_misc_coverage.jl`
- Test: `test/unit/test_multitrait_mcmc.jl`

**Step 1: Make RNG handling explicit**

Before threading, refactor independent block helpers so block updates can use block-local RNG streams.

Do not rely on global `rand()` or `randn()` inside threaded block loops.

**Step 2: Add tests for deterministic execution shape**

Use block-diagonal toy data and controlled RNG setup to compare:

- sequential independent-block helper
- threaded independent-block helper

Only require strict equality if RNG order is explicitly controlled. Otherwise require valid execution and use exactness tests for statistical agreement.

**Step 3: Add threaded execution**

Use `Threads.@threads` over blocks for `independent_blocks=true`.

Implementation requirements:

- thread-local buffers
- block-local RNGs
- no shared writes during block updates
- one reduction/barrier step that applies block deltas after all blocks finish

**Step 4: Run targeted tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl"); include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: threaded independent-block tests pass.

**Step 5: Commit**

```bash
git add src/1.JWAS/src/markers/BayesianAlphabet/BayesABC.jl src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl test/unit/test_misc_coverage.jl test/unit/test_multitrait_mcmc.jl
git commit -m "feat: parallelize independent block sweeps"
```

### Task 8: Add exactness and approximation tests

**Files:**
- Modify: `test/unit/test_misc_coverage.jl`
- Modify: `test/unit/test_multitrait_mcmc.jl`

**Step 1: Add block-diagonal synthetic tests**

Construct synthetic data with:

```text
X_b' W X_c = 0  for all b != c
```

Cover:

- single-trait BayesC
- single-trait BayesR
- multi-trait BayesC

Expected:

- exact sequential fast blocks and independent blocks agree up to Monte Carlo tolerance

**Step 2: Add coupled-block smoke tests**

Use data where cross-block leakage is nonzero.

Expected:

- independent blocks run successfully
- tests do not assert equality with exact sequential fast blocks
- comments state that this is the approximate regime

**Step 3: Run targeted tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl"); include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: tests pass.

**Step 4: Commit**

```bash
git add test/unit/test_misc_coverage.jl test/unit/test_multitrait_mcmc.jl
git commit -m "test: cover independent block assumptions"
```

### Task 9: Benchmark production behavior and server scaling

**Files:**
- Create: `benchmarks/reports/2026-04-15-independent-blocks-fast-block-parallelism-report.md`
- Create: `docs/plans/2026-04-15-independent-blocks-fast-block-parallelism-implementation.md`
- Modify or create: benchmark driver under `benchmarks/`

**Step 1: Write production benchmark driver**

Benchmark production `runMCMC`, not only helper code.

Include:

- single-trait BayesC / BayesR
- annotated single-trait BayesC / BayesR
- plain multi-trait BayesC
- annotated 2-trait BayesC

Compare:

- `independent_blocks=false`
- `independent_blocks=true`

Use:

- multiple seeds
- longer chains
- supplied block starts
- at least two thread counts if the server permits

**Step 2: Enforce marker alignment in comparisons**

For marker-level output comparisons:

- normalize marker IDs
- join by explicit marker key
- assert expected marker count
- assert key sets match
- assert no duplicated or dropped markers
- save comparison-ready aligned tables

**Step 3: Run benchmarks**

Recommended environment:

```bash
export JULIA_NUM_THREADS=<num_cores>
export OPENBLAS_NUM_THREADS=1
julia --project=. benchmarks/<driver>.jl
```

**Step 4: Write report**

The report must separate:

- block-diagonal exactness checks
- realistic-data approximation behavior
- runtime and thread scaling
- practical recommendations for server use

**Step 5: Commit**

```bash
git add benchmarks docs/plans
git commit -m "bench: evaluate independent block parallelism"
```

### Task 10: Update manual and docs

**Files:**
- Modify: `docs/src/manual/workflow.md`
- Modify: `docs/src/manual/annotated_bayesc.md`
- Modify: `docs/src/manual/multitrait_annotated_bayesc.md`
- Modify: `src/1.JWAS/src/JWAS.jl`

**Step 1: Document API**

Document:

- `independent_blocks=false` as the exact default
- `independent_blocks=true` as an explicit approximate mode
- requirement that `fast_blocks != false`
- user-provided block starts through `fast_blocks=[...]`
- server environment recommendation:
  - `JULIA_NUM_THREADS`
  - `OPENBLAS_NUM_THREADS=1`

**Step 2: Document scientific interpretation**

Explain:

- exact if off-block weighted crossproducts are zero
- approximate otherwise
- pedigree can help propose blocks but does not replace genotype leakage checks

**Step 3: Build docs**

Run:

```bash
julia --project=docs --startup-file=no docs/make.jl
```

Expected: successful Documenter build.

**Step 4: Commit**

```bash
git add docs src/1.JWAS/src/JWAS.jl
git commit -m "docs: add independent block parallelism guidance"
```

### Task 11: Full verification

**Files:**
- No code changes expected

**Step 1: Run full tests**

Run:

```bash
julia --project=. --startup-file=no test/runtests.jl
```

Expected: all tests pass.

**Step 2: Run docs build**

Run:

```bash
julia --project=docs --startup-file=no docs/make.jl
```

Expected: successful build.

**Step 3: Inspect git state**

Run:

```bash
git status --short
```

Expected: no unintended tracked changes.
