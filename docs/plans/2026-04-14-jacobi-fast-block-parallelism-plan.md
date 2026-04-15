# Jacobi Fast-Block Parallelism Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add an explicit `block_sweep=:gauss_seidel|:jacobi` option so JWAS can run approximate independent-block Jacobi sweeps for fast-block single-trait and multi-trait BayesianAlphabet samplers, with user-supplied block partitions and thread-parallel execution.

**Architecture:** Keep the current fast-block sweep as the exact default (`:gauss_seidel`) and add a second explicit path (`:jacobi`) that freezes the sweep-level corrected phenotype / residual, samples blocks independently, and applies block deltas at the barrier. Reuse existing within-block BayesABC, BayesR, and MTBayesABC logic as much as possible; change only the inter-block sweep order and the API plumbing needed to expose the new mode.

**Tech Stack:** Julia, existing JWAS BayesianAlphabet samplers, `Threads.@threads`, existing `fast_blocks` / block matrix infrastructure, `test/runtests.jl`, Documenter.

---

### Task 1: Add the new block-sweep option to the public MCMC API

**Files:**
- Modify: `src/1.JWAS/src/JWAS.jl`
- Modify: `src/1.JWAS/src/types.jl`
- Test: `test/unit/test_misc_coverage.jl`

**Step 1: Write the failing API tests**

Add tests that expect:

- `runMCMC(...; block_sweep=:gauss_seidel)` to be accepted
- `runMCMC(...; block_sweep=:jacobi)` to be accepted
- invalid symbols such as `:foo` to error clearly

**Step 2: Run the targeted test file and confirm failure**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: failure because `block_sweep` is not a recognized keyword yet.

**Step 3: Add minimal API plumbing**

Implement:

- `block_sweep::Symbol = :gauss_seidel` in both public `runMCMC` signatures
- `block_sweep` field on `MCMCinfo`
- validation that only `:gauss_seidel` and `:jacobi` are accepted

**Step 4: Run the targeted test again**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: the new API tests pass.

**Step 5: Commit**

```bash
git add src/1.JWAS/src/JWAS.jl src/1.JWAS/src/types.jl test/unit/test_misc_coverage.jl
git commit -m "feat: add block_sweep MCMC option"
```

### Task 2: Allow explicit user-provided block partitions

**Files:**
- Modify: `src/1.JWAS/src/JWAS.jl`
- Modify: `src/1.JWAS/src/input_data_validation.jl`
- Test: `test/unit/test_misc_coverage.jl`

**Step 1: Write the failing partition tests**

Add tests that expect:

- `fast_blocks=[1, 101, 201]` to be accepted
- unsorted vectors to error
- vectors not starting at `1` to error
- duplicate starts to error

**Step 2: Run the targeted test file and confirm failure**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: failure because `fast_blocks` currently only accepts `false`, `true`, or a number.

**Step 3: Implement explicit block-start support**

In `runMCMC(...)`, extend the parsing logic:

- `fast_blocks === true` keeps the current heuristic
- `fast_blocks isa Number` keeps the current fixed-width logic
- `fast_blocks isa AbstractVector{<:Integer}` uses the supplied block starts directly after validation

Also validate:

- sorted ascending
- unique
- first start is `1`
- all starts are within marker bounds

**Step 4: Run the targeted test again**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: the partition tests pass.

**Step 5: Commit**

```bash
git add src/1.JWAS/src/JWAS.jl src/1.JWAS/src/input_data_validation.jl test/unit/test_misc_coverage.jl
git commit -m "feat: allow explicit fast block partitions"
```

### Task 3: Refactor single-trait BayesABC block code to separate sweep policy from within-block updates

**Files:**
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/BayesABC.jl`
- Test: `test/unit/test_misc_coverage.jl`

**Step 1: Write a failing regression test for block-sweep dispatch**

Add a small test that:

- exercises `BayesABC_block!` with `block_sweep=:gauss_seidel`
- exercises `BayesABC_block!` with `block_sweep=:jacobi`
- confirms both paths execute on a toy example

**Step 2: Run the targeted test and confirm failure**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: failure because `BayesABC_block!` has only one sweep policy today.

**Step 3: Extract the current sweep into an explicit Gauss-Seidel helper**

Refactor the current implementation into something like:

- `_BayesABC_block_gauss_seidel!(...)`
- `_BayesABC_block_jacobi!(...)`

Keep the current code path behaviorally unchanged in the Gauss-Seidel helper.

**Step 4: Add the Jacobi helper without threading first**

Implement a first Jacobi version that:

- freezes `yCorr_old` at the start of the sweep
- computes each block's local RHS from that frozen snapshot
- updates each block independently
- applies the combined block delta after all blocks finish

Do this sequentially first; keep threading for a later task.

**Step 5: Run the targeted tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: both Gauss-Seidel and Jacobi paths execute on the toy case.

**Step 6: Commit**

```bash
git add src/1.JWAS/src/markers/BayesianAlphabet/BayesABC.jl test/unit/test_misc_coverage.jl
git commit -m "feat: add jacobi sweep for BayesABC blocks"
```

### Task 4: Refactor single-trait BayesR block code the same way

**Files:**
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`
- Test: `test/unit/test_misc_coverage.jl`
- Test: `test/unit/test_annotated_bayesr.jl`

**Step 1: Write failing BayesR block-sweep tests**

Add tests that cover:

- plain BayesR block Gauss-Seidel
- plain BayesR block Jacobi
- annotated BayesR block Jacobi

**Step 2: Run the targeted tests and confirm failure**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl"); include("test/unit/test_annotated_bayesr.jl")'
```

Expected: failure because BayesR block code does not understand `block_sweep`.

**Step 3: Extract sweep policy helpers**

Create:

- `_BayesR_block_gauss_seidel!(...)`
- `_BayesR_block_jacobi!(...)`

Reuse the existing within-block BayesR class-update logic. Only change the inter-block sweep order.

**Step 4: Run the targeted tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl"); include("test/unit/test_annotated_bayesr.jl")'
```

Expected: BayesR block Jacobi tests pass.

**Step 5: Commit**

```bash
git add src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl test/unit/test_misc_coverage.jl test/unit/test_annotated_bayesr.jl
git commit -m "feat: add jacobi sweep for BayesR blocks"
```

### Task 5: Refactor multi-trait BayesC block code to support Jacobi sweeps

**Files:**
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl`
- Modify: `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`
- Test: `test/unit/test_multitrait_mcmc.jl`
- Test: `test/unit/test_annotated_bayesc.jl`

**Step 1: Write failing multi-trait block-sweep tests**

Add tests covering:

- plain multi-trait BayesC block sampler I with `block_sweep=:jacobi`
- plain multi-trait BayesC block sampler II with `block_sweep=:jacobi`
- annotated 2-trait BayesC block sampler I with `block_sweep=:jacobi`
- annotated 2-trait BayesC block sampler II with `block_sweep=:jacobi`

**Step 2: Run the targeted tests and confirm failure**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl"); include("test/unit/test_annotated_bayesc.jl")'
```

Expected: failure because `MTBayesABC_block!` currently only has the Gauss-Seidel sweep.

**Step 3: Split the current block code by sweep policy**

Refactor `MTBayesABC_block!` into:

- `_MTBayesABC_block_samplerI_gauss_seidel!(...)`
- `_MTBayesABC_block_samplerI_jacobi!(...)`
- `_MTBayesABC_block_samplerII_gauss_seidel!(...)`
- `_MTBayesABC_block_samplerII_jacobi!(...)`

Reuse:

- existing `GlobalPiPrior` / `MarkerSpecificPiPrior`
- existing sampler I vs sampler II math
- existing `mt_bayesc_sampler_mode(...)`

Only change the inter-block sweep order.

**Step 4: Implement Jacobi residual freezing**

For the Jacobi helpers:

- freeze the trait-wise corrected phenotype vectors at sweep start
- compute each block's local RHS from that frozen snapshot
- store block-local effect deltas
- apply the combined delta after all blocks finish

**Step 5: Run the targeted tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl"); include("test/unit/test_annotated_bayesc.jl")'
```

Expected: plain and annotated multi-trait Jacobi tests pass.

**Step 6: Commit**

```bash
git add src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl test/unit/test_multitrait_mcmc.jl test/unit/test_annotated_bayesc.jl
git commit -m "feat: add jacobi sweep for multitrait BayesC blocks"
```

### Task 6: Add thread-parallel execution for Jacobi block sweeps

**Files:**
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/BayesABC.jl`
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl`
- Test: `test/unit/test_misc_coverage.jl`
- Test: `test/unit/test_multitrait_mcmc.jl`

**Step 1: Add a failing test for deterministic Jacobi execution shape**

Write a small regression that checks:

- Jacobi mode executes with multiple blocks and returns valid updated states
- the result is not sensitive to block execution order on an exactly block-diagonal toy problem

This can be implemented by comparing:

- sequential Jacobi
- threaded Jacobi

on the same toy data and same RNG seeding scheme.

**Step 2: Run the targeted tests and confirm failure**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl"); include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: failure because Jacobi is still sequential.

**Step 3: Implement thread-parallel block execution**

For Jacobi helpers only:

- use `Threads.@threads` over blocks
- allocate thread-local buffers
- use thread-local RNG state
- return block-local deltas
- reduce those deltas once at the barrier

Set the implementation expectation in code comments:

- users should set `OPENBLAS_NUM_THREADS=1`
- parallelism is across blocks, not through BLAS

**Step 4: Run the targeted tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl"); include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: threaded Jacobi tests pass.

**Step 5: Commit**

```bash
git add src/1.JWAS/src/markers/BayesianAlphabet/BayesABC.jl src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl test/unit/test_misc_coverage.jl test/unit/test_multitrait_mcmc.jl
git commit -m "feat: parallelize jacobi block sweeps"
```

### Task 7: Add exactness tests on block-diagonal synthetic data

**Files:**
- Modify: `test/unit/test_multitrait_mcmc.jl`
- Modify: `test/unit/test_misc_coverage.jl`

**Step 1: Write exactness tests**

Add synthetic block-diagonal cases where off-block crossproducts are exactly zero:

- single-trait BayesC/BayesR
- multi-trait BayesC

Test that Gauss-Seidel and Jacobi agree up to Monte Carlo tolerance on:

- posterior state frequencies in tiny exact cases
- or identical block-local updates under controlled RNG when the block conditionals decouple exactly

**Step 2: Run the tests and confirm failure first if needed**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl"); include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: if not already satisfied, fail until the synthetic exactness checks are wired correctly.

**Step 3: Finalize the exactness regressions**

Keep the tests small and deterministic. These are the core mathematical safety checks for the new mode.

**Step 4: Run the tests again**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl"); include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: exactness tests pass.

**Step 5: Commit**

```bash
git add test/unit/test_misc_coverage.jl test/unit/test_multitrait_mcmc.jl
git commit -m "test: add exactness checks for jacobi block sweeps"
```

### Task 8: Benchmark approximation quality and server scaling

**Files:**
- Create: `benchmarks/reports/2026-04-14-jacobi-fast-block-parallelism-report.md`
- Create: `docs/plans/2026-04-14-jacobi-fast-block-parallelism-implementation.md`
- Modify or Create: benchmark driver under `benchmarks/`

**Step 1: Write the benchmark driver**

Benchmark the production JWAS path, not helper functions only.

Include:

- single-trait BayesC / BayesR
- multi-trait BayesC plain / annotated where relevant
- `block_sweep=:gauss_seidel`
- `block_sweep=:jacobi`

Use:

- realistic user-supplied block partitions
- multiple seeds
- longer chains, not a single short run

**Step 2: Run the production benchmarks**

Record:

- runtime
- strong scaling across thread counts
- EBV / prediction summaries
- key posterior summaries

**Step 3: Write the report**

In the report, explicitly separate:

- exactness on block-diagonal synthetic data
- approximation behavior on realistic data
- server scaling behavior

Also document:

- that the main speedup comes from parallel block sweeps
- that Julia-vs-OpenMP is secondary to the algorithmic change

**Step 4: Save implementation notes**

Record:

- commands
- environment variables (`JULIA_NUM_THREADS`, `OPENBLAS_NUM_THREADS`)
- seed sets
- benchmark hardware assumptions

**Step 5: Commit**

```bash
git add benchmarks docs/plans
git commit -m "bench: evaluate jacobi fast block parallelism"
```

### Task 9: Update manual and API docs

**Files:**
- Modify: `docs/src/manual/workflow.md`
- Modify: `docs/src/manual/annotated_bayesc.md`
- Modify: `docs/src/manual/multitrait_annotated_bayesc.md`
- Modify: docstrings in `src/1.JWAS/src/JWAS.jl`

**Step 1: Document the new option**

Explain:

- `block_sweep=:gauss_seidel|:jacobi`
- user-provided block partitions through `fast_blocks`
- exact vs approximate interpretation
- recommendation to use user-chosen contiguous blocks
- pedigree-informed partitioning as external guidance, not built-in functionality

**Step 2: Run the docs build**

Run:

```bash
julia --project=docs --startup-file=no docs/make.jl
```

Expected: successful Documenter build.

**Step 3: Commit**

```bash
git add docs src/1.JWAS/src/JWAS.jl
git commit -m "docs: add jacobi fast block parallelism guidance"
```

### Task 10: Run full regression verification

**Files:**
- No code changes expected

**Step 1: Run the full test suite**

Run:

```bash
julia --project=. --startup-file=no test/runtests.jl
```

Expected: all tests pass.

**Step 2: Run the docs build again**

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

Expected: only intended tracked changes remain.

**Step 4: Final verification commit if needed**

If any last-minute fixes were made after earlier commits:

```bash
git add -A
git commit -m "chore: finalize jacobi fast block parallelism"
```
