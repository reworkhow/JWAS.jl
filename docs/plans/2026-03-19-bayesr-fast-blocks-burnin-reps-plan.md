# BayesR Fast Blocks Burnin Repetition Schedule Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make BayesR `fast_blocks` use one within-block sweep per outer iteration during burnin, then switch to `nreps = block_size` after burnin.

**Architecture:** Keep the behavior change BayesR-specific and production-local. Compute the desired repetition count in the MCMC loop, pass it into `BayesR_block!`, and add regression coverage that proves the schedule is applied correctly.

**Tech Stack:** Julia, JWAS MCMC code, existing `Test` framework.

---

### Task 1: Add scheduling regression tests

**Files:**
- Modify: `test/unit/test_bayesr.jl`
- Reference: `src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`
- Reference: `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`

**Step 1: Write the failing test**

Add a unit test for a small BayesR block case that checks:

- BayesR block scheduling uses `nreps = 1` during burnin
- BayesR block scheduling uses `nreps = block_size` after burnin

Prefer testing the scheduling helper directly if one is introduced.

**Step 2: Run the targeted test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected: the new test fails because the current code always uses `nreps = block_size`.

**Step 3: Commit the failing test only if useful**

If you want strict TDD history:

```bash
git add test/unit/test_bayesr.jl
git commit -m "test: add BayesR fast-block burnin scheduling regression"
```

Otherwise continue directly.

### Task 2: Implement the BayesR repetition schedule

**Files:**
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`
- Modify: `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`

**Step 1: Add minimal scheduling logic**

Introduce a small helper or inline logic in `MCMC_BayesianAlphabet.jl` that computes:

- `block_reps = 1` during effective burnin
- `block_reps = block_size` after effective burnin

for BayesR fast-block runs only.

**Step 2: Pass the repetition count into `BayesR_block!`**

Update the BayesR block call path so `BayesR_block!` takes an explicit repetition count.

**Step 3: Use the passed repetition count in the block kernel**

Replace the current hardcoded:

```julia
nreps = block_size
```

with the passed value.

**Step 4: Keep non-BayesR behavior unchanged**

Do not modify:

- dense BayesR
- BayesABC block behavior
- public `runMCMC` API

### Task 3: Verify targeted behavior

**Files:**
- Modify if needed: `test/unit/test_bayesr.jl`

**Step 1: Run the targeted BayesR tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected: PASS.

**Step 2: If the test design is weak, strengthen it**

If the first test only checks helper output, add a second regression that exercises the actual BayesR block call path enough to prove the repetition count is consumed.

### Task 4: Full verification

**Files:**
- No new files unless failure triage requires edits

**Step 1: Run the full suite**

Run:

```bash
julia --project=. --startup-file=no test/runtests.jl
```

Expected: PASS with no regressions.

**Step 2: Summarize remaining local-only benchmark edits**

Before claiming completion, explicitly report whether these earlier exploratory files are still dirty:

- `benchmarks/bayesr_fast_blocks_parity.jl`
- `benchmarks/bayesr_parity_common.jl`
- `test/unit/test_bayesr_parity.jl`

### Task 5: Record implementation and commit

**Files:**
- Create: `docs/plans/2026-03-19-bayesr-fast-blocks-burnin-reps-implementation.md`
- Modify: any production/test files changed above

**Step 1: Write the implementation record**

Document:

- what changed
- the final scheduling rule
- tests run
- whether any benchmark files were intentionally left uncommitted

**Step 2: Commit the narrow production change**

Example:

```bash
git add src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl \
        src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl \
        test/unit/test_bayesr.jl \
        docs/plans/2026-03-19-bayesr-fast-blocks-burnin-reps-implementation.md
git commit -m "feat: gate BayesR fast-block repetitions by burnin"
```
