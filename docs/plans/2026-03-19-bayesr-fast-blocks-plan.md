# BayesR Fast Blocks Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add single-trait dense `fast_blocks` support for BayesR in JWAS and validate it against the existing dense BayesR path by posterior-summary parity.

**Architecture:** Reuse the current dense BayesR model and add a block-update kernel in `BayesR.jl`, mirroring the existing BayesABC block structure. The MCMC entry point and validation changes should stay narrow: enable `fast_blocks` for BayesR, dispatch to the block kernel, and keep all other current BayesR restrictions in place.

**Tech Stack:** Julia, JWAS MCMC code, BLAS, `Test`, benchmark scripts under `benchmarks/`, git.

---

### Task 1: Add Failing Validation Tests For BayesR Fast Blocks

**Files:**
- Modify: `test/unit/test_bayesr.jl`

**Step 1: Write the failing tests**

Add tests that currently fail because BayesR rejects `fast_blocks`:
- one test that `runMCMC(...; fast_blocks=true)` is accepted for single-trait dense BayesR
- one test that unsupported BayesR combinations still error:
  - `storage=:stream`
  - multi-trait
  - annotations
  - RRM

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected:
- FAIL because BayesR currently errors on `fast_blocks`

**Step 3: Commit**

```bash
git add test/unit/test_bayesr.jl
git commit -m "test: add BayesR fast_blocks validation coverage"
```

### Task 2: Enable BayesR Fast Blocks In Validation

**Files:**
- Modify: `src/1.JWAS/src/input_data_validation.jl`
- Test: `test/unit/test_bayesr.jl`

**Step 1: Write minimal validation change**

Remove only the BayesR hard rejection of `fast_blocks` while preserving the current BayesR errors for:
- `storage=:stream`
- multi-trait
- annotations
- RRM

**Step 2: Run test to verify validation behavior**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected:
- validation tests for BayesR `fast_blocks` now pass
- unsupported BayesR combinations still fail as intended

**Step 3: Commit**

```bash
git add src/1.JWAS/src/input_data_validation.jl test/unit/test_bayesr.jl
git commit -m "feat: allow BayesR fast_blocks validation path"
```

### Task 3: Write BayesR Block Kernel Tests First

**Files:**
- Modify: `test/unit/test_bayesr.jl`

**Step 1: Write the failing tests**

Add focused tests for the new block path:
- BayesR single-trait run completes with `fast_blocks=true`
- output contains BayesR `pi` with 4 rows
- marker effects output exists
- class labels remain in `1:4`

Add one parity-style test on a small deterministic dataset:
- run dense BayesR and block BayesR with the same seed
- compare posterior-summary outputs with reasonable tolerances, not exact identity

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected:
- FAIL because BayesR has no block kernel yet

**Step 3: Commit**

```bash
git add test/unit/test_bayesr.jl
git commit -m "test: add BayesR fast_blocks kernel coverage"
```

### Task 4: Implement BayesR Block Wrapper And Kernel

**Files:**
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`
- Reference: `src/1.JWAS/src/markers/BayesianAlphabet/BayesABC.jl`

**Step 1: Add the wrapper**

Add:
- `BayesR_block!(genotypes, ycorr, vare, Rinv=...)`

The wrapper should pass:
- `genotypes.MArray`
- `genotypes.mpRinvm`
- `genotypes.genotypes`
- `genotypes.MpRinvM`
- `genotypes.α[1]`
- `genotypes.δ[1]`
- `genotypes.G.val`
- `genotypes.π`
- `BAYESR_GAMMA`

**Step 2: Add the block kernel**

Add the lower-level block function mirroring the BayesABC block path:
- cache `XpRinvycorr` per block with `block_rhs!`
- save `αold_block`
- use `nreps = block_size`
- for each marker in the block:
  - compute BayesR `rhs`
  - evaluate class log probabilities
  - sample class label
  - sample nonzero `α`
  - update cached block RHS
- at the end of the block, push the final block correction to `yCorr` with `mul!`

**Step 3: Run test to verify block kernel behavior**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected:
- targeted BayesR tests pass

**Step 4: Commit**

```bash
git add src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl test/unit/test_bayesr.jl
git commit -m "feat: add BayesR fast_blocks kernel"
```

### Task 5: Dispatch BayesR To The Block Path

**Files:**
- Modify: `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`
- Test: `test/unit/test_bayesr.jl`

**Step 1: Add the block dispatch**

Update the single-trait BayesR branch so that:
- dense BayesR uses `BayesR!` when `fast_blocks == false`
- dense BayesR uses `BayesR_block!` when `fast_blocks != false`

Do not change:
- multi-trait handling
- `pi` updating
- marker-variance updating

**Step 2: Run targeted tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected:
- BayesR dense and block tests pass

**Step 3: Commit**

```bash
git add src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl test/unit/test_bayesr.jl
git commit -m "feat: dispatch BayesR to fast_blocks path"
```

### Task 6: Add Dense-Vs-Block Posterior-Summary Benchmark

**Files:**
- Create or Modify: `benchmarks/bayesr_fast_blocks_parity.jl`
- Create: `benchmarks/reports/2026-03-19-bayesr-fast-blocks-parity.md`

**Step 1: Write the benchmark script**

The script should:
- build one small or moderate dense BayesR dataset
- run dense BayesR
- run block BayesR with the same seed and chain settings
- compare posterior summaries:
  - `sigmaSq`
  - residual variance
  - `pi`
  - nonzero frequency
  - marker-effect posterior mean correlation

**Step 2: Run the benchmark**

Run:

```bash
julia --project=. --startup-file=no benchmarks/bayesr_fast_blocks_parity.jl
```

Expected:
- script completes and writes a concise parity report

**Step 3: Commit**

```bash
git add benchmarks/bayesr_fast_blocks_parity.jl benchmarks/reports/2026-03-19-bayesr-fast-blocks-parity.md
git commit -m "benchmarks: add BayesR fast_blocks parity report"
```

### Task 7: Run Full Verification

**Files:**
- No code changes unless verification fails

**Step 1: Run targeted BayesR tests**

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected:
- PASS

**Step 2: Run full test suite**

```bash
julia --project=. --startup-file=no test/runtests.jl
```

Expected:
- PASS

**Step 3: If docs changed, build docs**

```bash
julia --project=docs --startup-file=no docs/make.jl
```

Expected:
- PASS

**Step 4: Commit any verification-driven fixes**

```bash
git add <relevant files>
git commit -m "fix: address BayesR fast_blocks verification issues"
```

### Task 8: Write Final Implementation Record

**Files:**
- Create: `docs/plans/2026-03-19-bayesr-fast-blocks-implementation.md`

**Step 1: Write the implementation record**

Document:
- final architecture
- files changed
- any deviations from the approved design
- benchmark/parity outcome
- verification commands and results

**Step 2: Commit**

```bash
git add docs/plans/2026-03-19-bayesr-fast-blocks-implementation.md
git commit -m "docs: add BayesR fast_blocks implementation record"
```
