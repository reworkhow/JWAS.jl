# BayesR Fast Blocks Fixed-Hyperparameters Diagnostic Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extend the BayesR `fast_blocks` benchmark so it can compare dense vs block BayesR with fixed `pi` and fixed `sigmaSq`, and quantify baseline dense-BayesR seed-to-seed variation under the same fixed-hyperparameter setup.

**Architecture:** Keep all changes in the benchmark layer. Reuse `benchmarks/bayesr_fast_blocks_parity.jl` and `benchmarks/bayesr_parity_common.jl`, and use existing JWAS flags (`estimatePi=false`, fixed marker variance input, `estimate_variance=false`) rather than adding benchmark-only controls to production sampler code.

**Tech Stack:** Julia, JWAS, CSV.jl, DataFrames.jl, benchmark scripts under `benchmarks/`, parity helper utilities, `Test`.

---

### Task 1: Add Benchmark Helper Coverage For The New Diagnostic Modes

**Files:**
- Modify: `test/unit/test_bayesr_parity.jl`
- Reference: `benchmarks/bayesr_parity_common.jl`

**Step 1: Write the failing tests**

Add helper-level tests for the new benchmark outputs:

- a test that runtime metadata can record requested vs effective burnin and block
  size
- a test that a fixed-hyperparameter summary can be written and read without
  relying on estimated `pi`
- a test that dense multiseed summary aggregation works for within-method runs

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected:

- FAIL because the new helper behavior does not exist yet

**Step 3: Commit**

```bash
git add test/unit/test_bayesr_parity.jl
git commit -m "test: add BayesR fixed-hyperparameter benchmark coverage"
```

### Task 2: Extend The Benchmark Helpers For Fixed-Hyperparameter Reporting

**Files:**
- Modify: `benchmarks/bayesr_parity_common.jl`
- Test: `test/unit/test_bayesr_parity.jl`

**Step 1: Add minimal helper support**

Extend helper functions as needed so the benchmark can:

- write/read scalar summaries for fixed-hyperparameter runs
- normalize fixed `pi` output cleanly
- summarize within-method multiseed runs without assuming a reference method

Keep the helper API narrow and benchmark-focused.

**Step 2: Run tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected:

- new helper tests pass

**Step 3: Commit**

```bash
git add benchmarks/bayesr_parity_common.jl test/unit/test_bayesr_parity.jl
git commit -m "feat: extend BayesR benchmark helpers for fixed-hyperparameter diagnostics"
```

### Task 3: Add Fixed-Hyperparameter Dense-Vs-Block Mode

**Files:**
- Modify: `benchmarks/bayesr_fast_blocks_parity.jl`
- Reference: `src/1.JWAS/src/JWAS.jl`

**Step 1: Add the benchmark mode**

Extend the script so it can run:

- dense BayesR with fixed `pi` and fixed `sigmaSq`
- block BayesR with fixed `pi` and fixed `sigmaSq`

Use existing JWAS controls:

- `estimatePi=false`
- `estimate_variance=false` for marker variance
- `Pi=start_pi`
- marker variance input = `start_sigma_sq`

Residual variance should still be estimated.

**Step 2: Keep the current block-aware burnin handling**

Retain the benchmark-only burnin adjustment for block runs so the script always
produces post-burnin samples.

**Step 3: Run the script once**

Run:

```bash
julia --project=. --startup-file=no benchmarks/bayesr_fast_blocks_parity.jl
```

Expected:

- fixed-hyperparameter dense and block summaries are written
- no all-zero summary files caused by burnin/output mismatch

**Step 4: Commit**

```bash
git add benchmarks/bayesr_fast_blocks_parity.jl
git commit -m "feat: add fixed-hyperparameter BayesR fast_blocks benchmark mode"
```

### Task 4: Add Dense Multiseed Baseline Mode

**Files:**
- Modify: `benchmarks/bayesr_fast_blocks_parity.jl`
- Possibly modify: `benchmarks/bayesr_parity_common.jl`
- Test: `test/unit/test_bayesr_parity.jl`

**Step 1: Add multiseed dense-only execution**

Extend the benchmark so it can:

- run dense BayesR only
- reuse one dataset
- vary only the seed
- write:
  - `multiseed_runs.csv`
  - `multiseed_summary.csv`

Keep this in the same script unless that becomes messy.

**Step 2: Run the helper tests again**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected:

- helper tests still pass

**Step 3: Commit**

```bash
git add benchmarks/bayesr_fast_blocks_parity.jl benchmarks/bayesr_parity_common.jl test/unit/test_bayesr_parity.jl
git commit -m "feat: add dense multiseed baseline for BayesR fast_blocks benchmark"
```

### Task 5: Run The Diagnostic Benchmark And Save The Report

**Files:**
- Run: `benchmarks/bayesr_fast_blocks_parity.jl`
- Create: `benchmarks/reports/2026-03-19-bayesr-fast-blocks-fixed-hyperparameters.md`

**Step 1: Run fixed-hyperparameter dense-vs-block**

Run:

```bash
julia --project=. --startup-file=no benchmarks/bayesr_fast_blocks_parity.jl
```

Capture:

- residual variance differences
- mean nonzero frequency differences
- marker-effect summary differences
- runtime

**Step 2: Run dense multiseed baseline**

Run the same script in dense-multiseed mode with a small seed set, for example:

```bash
julia --project=. --startup-file=no benchmarks/bayesr_fast_blocks_parity.jl /tmp/bayesr_fast_blocks_diag
```

or the exact mode-specific invocation added in Task 4.

Capture:

- within-dense seed variation for the same summaries

**Step 3: Write the report**

Document:

- benchmark setup
- fixed `pi` / fixed `sigmaSq` configuration
- dense-vs-block results
- dense-vs-dense multiseed variation
- conclusion:
  - block-path problem
  - or mainly `pi`/`sigmaSq` feedback problem

**Step 4: Commit**

```bash
git add benchmarks/reports/2026-03-19-bayesr-fast-blocks-fixed-hyperparameters.md benchmarks/bayesr_fast_blocks_parity.jl
git commit -m "benchmarks: add BayesR fixed-hyperparameter fast_blocks diagnostic report"
```

### Task 6: Write The Implementation Record And Run Final Verification

**Files:**
- Create: `docs/plans/2026-03-19-bayesr-fast-blocks-fixed-hyperparameters-implementation.md`

**Step 1: Write the implementation record**

Include:

- files changed
- benchmark modes added
- key results
- deviations from the earlier fast-blocks parity workflow
- whether the diagnosis points to the block kernel or to the `pi`/`sigmaSq`
  feedback loop

**Step 2: Run verification**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
julia --project=. --startup-file=no test/runtests.jl
```

Expected:

- benchmark-helper tests pass
- full suite passes

**Step 3: Commit**

```bash
git add docs/plans/2026-03-19-bayesr-fast-blocks-fixed-hyperparameters-implementation.md
git commit -m "docs: add BayesR fixed-hyperparameter diagnostic implementation record"
```
