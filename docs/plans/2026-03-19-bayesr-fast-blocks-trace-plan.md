# BayesR Fast Blocks Trace Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a short fixed-hyperparameter trace benchmark for dense vs BayesR `fast_blocks` so the remaining residual-variance mismatch can be localized.

**Architecture:** Extend the existing `benchmarks/bayesr_fast_blocks_parity.jl` driver with a new trace mode that uses a benchmark-local outer Gibbs loop, fixed `pi`, fixed `sigmaSq`, and the production BayesR dense/block marker kernels. Save dense and block per-iteration traces plus one aligned comparison file.

**Tech Stack:** Julia, JWAS BayesR kernels, CSV, DataFrames, benchmark reports under `benchmarks/reports/`.

---

### Task 1: Add Trace Helper Tests First

**Files:**
- Modify: `test/unit/test_bayesr_parity.jl`

**Step 1: Write failing tests**

Add coverage for the trace benchmark helpers you will need:

- trace row writer / schema
- trace comparison table generator

The tests should check:

- expected trace columns exist
- dense/block comparison aligns by iteration
- absolute difference columns are computed correctly

**Step 2: Run targeted test**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected:
- FAIL because the new helper functions do not exist yet

**Step 3: Commit**

```bash
git add test/unit/test_bayesr_parity.jl
git commit -m "test: add BayesR fast_blocks trace helper coverage"
```

### Task 2: Add Trace Helper Functions

**Files:**
- Modify: `benchmarks/bayesr_parity_common.jl`
- Test: `test/unit/test_bayesr_parity.jl`

**Step 1: Add the minimal helpers**

Add helpers for:

- writing trace rows to a `DataFrame`
- building an aligned dense-vs-block comparison table

Keep the helpers generic and benchmark-local.

**Step 2: Run targeted test**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected:
- PASS for the new helper coverage

**Step 3: Commit**

```bash
git add benchmarks/bayesr_parity_common.jl test/unit/test_bayesr_parity.jl
git commit -m "feat: add BayesR trace benchmark helpers"
```

### Task 3: Add Fixed-Hyperparameters Trace Mode

**Files:**
- Modify: `benchmarks/bayesr_fast_blocks_parity.jl`

**Step 1: Add the new mode**

Extend the benchmark driver with:

- `JWAS_BAYESR_BLOCK_MODE=fixed_hyperparameters_trace`

The mode should:

- build the same synthetic dataset as the current block benchmark
- fix `pi`
- fix `sigmaSq`
- use `chain_length = 100`, `burnin = 0` by default
- run dense and block on the same seed

**Step 2: Add the benchmark-local Gibbs loop**

Implement a local outer loop that:

- calls the production dense or block BayesR marker kernel
- keeps `pi` fixed
- keeps `sigmaSq` fixed
- updates residual variance locally
- records one trace row per iteration

Trace columns:

- `iter`
- `residual_variance`
- `ycorr_norm`
- `alpha_norm`
- `alpha_abs_mean`
- `nnz`
- `max_abs_alpha`

**Step 3: Write outputs**

Write:

- `dense_trace.csv`
- `fast_blocks_trace.csv`
- `trace_comparison.csv`

**Step 4: Run the benchmark**

Run:

```bash
julia --project=. --startup-file=no benchmarks/bayesr_fast_blocks_parity.jl /tmp/bayesr_fast_blocks_trace_run
```

Expected:
- benchmark completes
- trace files are written

**Step 5: Commit**

```bash
git add benchmarks/bayesr_fast_blocks_parity.jl
git commit -m "benchmarks: add BayesR fast_blocks trace mode"
```

### Task 4: Write Trace Report

**Files:**
- Create: `benchmarks/reports/2026-03-19-bayesr-fast-blocks-trace.md`
- Create: `docs/plans/2026-03-19-bayesr-fast-blocks-trace-implementation.md`

**Step 1: Summarize the benchmark**

Document:

- setup
- trace fields
- first visible divergence point
- interpretation of whether drift appears first in:
  - `yCorr`
  - `alpha`
  - or `residual_variance`

**Step 2: Save the implementation record**

Document:

- what changed
- what was verified
- what the trace concluded

**Step 3: Commit**

```bash
git add benchmarks/reports/2026-03-19-bayesr-fast-blocks-trace.md docs/plans/2026-03-19-bayesr-fast-blocks-trace-implementation.md
git commit -m "docs: record BayesR fast_blocks trace diagnostic"
```

### Task 5: Run Full Verification

**Files:**
- No code changes unless verification fails

**Step 1: Run targeted parity tests**

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected:
- PASS

**Step 2: Run the trace benchmark**

```bash
julia --project=. --startup-file=no benchmarks/bayesr_fast_blocks_parity.jl /tmp/bayesr_fast_blocks_trace_run
```

Expected:
- PASS

**Step 3: Run full test suite**

```bash
julia --project=. --startup-file=no test/runtests.jl
```

Expected:
- PASS

**Step 4: Commit if verification changes were needed**

If verification required code or report changes:

```bash
git add <files>
git commit -m "test: finalize BayesR fast_blocks trace diagnostic"
```
