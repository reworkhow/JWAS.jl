# BayesC and BayesR Fast-Blocks Comparison Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Benchmark the production dense and `fast_blocks` paths for single-trait BayesC and BayesR under one long-chain multiseed protocol.

**Architecture:** Add a benchmark-only script that runs six production JWAS cases on one synthetic dataset, extracts comparable posterior summaries, and writes pairwise dense-vs-block comparison tables plus a markdown report. Keep production sampler code untouched.

**Tech Stack:** Julia, CSV.jl, DataFrames.jl, JWAS production `runMCMC`

---

### Task 1: Add a failing subprocess test for the new benchmark

**Files:**
- Modify: `test/unit/test_bayesr_parity.jl`

**Step 1: Write the failing test**

Add a subprocess test that runs the new benchmark script with a short multiseed
configuration and expects:

- `comparison_runs.csv`
- `comparison_pairwise_summary.csv`

It should also expect the six methods:

- `BayesC_dense`
- `BayesC_fast_blocks_default`
- `BayesC_fast_blocks_1`
- `BayesR_dense`
- `BayesR_fast_blocks_default`
- `BayesR_fast_blocks_1`

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected:

- failure because the benchmark script does not exist yet

### Task 2: Implement the benchmark script

**Files:**
- Create: `benchmarks/bayesc_bayesr_fast_blocks_comparison.jl`
- Modify: `benchmarks/bayesr_parity_common.jl` if a small shared helper is clearly reusable

**Step 1: Add the production run helper**

Implement a small helper that:

- builds a synthetic dataset
- runs `get_genotypes`
- builds the model
- runs `runMCMC`
- extracts posterior summaries for BayesC or BayesR

**Step 2: Add the six-case benchmark matrix**

Run:

- BayesC dense
- BayesC `fast_blocks=true`
- BayesC `fast_blocks=1`
- BayesR dense
- BayesR `fast_blocks=true`
- BayesR `fast_blocks=1`

Use the same seed list and same dataset for all cases.

**Step 3: Write output files**

Write:

- `comparison_runs.csv`
- `comparison_pairwise_summary.csv`

### Task 3: Re-run the subprocess test and full suite

**Files:**
- Test: `test/unit/test_bayesr_parity.jl`

**Step 1: Run targeted test**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected:

- new benchmark test passes

**Step 2: Run full suite**

Run:

```bash
julia --project=. --startup-file=no test/runtests.jl
```

Expected:

- full suite passes

### Task 4: Run the production benchmark and write the report

**Files:**
- Create: `benchmarks/reports/2026-03-20-bayesc-bayesr-fast-blocks-comparison-report.md`
- Create: `docs/plans/2026-03-20-bayesc-bayesr-fast-blocks-comparison-implementation.md`

**Step 1: Run the real benchmark**

Use a long-chain multiseed configuration through the new script and save the
CSV outputs under a temporary benchmark directory.

**Step 2: Write the report**

Summarize:

- BayesC dense vs block differences
- BayesR dense vs block differences
- whether `fast_blocks=1` is effectively identical to dense
- how default `fast_blocks=true` differs for each method

**Step 3: Write the implementation note**

Record:

- files added or modified
- benchmark protocol
- verification commands run

### Task 5: Commit

**Step 1: Commit the benchmark bundle**

Commit the benchmark script, test update, report, and implementation note after
verification.
