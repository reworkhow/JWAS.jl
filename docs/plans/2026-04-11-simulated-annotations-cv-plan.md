# Simulated Annotations Cross-Validation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add production K-fold cross-validation benchmarking on the packaged `simulated_annotations` fixture for the current multi-trait BayesC-family methods plus the requested single-trait BayesC/BayesR baselines.

**Architecture:** Extend `benchmarks/simulated_annotations_multitrait_comparison.jl` with deterministic fold creation, phenotype masking, out-of-fold EBV scoring, and CV-specific summaries while keeping the existing non-CV marker-recovery path intact. Reuse the current `run_case` production flow and add focused tests around fold assignment and held-out EBV joins.

**Tech Stack:** Julia, DataFrames.jl, CSV.jl, JWAS production `get_genotypes`/`build_model`/`runMCMC` path, existing benchmark harness.

---

### Task 1: Add CV benchmark controls and fold utilities

**Files:**
- Modify: `/Users/haocheng/Github/JWAS.jl/.worktrees/mt-bayesc-sampler-benchmark/benchmarks/simulated_annotations_multitrait_comparison.jl`
- Test: `/Users/haocheng/Github/JWAS.jl/.worktrees/mt-bayesc-sampler-benchmark/test/unit/test_misc_coverage.jl`

**Step 1: Write the failing test**

Add coverage for:

- deterministic fold assignment for the same `(IDs, nfolds, seed)`
- no duplicate or missing IDs across folds
- each fold has non-empty held-out IDs

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: failure because the CV helper functions do not exist yet.

**Step 3: Write minimal implementation**

Add to the benchmark harness:

- `env_int` / `env_bool` reads for CV mode, fold count, and CV focus mode
- `cv_fold_assignments(ids, nfolds, seed)`
- validation that folds partition the full ID set exactly once

**Step 4: Run test to verify it passes**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: CV fold tests pass.

### Task 2: Add phenotype masking and held-out EBV scoring helpers

**Files:**
- Modify: `/Users/haocheng/Github/JWAS.jl/.worktrees/mt-bayesc-sampler-benchmark/benchmarks/simulated_annotations_multitrait_comparison.jl`
- Test: `/Users/haocheng/Github/JWAS.jl/.worktrees/mt-bayesc-sampler-benchmark/test/unit/test_misc_coverage.jl`

**Step 1: Write the failing test**

Add coverage for:

- masking both traits for held-out multi-trait individuals
- masking a single target trait for single-trait runs
- joining EBV output back to held-out individuals only with the expected row count

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: failure because masking/scoring helpers do not exist yet.

**Step 3: Write minimal implementation**

Add helpers such as:

- `masked_phenotype_frame(pheno_df, case, heldout_ids)`
- `heldout_prediction_row(case, seed, fold, trait, ebv_df, pheno_df)`
- held-out join assertions mirroring the existing marker/ID join checks

Compute:

- held-out `cor(y, EBV)`
- held-out RMSE
- held-out count

**Step 4: Run test to verify it passes**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: CV masking and held-out join tests pass.

### Task 3: Add a CV method matrix and benchmark execution path

**Files:**
- Modify: `/Users/haocheng/Github/JWAS.jl/.worktrees/mt-bayesc-sampler-benchmark/benchmarks/simulated_annotations_multitrait_comparison.jl`
- Test: `/Users/haocheng/Github/JWAS.jl/.worktrees/mt-bayesc-sampler-benchmark/test/unit/test_misc_coverage.jl`

**Step 1: Write the failing test**

Add coverage for a new CV method-selection mode that includes:

- `MT_BayesC`
- `MT_BayesC_I`
- `MT_BayesC_II`
- `MT_Annotated_BayesC_I`
- `MT_Annotated_BayesC_II`
- `MT_EmptyAnnotated_BayesC_I`
- `MT_EmptyAnnotated_BayesC_II`
- `BayesC_single`
- `Annotated_BayesC_single`
- `BayesR_single`
- `Annotated_BayesR_single`

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: failure because the CV focus mode is not implemented yet.

**Step 3: Write minimal implementation**

Add:

- a CV-specific method selector
- `run_case(...; pheno_input=...)` support so masked phenotypes can be passed into the production path
- a `run_cv_benchmark(...)` loop over seeds, folds, and methods

Save:

- `cv_fold_assignments.csv`
- `cv_per_fold_summary.csv`
- raw fold outputs under method/seed/fold directories

**Step 4: Run test to verify it passes**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: method-matrix CV selection tests pass.

### Task 4: Add CV aggregation and summary outputs

**Files:**
- Modify: `/Users/haocheng/Github/JWAS.jl/.worktrees/mt-bayesc-sampler-benchmark/benchmarks/simulated_annotations_multitrait_comparison.jl`
- Test: `/Users/haocheng/Github/JWAS.jl/.worktrees/mt-bayesc-sampler-benchmark/test/unit/test_misc_coverage.jl`

**Step 1: Write the failing test**

Add coverage for aggregation helpers that produce:

- per-method trait-level CV summaries
- family-level single-trait means across `y1` and `y2`
- runtime summaries across folds and seeds

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: failure because the CV summary helpers are missing.

**Step 3: Write minimal implementation**

Add:

- `summarize_cv_results(per_fold_df)`
- output files:
  - `cv_per_fold_summary.csv`
  - `cv_method_summary.csv`
  - `cv_family_summary.csv`

Summaries should include:

- `heldout_correlation_mean/sd`
- `heldout_rmse_mean/sd`
- `runtime_seconds_mean/sd`
- held-out count

**Step 4: Run test to verify it passes**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: aggregation tests pass.

### Task 5: Add a focused smoke CV benchmark

**Files:**
- Modify: `/Users/haocheng/Github/JWAS.jl/.worktrees/mt-bayesc-sampler-benchmark/test/unit/test_misc_coverage.jl`
- Modify: `/Users/haocheng/Github/JWAS.jl/.worktrees/mt-bayesc-sampler-benchmark/benchmarks/simulated_annotations_multitrait_comparison.jl`

**Step 1: Write the failing test**

Add a contract-style smoke test or script path expectation for CV mode output files.

**Step 2: Run the smoke benchmark to verify current failure**

Run:

```bash
JWAS_SIMULATED_MT_BENCHMARK_MODE=cv \
JWAS_SIMULATED_MT_CV_FOLDS=2 \
JWAS_SIMULATED_MT_SEEDS=11 \
JWAS_SIMULATED_MT_CHAIN_LENGTH=150 \
JWAS_SIMULATED_MT_BURNIN=50 \
JWAS_SIMULATED_MT_OUTPUT_FREQ=50 \
JWAS_SIMULATED_MT_WARMUP=false \
julia --project=. --startup-file=no benchmarks/simulated_annotations_multitrait_comparison.jl /tmp/jwas_simulated_annotations_cv_smoke
```

Expected: current failure or missing CV summary files.

**Step 3: Write minimal implementation**

Ensure CV mode writes the expected outputs and returns without invoking the
non-CV marker-recovery-only summaries as the primary result.

**Step 4: Run the smoke benchmark to verify it passes**

Run the same command.

Expected:

- benchmark completes
- `/tmp/jwas_simulated_annotations_cv_smoke/cv_per_fold_summary.csv` exists
- `/tmp/jwas_simulated_annotations_cv_smoke/cv_method_summary.csv` exists

### Task 6: Run the production CV benchmark

**Files:**
- Modify: `/Users/haocheng/Github/JWAS.jl/.worktrees/mt-bayesc-sampler-benchmark/benchmarks/simulated_annotations_multitrait_comparison.jl` if runtime tuning is needed
- Create: `/Users/haocheng/Github/JWAS.jl/.worktrees/mt-bayesc-sampler-benchmark/benchmarks/reports/2026-04-11-simulated-annotations-cv-report.md`

**Step 1: Run the production benchmark**

Suggested first full command:

```bash
JWAS_SIMULATED_MT_BENCHMARK_MODE=cv \
JWAS_SIMULATED_MT_CV_FOLDS=5 \
JWAS_SIMULATED_MT_SEEDS=101,202,303 \
JWAS_SIMULATED_MT_CHAIN_LENGTH=4000 \
JWAS_SIMULATED_MT_BURNIN=1000 \
JWAS_SIMULATED_MT_OUTPUT_FREQ=40 \
JWAS_SIMULATED_MT_WARMUP=true \
julia --project=. --startup-file=no benchmarks/simulated_annotations_multitrait_comparison.jl /tmp/jwas_simulated_annotations_cv_20260411
```

**Step 2: Inspect outputs**

Check:

- `cv_per_fold_summary.csv`
- `cv_method_summary.csv`
- `cv_family_summary.csv`

Verify row counts, trait coverage, and held-out counts before writing the
report.

**Step 3: Write the report**

Summarize:

- out-of-fold `cor(y, EBV)` by method and trait
- out-of-fold RMSE by method and trait
- method-family means
- runtime tradeoffs

### Task 7: Add implementation notes and verify

**Files:**
- Create: `/Users/haocheng/Github/JWAS.jl/.worktrees/mt-bayesc-sampler-benchmark/docs/plans/2026-04-11-simulated-annotations-cv-implementation.md`
- Verify: `/Users/haocheng/Github/JWAS.jl/.worktrees/mt-bayesc-sampler-benchmark/test/unit/test_misc_coverage.jl`
- Verify: `/Users/haocheng/Github/JWAS.jl/.worktrees/mt-bayesc-sampler-benchmark/test/unit/test_multitrait_mcmc.jl`

**Step 1: Run targeted tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: both pass.

**Step 2: Run the full suite**

Run:

```bash
julia --project=. --startup-file=no test/runtests.jl
```

Expected: full suite passes.

**Step 3: Write implementation note**

Record:

- commands run
- output directories
- main CV results
- any limitations or follow-up recommendations
