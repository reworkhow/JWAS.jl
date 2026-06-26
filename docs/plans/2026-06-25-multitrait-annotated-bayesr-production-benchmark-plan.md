# Multi-Trait Annotated BayesR Production Benchmark Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extend the simulated-annotations benchmark so it evaluates causal-variant discovery accuracy and cross-chain stability for multi-trait annotated BayesR against practical baselines.

**Architecture:** Reuse `benchmarks/simulated_annotations_multitrait_comparison.jl` because it already runs the production JWAS path and joins marker outputs to truth by explicit marker IDs. Add metric helpers for average precision, top-k precision/count, and seed-to-seed stability from the joined marker tables. Keep raw MCMC output under `benchmarks/out/` and save the production report under `benchmarks/reports/`.

**Tech Stack:** Julia, JWAS, CSV.jl, DataFrames.jl, Statistics, bundled `simulated_annotations` dataset.

---

### Task 1: Plan Record

**Files:**
- Create: `docs/plans/2026-06-25-multitrait-annotated-bayesr-production-benchmark-plan.md`

**Steps:**
1. Write this implementation plan.
2. Commit it with the later benchmark/report changes after verification.

### Task 2: Test New Metric Contracts

**Files:**
- Modify: `test/unit/test_misc_coverage.jl`

**Steps:**
1. Add tests for `average_precision`, `top_k_precision`, and the new benchmark focus mode.
2. Add tests for a small two-seed marker-stability summary using hand-built marker tables.
3. Run `julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'`.
4. Expected first result: FAIL because the new metric helpers and stability function do not exist yet.

### Task 3: Implement Metrics

**Files:**
- Modify: `benchmarks/simulated_annotations_multitrait_comparison.jl`

**Steps:**
1. Add `top_k_precision`, `average_precision`, and `jaccard_index` helpers.
2. Extend `summarize_case` with per-trait causal metrics:
   - `topk_precision`
   - `topk_true_positive_count`
   - `average_precision`
   - `mean_pip_active`
   - `mean_pip_inactive`
3. Extend multi-trait any-active summaries with:
   - `any_active_topk_precision`
   - `any_active_topk_true_positive_count`
   - `any_active_average_precision`
4. Extend `method_summary` aggregation with mean and standard deviation for the new metrics.
5. Add a focused production selector, `:bayesr_quality`, containing:
   - `MT_Annotated_BayesR`
   - `MT_Annotated_BayesC_I`
   - `Annotated_BayesR_y1`
   - `Annotated_BayesR_y2`
   - `BayesR_y1`
   - `BayesR_y2`

### Task 4: Implement Stability Summaries

**Files:**
- Modify: `benchmarks/simulated_annotations_multitrait_comparison.jl`

**Steps:**
1. Add `summarize_marker_stability`.
2. For each method/family and trait with at least two seeds, compute all pairwise:
   - PIP correlation
   - marker-effect estimate correlation
   - top-k Jaccard overlap
3. For two-trait families, also compute any-active PIP stability from `max(pip_y1, pip_y2)`.
4. Save:
   - `marker_stability_pairwise.csv`
   - `marker_stability_summary.csv`

### Task 5: Verify Code

**Commands:**
- `julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'`
- `julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'`

### Task 6: Run Production-Style Benchmark

**Baseline comparison command:**

```bash
JWAS_SIMULATED_MT_FOCUS_MODE=bayesr_quality \
JWAS_SIMULATED_MT_SEEDS=101,202,303,404,505 \
JWAS_SIMULATED_MT_CHAIN_LENGTH=10000 \
JWAS_SIMULATED_MT_BURNIN=2000 \
JWAS_SIMULATED_MT_OUTPUT_FREQ=20 \
JWAS_SIMULATED_MT_WARMUP=true \
julia --project=. --startup-file=no \
  benchmarks/simulated_annotations_multitrait_comparison.jl \
  benchmarks/out/mt_annotated_bayesr_quality_20260625
```

**Long-chain confirmatory command for the new method:**

```bash
JWAS_SIMULATED_MT_FOCUS_MODE=mt_annotated_bayesr \
JWAS_SIMULATED_MT_SEEDS=101,202,303,404,505 \
JWAS_SIMULATED_MT_CHAIN_LENGTH=50000 \
JWAS_SIMULATED_MT_BURNIN=10000 \
JWAS_SIMULATED_MT_OUTPUT_FREQ=50 \
JWAS_SIMULATED_MT_WARMUP=true \
julia --project=. --startup-file=no \
  benchmarks/simulated_annotations_multitrait_comparison.jl \
  benchmarks/out/mt_annotated_bayesr_quality_long_20260625
```

### Task 7: Report

**Files:**
- Create: `benchmarks/reports/2026-06-25-multitrait-annotated-bayesr-production-quality-report.md`
- Modify: `docs/plans/2026-06-24-multitrait-annotated-bayesr-implementation.md`

**Report contents:**
1. Exact benchmark command and dataset.
2. Method table.
3. Causal-variant discovery metrics:
   - average precision
   - top-k recall/precision
   - PIP gap
   - any-active and shared-effect metrics
4. Cross-seed stability metrics.
5. Limitations and recommendation.

### Task 8: Final Verification and Commit

**Commands:**
- `julia --project=docs --startup-file=no docs/make.jl`
- `git status --short`
- `git add ...`
- `git commit -m "benchmarks: add production BayesR quality simulation"`
