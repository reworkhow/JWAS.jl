# Multi-Trait Annotated BayesC Sampler Benchmark Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

## Goal

Extend the production benchmark harness to compare multi-trait annotated BayesC sampler I versus sampler II on the packaged `simulated_annotations` fixture, alongside conventional multi-trait BayesC and the requested single-trait BayesC/BayesR baselines.

## Architecture

Reuse the existing benchmark harness in `benchmarks/simulated_annotations_multitrait_comparison.jl`, expand the method matrix, keep all marker/truth comparisons keyed by normalized marker IDs, and compute:

- true `P11`-based shared recovery for multi-trait methods
- overlap-based shared recovery for paired single-trait methods
- genomic prediction and any-active recovery summaries

## Tech Stack

- Julia
- existing JWAS production path
- CSV/DataFrames benchmark summaries
- regression coverage in `test/unit/test_misc_coverage.jl`

---

### Task 1: Expand the benchmark method matrix

**Files**

- Modify: `benchmarks/simulated_annotations_multitrait_comparison.jl`

**Work**

- Add explicit cases for:
  - `MT_BayesC`
  - `MT_Annotated_BayesC_I`
  - `MT_Annotated_BayesC_II`
  - `BayesC_y1`, `BayesC_y2`
  - `Annotated_BayesC_y1`, `Annotated_BayesC_y2`
  - `BayesR_y1`, `BayesR_y2`
  - `Annotated_BayesR_y1`, `Annotated_BayesR_y2`
- Pass `multi_trait_sampler=:I` and `:II` into the annotated multi-trait BayesC cases.
- Preserve the production JWAS run path in `run_case`.

### Task 2: Add regression coverage for the benchmark contract

**Files**

- Modify: `test/unit/test_misc_coverage.jl`

**Work**

- Add a lightweight source-contract test that checks the benchmark file for:
  - the new method variants
  - the two annotated multi-trait sampler settings
  - the new shared-summary functions
  - removal of the old `min(PIP_y1, PIP_y2)` single-trait proxy

### Task 3: Implement the new shared-recovery summaries

**Files**

- Modify: `benchmarks/simulated_annotations_multitrait_comparison.jl`

**Work**

- Keep the existing joined-marker and explicit-key alignment checks.
- For multi-trait methods:
  - reconstruct `P11` from saved marker-effect samples
  - rank markers by `P11`
  - declare top `k_shared_true`
  - write per-seed and summarized shared recovery tables
- For single-trait families:
  - pair the trait-1 and trait-2 outputs by normalized marker ID
  - detect top `k_y1_true` and top `k_y2_true`
  - declare shared markers from the overlap of detected sets
  - write per-seed and summarized shared recovery tables
- Combine the multi-trait and single-trait shared summaries into one benchmark-level pleiotropy summary.

### Task 4: Run smoke and production benchmarks

**Files**

- Modify or create under `benchmarks/reports/`

**Work**

- Run a smoke matrix with short chains and warmup disabled to validate the full method set.
- Run a production benchmark with longer chains and multiple seeds.
- Save all generated outputs under a dedicated benchmark output directory.

### Task 5: Write implementation and report notes

**Files**

- Create: `docs/plans/2026-04-10-multitrait-annotated-bayesc-sampler-benchmark-implementation.md`
- Create: `benchmarks/reports/2026-04-10-multitrait-annotated-bayesc-sampler-benchmark-report.md`

**Work**

- Record what changed in the benchmark harness.
- Summarize the production results, with emphasis on:
  - sampler I versus sampler II
  - shared recovery versus single-trait overlap baselines
  - prediction accuracy for genomic prediction

### Task 6: Verification

**Commands**

Run, at minimum:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
julia --project=. --startup-file=no test/runtests.jl
```

If report/docs are touched further, also run any relevant build or smoke verification needed by those changes.
