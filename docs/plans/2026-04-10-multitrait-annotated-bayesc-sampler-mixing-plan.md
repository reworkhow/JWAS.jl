# Multi-Trait Annotated BayesC Sampler Mixing Benchmark Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extend the multi-trait sampler benchmark so the sampler `:I` versus `:II` comparison uses longer production-path chains and explicit shared-state mixing diagnostics, with conventional multi-trait BayesC as the reference.

**Architecture:** Reuse the existing benchmark harness and output structure, add a focused rerun mode for the three multi-trait methods, preserve keyed marker/truth joins, and add per-seed plus summarized `P11`-based diagnostics. The earlier exact-posterior regression test remains the correctness backstop; this benchmark should measure mixing and efficiency.

**Tech Stack:** Julia, JWAS production MCMC path, CSV/DataFrames benchmark summaries, unit coverage in `test/unit`.

---

### Task 1: Record the new exact-posterior regression coverage

**Files:**
- Modify: `test/unit/test_multitrait_mcmc.jl`

**Step 1: Keep the exact tiny-case helper and sampler-equivalence test**

- Retain the one-marker two-trait helpers that compute the exact posterior and empirical sampler frequencies.
- Keep the regression asserting sampler I and sampler II both match the same target posterior within finite Monte Carlo tolerance.

**Step 2: Verify the focused multi-trait test file**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
```

Expected:
- all multi-trait tests pass
- the new sampler-equivalence test passes

### Task 2: Add a focused multi-trait rerun mode to the benchmark harness

**Files:**
- Modify: `benchmarks/simulated_annotations_multitrait_comparison.jl`

**Step 1: Add a method filter for the longer-chain rerun**

- Add a small, explicit selection path so the harness can run only:
  - `MT_BayesC`
  - `MT_Annotated_BayesC_I`
  - `MT_Annotated_BayesC_II`
- Keep the full matrix as the default behavior unless the focused mode is requested.

**Step 2: Preserve current keyed summary behavior**

- Do not change marker/truth alignment logic.
- Reuse the existing `P11` reconstruction and explicit marker-key joins.

### Task 3: Add explicit shared-state mixing summaries

**Files:**
- Modify: `benchmarks/simulated_annotations_multitrait_comparison.jl`

**Step 1: Extend per-seed multi-trait summaries**

- Ensure per-seed outputs already include or now include:
  - shared precision
  - shared recall
  - shared F1
  - mean `P11` on true shared loci
  - mean `P11` on nonshared loci
  - runtime

**Step 2: Add summarized variability across seeds**

- Produce a compact summary table for the focused multi-trait rerun that includes:
  - mean and standard deviation across seeds for shared precision/recall/F1
  - mean and standard deviation across seeds for true-shared `P11`
  - mean runtime

### Task 4: Run smoke and production focused reruns

**Files:**
- Write outputs under a new `/tmp/...` benchmark directory

**Step 1: Smoke focused rerun**

Run a short benchmark using only the three multi-trait methods to verify the focused mode and new summaries.

**Step 2: Production focused rerun**

Run a longer-chain comparison for:
- `MT_BayesC`
- `MT_Annotated_BayesC_I`
- `MT_Annotated_BayesC_II`

Use multiple seeds and enough chain length to reduce the chance that the
sampler-I versus sampler-II gap is just short-chain noise.

### Task 5: Write the follow-up report

**Files:**
- Create: `benchmarks/reports/2026-04-10-multitrait-annotated-bayesc-sampler-mixing-report.md`
- Create: `docs/plans/2026-04-10-multitrait-annotated-bayesc-sampler-mixing-implementation.md`

**Step 1: Summarize the exact-versus-mixing distinction**

- State clearly that the new exact tiny-case regression supports a common target posterior for both samplers.
- Frame the benchmark as a finite-chain mixing and runtime comparison.

**Step 2: Summarize the focused rerun results**

- Compare sampler I versus II on:
  - shared recovery
  - seed-to-seed variability
  - `P11` separation
  - runtime
- Use conventional multi-trait BayesC as the reference row.

### Task 6: Verification

**Commands:**

Run, at minimum:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
julia --project=. --startup-file=no test/runtests.jl
```

If the report or docs are updated in-tree, also run:

```bash
julia --project=docs --startup-file=no docs/make.jl
```
