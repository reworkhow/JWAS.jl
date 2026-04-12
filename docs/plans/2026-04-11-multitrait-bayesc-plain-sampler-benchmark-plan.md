# Multi-Trait BayesC Plain vs Annotated Sampler Benchmark Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extend the production sampler benchmark so it compares explicit sampler-I versus sampler-II behavior in both plain and annotated multi-trait BayesC on the packaged two-trait fixture.

**Architecture:** Reuse the existing benchmark harness, keep the packaged `simulated_annotations` fixture, add explicit plain multi-trait BayesC `:I` and `:II` method cases, retain keyed `P11` summaries, and add a focused four-method sampler-comparison mode. Pair that production rerun with a plain multi-trait tiny-case exact regression so both the plain and annotated paths have the same low-level backstop.

**Tech Stack:** Julia, JWAS production MCMC path, CSV/DataFrames benchmark summaries, unit coverage in `test/unit`.

---

### Task 1: Add a plain multi-trait exact-posterior regression

**Files:**
- Modify: `test/unit/test_multitrait_mcmc.jl`

**Step 1: Reuse the exact helper for a global-`Pi` prior**

- Keep the existing exact posterior and empirical sampler helpers.
- Add a second regression test that uses:
  - `JWAS.GlobalPiPrior`
  - the same one-marker two-trait setup
  - the same state labels `00`, `10`, `01`, `11`

**Step 2: Assert sampler-I and sampler-II agreement for plain multi-trait BayesC**

- Verify:
  - sampler I matches the exact posterior within Monte Carlo tolerance
  - sampler II matches the exact posterior within Monte Carlo tolerance
  - sampler I and sampler II are close to each other

**Step 3: Verify the focused multi-trait test file**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
```

Expected:
- the new plain multi-trait regression passes
- the existing annotated regression still passes

### Task 2: Extend the benchmark method matrix with explicit plain sampler cases

**Files:**
- Modify: `benchmarks/simulated_annotations_multitrait_comparison.jl`
- Modify: `test/unit/test_misc_coverage.jl`

**Step 1: Add explicit plain multi-trait cases**

- Extend `method_cases()` with:
  - `MT_BayesC_I`
  - `MT_BayesC_II`

Both should be:

- `method = "BayesC"`
- `annotated = false`
- `multitrait = true`
- `trait = ""`
- `multi_trait_sampler = :I` or `:II`

**Step 2: Keep backward compatibility**

- Keep the older `MT_BayesC` auto row if it is still useful for the broader
  report.
- Do not break the existing full benchmark unless the user explicitly wants
  the old auto row removed later.

**Step 3: Update contract coverage**

- Extend `test_misc_coverage.jl` so it checks:
  - the new plain `:I` and `:II` method cases exist
  - only the intended multi-trait cases carry non-`auto` sampler settings

### Task 3: Add a focused four-method sampler-comparison selector

**Files:**
- Modify: `benchmarks/simulated_annotations_multitrait_comparison.jl`
- Modify: `test/unit/test_misc_coverage.jl`

**Step 1: Add an explicit selection path**

- Add a focused selector for:
  - `MT_BayesC_I`
  - `MT_BayesC_II`
  - `MT_Annotated_BayesC_I`
  - `MT_Annotated_BayesC_II`

This can be:

- a new environment flag
- or an extension of the current focused selector

Use the smallest change that keeps the harness readable.

**Step 2: Preserve keyed summary behavior**

- Do not change any marker/truth join logic.
- Reuse the existing `P11` reconstruction, per-seed shared summaries, and
  mixing summary path.

**Step 3: Update contract coverage for the selector**

- Check that the focused sampler comparison returns the intended four variants
  only.

### Task 4: Run smoke and production sampler-comparison reruns

**Files:**
- Write outputs under a new `/tmp/...` directory

**Step 1: Smoke rerun**

Run a short focused benchmark for the four-method sampler comparison to verify:

- method selection
- plain and annotated sampler metadata
- summary file generation

**Step 2: Production rerun**

Run a longer-chain focused comparison for:

- `MT_BayesC_I`
- `MT_BayesC_II`
- `MT_Annotated_BayesC_I`
- `MT_Annotated_BayesC_II`

Use multiple seeds and enough chain length to reduce short-chain noise.

### Task 5: Write the follow-up report

**Files:**
- Create: `benchmarks/reports/2026-04-11-multitrait-bayesc-plain-sampler-benchmark-report.md`
- Create: `docs/plans/2026-04-11-multitrait-bayesc-plain-sampler-benchmark-implementation.md`

**Step 1: Summarize the exact-versus-production distinction**

- State that both the plain and annotated tiny-case regressions support the
  same target posterior for sampler I and sampler II.
- Frame the production rerun as a finite-chain mixing comparison.

**Step 2: Summarize the production rerun**

- Compare plain `:I` versus plain `:II`
- compare annotated `:I` versus annotated `:II`
- compare plain versus annotated runtime cost and shared-state recovery

**Step 3: Answer the main method question explicitly**

- Say whether the sampler-I versus sampler-II gap:
  - already appears in plain multi-trait BayesC
  - gets stronger under annotation
  - or is mostly annotation-specific

### Task 6: Verification

**Commands:**

Run, at minimum:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
julia --project=. --startup-file=no test/runtests.jl
```

If any docs under `docs/src` are changed, also run:

```bash
julia --project=docs --startup-file=no docs/make.jl
```
