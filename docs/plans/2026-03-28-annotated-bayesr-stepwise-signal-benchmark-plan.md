# Annotated BayesR Stepwise-Signal Benchmark Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a third annotated BayesR production benchmark scenario with explicit step-2 and step-3 annotation signal, run it, and update the report.

**Architecture:** Extend the existing benchmark script with a new `stepwise_annotation_signal` scenario whose truth is defined in the sequential conditional probabilities rather than only joint class probabilities. Reuse the current production benchmark pipeline, smoke test, and reporting outputs so the new scenario is directly comparable with the first two.

**Tech Stack:** Julia, DataFrames.jl, CSV.jl, JWAS production `runMCMC`

---

### Task 1: Add the failing smoke-test expectation

**Files:**
- Modify: `test/unit/test_bayesr_parity.jl`
- Test: `test/unit/test_bayesr_parity.jl`

**Step 1: Write the failing test**

Extend the annotated BayesR benchmark subprocess test to run with:

- `JWAS_ANNOT_BENCH_SCENARIO=stepwise_annotation_signal`

and assert:

- `comparison_summary.csv` contains `scenario`
- all summary rows have `scenario == "stepwise_annotation_signal"`
- `truth_metadata.csv` contains `scenario`
- all truth rows have `scenario == "stepwise_annotation_signal"`

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected:

- the annotated benchmark subprocess test fails because the benchmark script does
  not yet support the new scenario name

**Step 3: Commit**

Do not commit yet.

### Task 2: Implement the new scenario in the benchmark script

**Files:**
- Modify: `benchmarks/annotated_bayesr_comparison.jl`
- Test: `test/unit/test_bayesr_parity.jl`

**Step 1: Add scenario truth support**

Update `scenario_class_probabilities` so it supports:

- `stepwise_annotation_signal`

Implement it by converting:

- baseline: `p1=0.10`, `p2=0.20`, `p3=0.20`
- enriched: `p1=0.30`, `p2=0.60`, `p3=0.60`

into joint class probabilities before sampling marker classes.

**Step 2: Keep scenario metadata threaded through outputs**

Ensure the new scenario label is carried into:

- truth metadata
- run rows
- summary rows
- PIP group rows
- annotation coefficient rows

**Step 3: Run the smoke test again**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected:

- the subprocess benchmark test passes

### Task 3: Run the full production benchmark

**Files:**
- Modify: none
- Output: `/tmp/annotated_bayesr_benchmark_20260328_stepwise_signal`

**Step 1: Run the benchmark**

Run:

```bash
JWAS_ANNOT_BENCH_SCENARIO=stepwise_annotation_signal \
  julia --project=. --startup-file=no \
  benchmarks/annotated_bayesr_comparison.jl \
  /tmp/annotated_bayesr_benchmark_20260328_stepwise_signal
```

**Step 2: Inspect outputs**

Collect:

- truth class counts
- method summary table
- true-class PIP ordering
- annotated BayesR annotation-coefficient summary

Expected:

- all three methods complete
- stepwise signal appears in truth metadata

### Task 4: Update the report

**Files:**
- Modify: `benchmarks/reports/2026-03-28-annotated-bayesr-benchmark-report.md`

**Step 1: Extend the report**

Add the third scenario:

- truth summary
- method-level summary table
- step-specific interpretation

Explicitly compare:

- sparse scenario
- less-sparse scenario
- stepwise-annotation-signal scenario

**Step 2: State the conclusion clearly**

Document whether:

- annotated BayesR step-2 and step-3 coefficients became identifiable
- EBV correlation improved or not
- the main weakness still appears to be model structure vs benchmark design

### Task 5: Add the implementation note and verify everything

**Files:**
- Create: `docs/plans/2026-03-28-annotated-bayesr-stepwise-signal-benchmark-implementation.md`
- Test: `test/unit/test_bayesr_parity.jl`
- Test: `test/runtests.jl`
- Docs: `docs/make.jl`

**Step 1: Write the implementation note**

Summarize:

- code changes
- benchmark command
- main results
- verification

**Step 2: Run full verification**

Run:

```bash
julia --project=. --startup-file=no test/runtests.jl
julia --project=docs --startup-file=no docs/make.jl
```

Expected:

- full suite passes
- docs build passes

**Step 3: Commit**

```bash
git add benchmarks/annotated_bayesr_comparison.jl \
        benchmarks/reports/2026-03-28-annotated-bayesr-benchmark-report.md \
        test/unit/test_bayesr_parity.jl \
        docs/plans/2026-03-28-annotated-bayesr-stepwise-signal-benchmark-design.md \
        docs/plans/2026-03-28-annotated-bayesr-stepwise-signal-benchmark-plan.md \
        docs/plans/2026-03-28-annotated-bayesr-stepwise-signal-benchmark-implementation.md
git commit -m "benchmarks: add stepwise annotated BayesR scenario"
```
