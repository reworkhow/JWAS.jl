# Annotated BayesR Less-Sparse Benchmark Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extend the annotated BayesR production benchmark with a second, less-sparse upper-class scenario and update the report to compare both scenarios.

**Architecture:** Add a scenario selector to the existing production benchmark script, update the smoke test to require the new scenario metadata, run the second full benchmark, and revise the report so the first and second scenarios are interpreted together.

**Tech Stack:** Julia, JWAS, CSV.jl, DataFrames.jl, Statistics stdlib

---

### Task 1: Add the failing scenario-selection smoke test

**Files:**
- Modify: `test/unit/test_bayesr_parity.jl`
- Test: `test/unit/test_bayesr_parity.jl`

**Step 1: Write the failing test**

Update the annotated BayesR benchmark subprocess test to set:

- `JWAS_ANNOT_BENCH_SCENARIO=less_sparse_upper_classes`

and assert:

- `comparison_summary.csv` has a `scenario` column
- all rows in that column equal `less_sparse_upper_classes`

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected: FAIL because the current benchmark script does not emit scenario
metadata.

### Task 2: Implement the scenario selector in the benchmark script

**Files:**
- Modify: `benchmarks/annotated_bayesr_comparison.jl`
- Test: `test/unit/test_bayesr_parity.jl`

**Step 1: Add scenario configuration**

Implement a helper that maps scenario name to:

- baseline class probabilities
- enriched class probabilities

Support:

- `sparse_upper_classes`
- `less_sparse_upper_classes`

and error on unknown values.

**Step 2: Thread scenario metadata through outputs**

Add `scenario` to:

- `comparison_runs.csv`
- `comparison_summary.csv`
- `pip_group_summary.csv`
- `annotation_coefficients.csv`
- `truth_metadata.csv`

**Step 3: Re-run the smoke test**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected: PASS.

### Task 3: Run the second full benchmark

**Files:**
- Output: `/tmp/annotated_bayesr_benchmark_20260328_less_sparse`

**Step 1: Run the benchmark**

Run:

```bash
JWAS_ANNOT_BENCH_SCENARIO=less_sparse_upper_classes \
julia --project=. --startup-file=no benchmarks/annotated_bayesr_comparison.jl /tmp/annotated_bayesr_benchmark_20260328_less_sparse
```

using the same default chain settings as the first benchmark.

**Step 2: Inspect outputs**

Read:

- `comparison_summary.csv`
- `pip_group_summary.csv`
- `annotation_coefficients.csv`
- `truth_metadata.csv`

and compare them against the first scenario.

### Task 4: Update report and implementation note

**Files:**
- Modify: `benchmarks/reports/2026-03-28-annotated-bayesr-benchmark-report.md`
- Create: `docs/plans/2026-03-28-annotated-bayesr-less-sparse-benchmark-implementation.md`

**Step 1: Update the report**

Revise the report to include:

- scenario 1 results
- scenario 2 results
- direct comparison of what changed

**Step 2: Write the implementation note**

Record:

- files changed
- commands run
- output directories
- what the second scenario changed scientifically

### Task 5: Final verification

**Files:**
- Verify: `test/unit/test_bayesr_parity.jl`
- Verify: `test/runtests.jl`

**Step 1: Re-run the smoke test**

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

**Step 2: Re-run the full suite**

```bash
julia --project=. --startup-file=no test/runtests.jl
```

**Step 3: Commit**

Commit the scenario extension and updated report with a benchmark-focused
message.
