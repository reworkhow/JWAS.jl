# Annotated BayesR Benchmark Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a production-path benchmark that compares dense BayesR, dense annotated BayesR, and dense annotated BayesC on a synthetic dataset with informative annotations, then save a benchmark report.

**Architecture:** Create a new dedicated benchmark script that generates synthetic production inputs, runs the three methods through `runMCMC`, and writes truth-aware summary CSVs. Add one subprocess smoke test first, watch it fail, then implement the script, run the long benchmark, and write the report.

**Tech Stack:** Julia, JWAS, CSV.jl, DataFrames.jl, Statistics stdlib

---

### Task 1: Add the failing benchmark smoke test

**Files:**
- Modify: `test/unit/test_bayesr_parity.jl`
- Test: `test/unit/test_bayesr_parity.jl`

**Step 1: Write the failing test**

Add a new subprocess testset that runs:

```julia
benchmark_script = joinpath(repo_root, "benchmarks", "annotated_bayesr_comparison.jl")
```

with small environment overrides:

- `JWAS_ANNOT_BENCH_CHAIN_LENGTH=60`
- `JWAS_ANNOT_BENCH_BURNIN=20`
- `JWAS_ANNOT_BENCH_N_OBS=30`
- `JWAS_ANNOT_BENCH_N_MARKERS=40`
- `JWAS_ANNOT_BENCH_SEEDS=2026,2027`

and asserts that the script writes:

- `comparison_runs.csv`
- `comparison_summary.csv`
- `pip_group_summary.csv`
- `annotation_coefficients.csv`
- `truth_metadata.csv`

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected: FAIL because `benchmarks/annotated_bayesr_comparison.jl` does not exist yet.

### Task 2: Implement the production benchmark script

**Files:**
- Create: `benchmarks/annotated_bayesr_comparison.jl`
- Reuse: `benchmarks/bayesr_parity_common.jl`
- Test: `test/unit/test_bayesr_parity.jl`

**Step 1: Implement synthetic data generation**

Add helpers to:

- parse environment overrides
- simulate genotype matrix
- simulate two annotation columns
- generate true BayesR classes using annotation-driven enrichment
- sample marker effects from the truth classes
- build phenotype with moderate heritability
- write production CSV inputs for `get_genotypes`

**Step 2: Implement method runners**

Add production runners for:

- `BayesR`
- annotated `BayesR`
- annotated `BayesC`

Each runner should:

- call `get_genotypes(...)`
- build the model
- call `runMCMC(...)`
- collect EBV correlation and posterior summaries

**Step 3: Implement summary writers**

Write:

- `comparison_runs.csv`
- `comparison_summary.csv`
- `pip_group_summary.csv`
- `annotation_coefficients.csv`
- `truth_metadata.csv`

Keep the format compact and stable enough for repeated reruns.

**Step 4: Run the smoke test to verify it passes**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected: the new subprocess test passes.

### Task 3: Run the full benchmark and write the report

**Files:**
- Create: `benchmarks/reports/2026-03-28-annotated-bayesr-benchmark-report.md`
- Create: `docs/plans/2026-03-28-annotated-bayesr-benchmark-implementation.md`
- Output: benchmark directory under `/tmp` or `benchmarks/out`

**Step 1: Run the full benchmark**

Run:

```bash
julia --project=. --startup-file=no benchmarks/annotated_bayesr_comparison.jl <outdir>
```

with defaults:

- `n_obs = 200`
- `n_markers = 1000`
- `chain_length = 10000`
- `burnin = 2000`
- seeds `2026,2027,2028,2029,2030`

**Step 2: Inspect output CSVs**

Read the written CSVs and compute the report narrative around:

- EBV correlation by method
- PIP separation for causal vs null SNPs
- PIP separation for enriched vs non-enriched SNPs
- annotation-coefficient direction and magnitude
- BayesR vs annotated BayesR vs annotated BayesC tradeoffs

**Step 3: Write the benchmark report**

Save the report to:

- `benchmarks/reports/2026-03-28-annotated-bayesr-benchmark-report.md`

and include:

- benchmark protocol
- data-generating assumptions
- method comparison tables
- key conclusions
- limitations

**Step 4: Write the implementation note**

Save:

- `docs/plans/2026-03-28-annotated-bayesr-benchmark-implementation.md`

and record:

- files changed
- commands run
- where outputs were written
- key findings

### Task 4: Final verification

**Files:**
- Verify: `test/unit/test_bayesr_parity.jl`
- Verify: `test/runtests.jl`

**Step 1: Re-run the smoke test**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected: PASS.

**Step 2: Run the full test suite**

Run:

```bash
julia --project=. --startup-file=no test/runtests.jl
```

Expected: PASS.

**Step 3: Commit**

Commit the benchmark script, test update, report, and implementation note with a benchmark-focused message.
